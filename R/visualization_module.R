
create_volcano_plot <- function(limma_result, 
                                p_type = "adjusted", 
                                p_cutoff = 0.05, 
                                logFC_cutoff = 1,
                                gene_col = "Protein",
                                group_names = c("Group1", "Group2"),
                                colors = c(Up = "#E64B35", Down = "#4DBBD5", Not = "grey80"),
                                base_size = 14,        # 新增基础字号参数
                                label_size = 4.5,      # 新增标签字号参数
                                auto_font = TRUE) {    # 新增自动调节开关
  # 检查必要包
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required but not installed.")
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' required but not installed.")
  }
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  # 校验参数
  if (!p_type %in% c("adjusted", "raw")) {
    stop("p_type参数必须为 'adjusted' 或 'raw'")
  }
  if (!all(c("logFC", "P.Value", "adj.P.Val", gene_col) %in% colnames(limma_result))) {
    stop("输入数据缺少必要列")
  }
  
  # 选择使用的p值列
  p_column <- ifelse(p_type == "adjusted", "adj.P.Val", "P.Value")
  
  # 准备绘图数据
  plot_data <- limma_result
  plot_data$gene <- plot_data[[gene_col]]
  
  # 创建显著性分类
  plot_data$sig <- "Not"
  plot_data$sig[plot_data[[p_column]] < p_cutoff & plot_data$logFC >= logFC_cutoff] <- "Up"
  plot_data$sig[plot_data[[p_column]] < p_cutoff & plot_data$logFC <= -logFC_cutoff] <- "Down"
  plot_data$sig <- factor(plot_data$sig, levels = c("Up", "Down", "Not"))
  
  # 智能字号调节逻辑
  if(auto_font) {
    n_points <- nrow(limma_result)
    base_size <- dplyr::case_when(
      n_points > 10000 ~ 12,
      n_points > 5000 ~ 14,
      TRUE ~ 16
    )
    label_size <- dplyr::case_when(
      n_points > 10000 ~ 3.5,
      n_points > 5000 ~ 4,
      TRUE ~ 4.5
    )
  }
  # 创建基础火山图
  volc <- ggplot2::ggplot(data = plot_data, 
                          aes(x = logFC, 
                              y = -log10(.data[[p_column]]), 
                              color = sig)) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::theme_classic(base_size = base_size) +  # 应用基础字号
    ggplot2::theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size * 1.2),
      axis.title = element_text(size = base_size * 1.1),
      axis.text = element_text(size = base_size * 0.9),
      legend.text = element_text(size = base_size * 0.9),
      legend.title = element_text(size = base_size)
    ) +
    ggplot2::scale_color_manual(
      values = colors,
      name = "Expression",
      labels = c(
        Up = paste0("Up (n=", sum(plot_data$sig == "Up"), ")"),
        Down = paste0("Down (n=", sum(plot_data$sig == "Down"), ")"),
        Not = "Not significant"
      )
    ) +
    ggplot2::labs(
      x = expression("log"[2]*" Fold Change"),
      y = ifelse(p_type == "adjusted",
                 expression("-log"[10]*" (Adj.P)"),
                 expression("-log"[10]*" (P)")),
      title = paste0(group_names[1], " vs ", group_names[2])
    )
  
  # 添加显著性阈值线
  volc <- volc +
    ggplot2::geom_hline(
      yintercept = -log10(p_cutoff),
      linetype = "dashed",
      color = "black",
      linewidth = 0.6,
      alpha = 0.8
    ) +
    ggplot2::geom_vline(
      xintercept = c(-logFC_cutoff, logFC_cutoff),
      linetype = "dashed",
      color = "black",
      linewidth = 0.6,
      alpha = 0.8
    )
  
  # 添加标签（显示前10个显著基因）
  top_genes <- plot_data %>%
    dplyr::filter(sig %in% c("Up", "Down")) %>%
    dplyr::arrange(desc(abs(logFC))) %>%
    head(10)
  
  # 添加标签（动态调节标签字号）
  if(nrow(top_genes) > 0) {
    volc <- volc +
      ggrepel::geom_text_repel(
        data = top_genes,
        aes(label = gene),
        size = label_size,          # 应用标签字号
        box.padding = 0.5,
        max.overlaps = Inf,
        segment.color = "grey50"
      )
  }
  return(volc)
}
venn_plot <- function(set1, set2, 
                      categories = c("Set A", "Set B"),
                      title = "Venn Diagram",
                      colors = c("#1b9e77", "#d95f02"),
                      alpha = 0.5,
                      print.mode = c( "raw","percent"),
                      save.plot = FALSE,
                      filename = "venn_diagram.png") {
  
  # 检查并加载必要包
  if (!require(VennDiagram)) {
    install.packages("VennDiagram")
    library(VennDiagram)
  }
  
  # 输入处理
  if (is.vector(set1) && is.vector(set2)) {
    # 将向量转换为集合列表
    universe <- unique(c(set1, set2))
    set_list <- list(
      set1 = unique(set1),
      set2 = unique(set2)
    )
  } else if (is.list(set1) && length(set1) == 2) {
    # 如果输入是列表的情况
    set_list <- set1
    universe <- unique(unlist(set1))
  } else {
    stop("输入应为两个向量或包含两个集合的列表")
  }
  
  # 创建Venn图对象
  venn <- VennDiagram::venn.diagram(
    x = set_list,
    category.names = categories,
    filename = NULL,  # 不直接保存文件
    output = TRUE,
    
    # 可视化参数
    fill = colors,
    alpha = alpha,
    lty = "blank",
    cex = 1.5,
    cat.cex = 1.3,
    cat.pos = c(-30, 30),
    cat.dist = c(0.05, 0.05),
    margin = 0.05,
    
    # 数值显示设置
    print.mode = print.mode,
    sigdigs = 2
  )
  
  # 绘制图形
  grid::grid.newpage()
  grid::grid.draw(venn)
  grid::grid.text(title, y = 0.95, gp = grid::gpar(fontsize = 16))
  
  # 可选保存功能
  if (save.plot) {
    ggplot2::ggsave(filename, plot = venn, 
                    width = 8, height = 6, 
                    dpi = 300)
    message("图形已保存为：", filename)
  }
  
  # 返回交集信息
  overlap <- intersect(set_list[[1]], set_list[[2]])
  return(invisible(list(
    total_unique = length(universe),
    overlap_count = length(overlap),
    overlap_elements = overlap
  )))
}
