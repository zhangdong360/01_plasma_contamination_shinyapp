create_volcano_plot <- function(limma_result, 
                                p_type = "adjusted", 
                                p_cutoff = 0.05, 
                                logFC_cutoff = 1,
                                gene_col = "Protein",
                                group_names = c("Group1", "Group2"),
                                colors = c(Up = "#E64B35", Down = "#4DBBD5", Not = "grey80")) {
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
  
  # 创建基础火山图
  volc <- ggplot2::ggplot(data = plot_data, 
                          aes(x = logFC, 
                              y = -log10(.data[[p_column]]), 
                              color = sig)) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
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
  
  if (nrow(top_genes) > 0) {
    volc <- volc +
      ggrepel::geom_text_repel(
        data = top_genes,
        aes(label = gene),
        size = 3,
        box.padding = 0.5,
        max.overlaps = Inf,
        segment.color = "grey50"
      )
  }
  
  return(volc)
}
