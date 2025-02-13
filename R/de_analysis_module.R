limma_proteomics_analysis <- function(expr_matrix, 
                                      group_matrix, 
                                      compare = NULL, 
                                      p_cutoff = 0.05, 
                                      logFC_cutoff = 1, 
                                      adjust_method = "BH", 
                                      apply_voom = FALSE,
                                      p_type = "adjusted") {
  # 检查依赖包
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' required but not installed.")
  }
  
  # 验证输入数据结构
  if (!all(c("id", "group") %in% colnames(group_matrix))) {
    stop("分组矩阵必须包含'id'和'group'两列")
  }
  # 校验p_type参数有效性
  if (!p_type %in% c("adjusted", "raw")) {
    stop("p_type参数必须为'adjusted'或'raw'")
  }
  # 确保样本顺序匹配
  sample_ids <- colnames(expr_matrix)
  if (!all(sample_ids %in% group_matrix$id)) {
    missing_samples <- setdiff(sample_ids, group_matrix$id)
    stop("以下样本缺少分组信息: ", paste(missing_samples, collapse = ", "))
  }
  
  # 提取与表达矩阵匹配的分组信息
  matched_groups <- group_matrix$group[match(sample_ids, group_matrix$id)]
  group <- factor(matched_groups)
  
  # 设置默认比较组
  if (is.null(compare)) {
    compare <- levels(group)[2:1] # 默认第二个组 vs 第一个组  
  } else if (!all(compare %in% levels(group))) {
    stop("对比组必须存在于分组信息中")
  }
  
  # 数据转换（根据需求选择）
  if (apply_voom) {
    expr_matrix <- limma::voom(expr_matrix)$E
  }
  
  # 构建线性模型
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  
  # 创建对比
  contrast_formula <- paste(compare[1], "-", compare[2])
  contrast_matrix <- limma::makeContrasts(
    contrasts = contrast_formula, 
    levels = design
  )
  
  # 模型拟合
  fit <- limma::lmFit(expr_matrix, design)
  fit <- limma::contrasts.fit(fit, contrast_matrix)
  fit <- limma::eBayes(fit)
  
  # 提取结果
  top_table <- limma::topTable(
    fit, 
    number = Inf, 
    adjust.method = adjust_method
  )
  # 动态选择p值列 --------------------------------------------------
  p_column <- ifelse(p_type == "adjusted", "adj.P.Val", "P.Value")
  
  # 添加显著性标记
  top_table$significant <- (abs(top_table$logFC) >= logFC_cutoff) & 
    (top_table[, p_column] <= p_cutoff)
  
  # 添加蛋白质ID列
  if (!"Protein" %in% colnames(top_table)) {
    top_table <- cbind(Protein = rownames(top_table), top_table)
    rownames(top_table) <- NULL
  }
  
  return(top_table)
}
