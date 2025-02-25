plot_pvalue_distribution <- function(data, type = "pearson", alpha = 0.05, 
                                     title = "P值分布直方图", xlab = "P值", ylab = "频数") {
  # 检查必要包安装
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("请先安装ggplot2包: install.packages('ggplot2')")
  }
  
 
  # 初始化p值矩阵
  n <- ncol(data)
  p_values <- matrix(NA, nrow = n, ncol = n)
  
  test_result <- Hmisc::rcorr(t(data),type = type)
  p_values <- test_result$P

  
  # 创建绘图数据
  p_df <- data.frame(p_value = p_values[upper.tri(p_values)])
  
  # 生成ggplot对象
  ggplot2::ggplot(p_df, ggplot2::aes(x = p_value)) +
    ggplot2::geom_histogram(
      ggplot2::aes(fill = p_value < alpha),
      binwidth = 0.05,
      color = "black",
      boundary = 0
    ) +
    ggplot2::scale_fill_manual(values = c("FALSE" = "#87CEEB", "TRUE" = "#FF6B6B"), guide = "none") +
    ggplot2::geom_vline(xintercept = alpha, color = "#FF5252", linetype = "dashed", linewidth = 0.8) +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = ylab,
      caption = sprintf("显著性阈值: α = %.2f (红色虚线)", alpha)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid.major.x = ggplot2::element_blank()
    )
}

# 示例用法
# data(mtcars)
# plot_pvalue_distribution(mtcars)
# plot_pvalue_distribution(iris, alpha = 0.01, title = "鸢尾花数据集P值分布")
