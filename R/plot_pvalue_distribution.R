plot_stat_distribution <- function(data, type = "pearson", statistic = "pvalue", alpha = 0.05,
                                   title = NULL, xlab = NULL, ylab = "频数",
                                   binwidth = NULL, boundary = NULL,
                                   add_fit = TRUE, fit_color = "#2E8B57", fit_size = 1.2) {
  # 检查必要包安装
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("请先安装ggplot2包: install.packages('ggplot2')")
  }
  if (!requireNamespace("Hmisc", quietly = TRUE)) {
    stop("请先安装Hmisc包: install.packages('Hmisc')")
  }
  
  # 设置默认参数
  if (is.null(title)) {
    title <- if (statistic == "pvalue") "P值分布直方图" else "相关系数分布直方图"
  }
  if (is.null(xlab)) {
    xlab <- if (statistic == "pvalue") "P值" else "相关系数"
  }
  if (is.null(binwidth)) {
    binwidth <- if (statistic == "pvalue") 0.05 else 0.1
  }
  if (is.null(boundary)) {
    boundary <- if (statistic == "pvalue") 0 else -1
  }
  
  # 计算相关矩阵和p值矩阵
  test_result <- Hmisc::rcorr(as.matrix(t(data)), type = type)
  values <- if (statistic == "pvalue") test_result$P else test_result$r
  
  # 创建绘图数据
  df <- data.frame(value = values[upper.tri(values)])
  
  # 创建基础绘图对象
  p <- ggplot2::ggplot(df, ggplot2::aes(x = value))
  
  # 根据统计量类型添加图层
  if (statistic == "pvalue") {
    p <- p +
      ggplot2::geom_histogram(
        ggplot2::aes(fill = value < alpha, y = after_stat(count)),
        binwidth = binwidth,
        color = "black",
        boundary = boundary
      ) +
      ggplot2::scale_fill_manual(
        values = c("FALSE" = "#87CEEB", "TRUE" = "#FF6B6B"), 
        guide = "none"
      ) +
      ggplot2::geom_vline(
        xintercept = alpha, 
        color = "#FF5252", 
        linetype = "dashed", 
        linewidth = 0.8
      )
    
    # 添加P值的核密度估计
    if (add_fit) {
      p <- p + 
        ggplot2::geom_line(
          ggplot2::aes(y = after_stat(count * binwidth * density)),
          stat = "density",
          color = fit_color,
          linewidth = fit_size
        )
    }
  } else {
    p <- p +
      ggplot2::geom_histogram(
        ggplot2::aes(y = after_stat(count)),
        fill = "#87CEEB",
        binwidth = binwidth,
        color = "black",
        boundary = boundary
      )
    
    # 添加相关系数的密度曲线
    if (add_fit) {
      p <- p + 
        ggplot2::geom_line(
          ggplot2::aes(y = after_stat(count * binwidth * density)),
          stat = "density",
          color = fit_color,
          linewidth = fit_size
        )
    }
  }
  
  # 添加通用元素
  caption_text <- if (statistic == "pvalue") {
    sprintf("显著性阈值: α = %.2f (红色虚线)\n绿色曲线：核密度估计", alpha)
  } else {
    "绿色曲线：核密度估计"
  }
  
  p <- p +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = ylab,
      caption = caption_text
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(hjust = 0)
    )
  
  return(p)
}
