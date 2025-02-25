#' 表达矩阵相关性分析及可视化
#'
#' @param exprMatrix 输入矩阵或相关系数矩阵
#' @param input_type 输入类型（"matrix"/"correlation"，默认"matrix"）
#' @param corMethod 相关系数计算方法（仅当input_type="matrix"时生效）
#' @param colors 颜色梯度（默认蓝白红渐变）
#' @param displayNumbers 是否显示相关系数值
#' @param title 图像标题
#' @param cluster_rows 是否对行聚类
#' @param cluster_cols 是否对列聚类
#' @return 返回相关系数矩阵和热图对象
#' @export
#'
#' @examples
#' # 示例1：输入数据矩阵
#' expr <- matrix(rnorm(1000), nrow=100, dimnames=list(paste0("Gene",1:100), paste0("Sample",1:10)))
#' plot_expression_correlation(expr, input_type = "matrix")
#'
#' # 示例2：直接输入相关系数矩阵
#' cor_matrix <- cor(expr)
#' plot_expression_correlation(cor_matrix, input_type = "correlation", title = "Pre-computed Correlation")

plot_expression_correlation <- function(exprMatrix,
                                        input_type = "matrix",
                                        corMethod = "pearson",
                                        colors = colorRampPalette(c("blue", "white", "red"))(100),
                                        displayNumbers = FALSE,
                                        title = "Expression Correlation Matrix",
                                        cluster_rows = TRUE,
                                        cluster_cols = TRUE) {
  
  # 检查必要包
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Please install the package first: install.packages('pheatmap')")
  }
  
  # 根据输入类型处理数据
  if (input_type == "matrix") {
    # 计算相关系数矩阵
    cor_matrix <- cor(exprMatrix, method = corMethod)
  } else if (input_type == "correlation") {
    # 直接使用输入作为相关系数矩阵
    cor_matrix <- exprMatrix
  } else {
    stop("Invalid input_type. Please use 'matrix' or 'correlation'.")
  }
  
  # 设置显示数值
  annotation <- if (displayNumbers) round(cor_matrix, 2) else FALSE
  
  # 绘制热图
  p <- pheatmap::pheatmap(
    mat = cor_matrix,
    color = colors,
    display_numbers = annotation,
    number_color = "black",
    main = title,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    border_color = NA,
    show_rownames = TRUE,
    show_colnames = TRUE
  )
  
  # 返回结果
  return(list(correlation_matrix = cor_matrix, plot = p))
}
