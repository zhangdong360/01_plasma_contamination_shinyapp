library(ggplot2)
library(tidyr)
library(dplyr)

plot_protein_by_sample <- function(data) {
  # data 为 matrix，行为 protein, 列为 sample id
  
  # 将矩阵转换为数据框
  data <- as.data.frame(t(data))
  data$sample_id <- rownames(data)
  
  # 将数据转换为长格式
  data_long <- tidyr::gather(data, key = "protein", value = "value", -sample_id)
  # 绘制箱线图
  p <- ggplot(data_long, aes(x = sample_id, y = log2(value + 1), fill = sample_id)) + 
    geom_boxplot() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black", size = 10),
          axis.text.y = element_text(hjust = 1, colour = "black", size = 10),
          plot.title = element_text(hjust = 0.5),
          legend.position = "right") 
  return(p)
}


