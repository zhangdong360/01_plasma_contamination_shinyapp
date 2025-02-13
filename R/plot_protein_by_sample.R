library(ggplot2)
library(tidyr)
library(dplyr)

plot_protein_by_sample <- function(data, protein = NULL, group = NULL) {
  # data 为 matrix，行为 protein, 列为 sample id
  # protein 默认为 NULL，即对所有 protein 进行绘制
  # 如果想要绘制特定的 protein，则输入要绘制的 protein 向量
  # 如果要进行分组比较，则输入 group，默认不进行分组比较
  # group 要求为两列，id 为 sample id，group 为分组
  
  # 如果有指定的 protein，过滤数据
  if (!is.null(protein)) {
    data <- data[rownames(data)%in%protein,]
  }
  # 将矩阵转换为数据框
  data <- as.data.frame(t(data))
  data$sample_id <- rownames(data)
  # 如果有分组信息，合并分组信息
  if (!is.null(group)) {
    data <- merge(data, group, by.x = "sample_id", by.y = "id")
  }
  # 将数据转换为长格式
    data_long <- tidyr::gather(data, key = "protein", value = "value", -sample_id, -group)
  # 绘制箱线图
  if (is.null(group)) {
    p <- ggplot(data_long, aes(x = sample_id, y = log2(value + 1))) + 
      geom_boxplot() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black", size = 10),
            axis.text.y = element_text(hjust = 1, colour = "black", size = 10),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none")
  } else {
    p <- ggplot(data_long, aes(x = sample_id, y = log2(value + 1), fill = group)) + 
      geom_boxplot() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black", size = 10),
            axis.text.y = element_text(hjust = 1, colour = "black", size = 10),
            plot.title = element_text(hjust = 0.5),
            legend.position = "right") 
  }
  return(p)
}
