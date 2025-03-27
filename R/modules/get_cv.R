# get_Cv ----
# raw_data要求行为蛋白，列为样本
get_cv <- function(raw_data,protein = rownames(raw_data),
                   data_group = NULL,
                   by_group = F){
  raw_data <- raw_data[rownames(raw_data)%in%protein,,drop = F]
  result_cv <- matrix(nrow = dim(raw_data)[1])
  result_cv <- as.data.frame(result_cv)
  if(by_group){
    for (group in unique(data_group$group)){
      data_group_one <- data_group[data_group$group%in% group,]
      data <- raw_data[,colnames(raw_data)%in%data_group_one$id]
      mean_sd <- apply(t(data), 2, function(x) c(mean = mean(x), sd = sd(x)))
      mean_sd <- as.data.frame(mean_sd)
      mean_sd <- as.data.frame(t(mean_sd))
      # 计算每个样本的变异系数（CV）
      result_cv[,paste0(group,"_CV")] <- mean_sd$sd / mean_sd$mean
    }
    result_cv <- result_cv[,-1]
  }else if(by_group == F){
    mean_sd <- apply(t(raw_data), 2, function(x) c(mean = mean(x), sd = sd(x)))
    mean_sd <- as.data.frame(mean_sd)
    mean_sd <- as.data.frame(t(mean_sd))
    # 计算每个样本的变异系数（CV）
    result_cv <- mean_sd
    result_cv[,"CV"] <- result_cv$sd / result_cv$mean
  }
  return(result_cv)
}
