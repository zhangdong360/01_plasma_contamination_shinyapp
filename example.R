source("./data_check_v02.R")
source("./data_correct_v02.R")
df <- read.csv("./raw_data_aggr.csv")
rownames(df) <- df[, 1]  # 将第一列设置为行名
df <- df[, -1, drop = FALSE]  # 删除第一列
df
data_group <- read.csv("./group.csv")
rownames(data_group) <- data_group[, 1]  # 将第一列设置为行名
result_check <- data_check(data = df,data_group = data_group,cutoff = 0.9)
result_check$data
result_correct <- data_correct(data = result_check,type = "all")
