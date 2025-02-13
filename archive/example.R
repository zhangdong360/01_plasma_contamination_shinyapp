source("./R/data_check.R")
source("./R/data_correct.R")
source("./R/plot_protein_by_sample.R")
source("./R/plot_protein.R")
df <- read.csv("./tests/raw_data_aggr.csv")
rownames(df) <- df[, 1]  # 将第一列设置为行名
df <- df[, -1, drop = FALSE]  # 删除第一列
df
data_group <- read.csv("./tests/group.csv")
rownames(data_group) <- data_group[, 1]  # 将第一列设置为行名
result_check <- data_check(data = df,data_group = data_group,cutoff = 0.9)
result_check$data
result_correct <- data_correct(data = result_check,type = "all")

plot_protein_by_sample(data = df,protein = result_check$marker_list$erythrocyte,group = data_group)
plot_protein(data = df,protein = result_check$marker_list$erythrocyte)
result <- limma_proteomics_analysis(expr_matrix = log2(df),
                                    group_matrix = data_group,
                                    compare = c("AD","Ctrl"),
                                    p_type = "raw")
plot_1 <- create_volcano_plot(limma_result = result,
                              p_type = "raw",gene_col = "plot_1",
                              group_names = c("AD","Ctrl"),logFC_cutoff = 0)
plot_1
