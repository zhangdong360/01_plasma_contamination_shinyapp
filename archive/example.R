source("./R/data_check.R")
source("./R/data_correct.R")
source("./R/plot_protein_by_sample.R")
source("./R/plot_protein.R")
source("./R/plot_expression_correlation.R")
source("./R/modules/get_cv.R")
library(ggsci)
# test
df <- read.csv("./tests/raw_data_aggr.csv")
rownames(df) <- df[, 1]  # 将第一列设置为行名
df <- df[, -1, drop = FALSE]  # 删除第一列
df
data_group <- read.csv("./tests/group.csv")
rownames(data_group) <- data_group[, 1]  # 将第一列设置为行名
# cov
# df <- read.csv("../01_plasma_contamination/01_rawdata/data_COVID_19/Plasma_matrix_after_seqknn.csv")
# colnames(df)[1] <- "From"
# library(readr)
# gene_name <- read_delim("../01_plasma_contamination/01_rawdata/data_COVID_19/idmapping_2025_02_26.tsv", 
#                            delim = "\t", escape_double = FALSE, 
#                            trim_ws = TRUE)
# gene_name <- as.data.frame(gene_name)
# df <- merge(gene_name,df, by = "From")
# rownames(df) <- df$To  # 将第一列设置为行名
# df <- subset(df,select = -c(From,To))  # 删除第一列
# df
# write.csv(df,file = "../01_plasma_contamination/01_rawdata/data_COVID_19/data_matrix_input.csv")

df <- read.csv("../01_plasma_contamination/01_rawdata/data_COVID_19/data_matrix_input.csv")
rownames(df) <- df[, 1]  # 将第一列设置为行名
df <- df[, -1, drop = FALSE]  # 删除第一列
df
data_group <- read.csv("../01_plasma_contamination/01_rawdata/data_COVID_19/group.csv")
rownames(data_group) <- data_group[, 1]  # 将第一列设置为行名

result_check <- data_check(data = df,data_group = data_group,cutoff = 0.5)

result_correct <- data_correct(data = result_check,
                               type = c("coagulation","erythrocyte","platelet"),
                               erythrocyte_marker = result_check$marker_list$erythrocyte,
                               coagulation_marker = result_check$marker_list$coagulation,
                               platelet_marker = result_check$marker_list$platelet)
# CV ----
result_cv_coa <- get_cv(raw_data = result_check$rawdata,
                        protein = result_check$marker_list$coagulation)
result_cv_ery <- get_cv(raw_data = result_check$rawdata,
                        protein = result_check$marker_list$erythrocyte)
result_cv_pla <- get_cv(raw_data = result_check$rawdata,
                        protein = result_check$marker_list$platelet)
result_cv_other <- get_cv(raw_data =  result_check$rawdata[!rownames(result_check$rawdata)%in%
                                                     c(result_check$marker_list$coagulation,
                                                       result_check$marker_list$erythrocyte,
                                                       result_check$marker_list$platelet),])
result_cv_coa$type <- "coagulation"
result_cv_ery$type <- "erythrocyte"
result_cv_pla$type <- "platelet"
result_cv_other$type <- "other protein"
result_cv <- rbind(result_cv_coa,result_cv_ery,
                   result_cv_pla,result_cv_other)
result_cv$type <- factor(result_cv$type,
                         levels = c("coagulation",
                                    "erythrocyte",
                                    "platelet",
                                    "other protein"))

ggplot(result_cv,aes(x = type , y = CV, fill = type)) +
  geom_violin() +
  geom_boxplot(fill = "white",width = 0.2) +
  scale_fill_manual(values = c("coagulation" = "#E64B35FF",
                               "erythrocyte" = "#F39B7FFF",
                               "platelet" = "#7E6148FF",
                               "other protein" = "#3C5488FF")) +
  stat_compare_means(comparisons = list(c("other protein","coagulation"),
                                        c("other protein","erythrocyte"),
                                        c("other protein","platelet"))) + 
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13))
# cor ----


result_check$data
result_cor <- test_result <- Hmisc::rcorr(as.matrix(t(result_check$rawdata)), type = "pearson")
result_cor <- as.data.frame(result_cor$r)

result_cor_erythrocyte <- abs(result_cor[,colnames(result_cor)%in%c(result_check$marker_list$erythrocyte)])
result_cor_coagulation <- abs(result_cor[,colnames(result_cor)%in%c(result_check$marker_list$coagulation)])
result_cor_platelet <- abs(result_cor[,colnames(result_cor)%in%c(result_check$marker_list$platelet)])
result_cor_erythrocyte$avg <- NA
result_cor_erythrocyte <- as.matrix(result_cor_erythrocyte)
for (i in 1:dim(result_cor_erythrocyte)[1]) {
  result_cor_erythrocyte[i,"avg"] <- mean(result_cor_erythrocyte[i,1:dim(result_cor_erythrocyte)[2]-1])
}

result_cor_coagulation$avg <- NA
result_cor_coagulation <- as.matrix(result_cor_coagulation)
for (i in 1:dim(result_cor_coagulation)[1]) {
  result_cor_coagulation[i,"avg"] <- mean(result_cor_coagulation[i,1:dim(result_cor_coagulation)[2]-1])
}

result_cor_platelet$avg <- NA
result_cor_platelet <- as.matrix(result_cor_platelet)
for (i in 1:dim(result_cor_platelet)[1]) {
  result_cor_platelet[i,"avg"] <- mean(result_cor_platelet[i,1:dim(result_cor_platelet)[2]-1])
}



## 相关性分析 ----
result <- plot_expression_correlation(exprMatrix = result_check$data$erythrocyte[,-1:-2],displayNumbers = T,corMethod = "pearson")
result <- plot_expr_corrplot(exprMatrix = result_check$data$erythrocyte[,-1:-2])
result <- plot_expression_correlation(exprMatrix = result_check$correlation$erythrocyte$r,displayNumbers = T,input_type = "correlation")
plot_pvalue_distribution(pvalue_matrix =  result_check$correlation$erythrocyte$P)

plot_stat_distribution(data =  result_check$rawdata,statistic = "correlation")
test_result <- Hmisc::rcorr(as.matrix(t( result_check$rawdata)), type = "pearson")
statistic <- "correlation"
values <- if (statistic == "pvalue") test_result$P else test_result$r

data_ggplot <- data.frame(value = values[upper.tri(values)])
ggplot(data_ggplot,aes(x = value)) + 
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black") + 
  geom_density(color = "red")

# erythrocyte
plot_protein_by_sample(data = result_correct$rawdata[rownames(result_correct$correct_data)%in%result_correct$marker_list$erythrocyte,])
plot_protein_by_sample(data = result_correct$correct_data[rownames(result_correct$correct_data)%in%result_correct$marker_list$erythrocyte,])
# coagulation
plot_protein_by_sample(data = result_correct$rawdata[rownames(result_correct$correct_data)%in%result_correct$marker_list$coagulation,])
plot_protein_by_sample(data = result_correct$correct_data[rownames(result_correct$correct_data)%in%result_correct$marker_list$coagulation,])
# platelet
plot_protein_by_sample(data = result_correct$rawdata[rownames(result_correct$correct_data)%in%result_correct$marker_list$platelet,])
plot_protein_by_sample(data = result_correct$correct_data[rownames(result_correct$correct_data)%in%result_correct$marker_list$platelet,])
plot_protein(data = df,protein = result_check$marker_list$erythrocyte)
result <- limma_proteomics_analysis(expr_matrix = log2(df),
                                    group_matrix = data_group,
                                    compare = c("AD","Ctrl"),
                                    p_type = "raw")
plot_1 <- create_volcano_plot(limma_result = result,
                              p_type = "raw",gene_col = "Protein",
                              group_names = c("AD","Ctrl"),logFC_cutoff = 0)
plot_1

result_cv_coa <- get_cv(raw_data = df,protein = list_coagulation$GN)
result_cv_ery <- get_cv(raw_data = df,protein = list_erythrocyte$GN)
result_cv_pla <- get_cv(raw_data = df,protein = list_platelet$GN)
result_cv_other <- get_cv(raw_data = df[!rownames(df)%in%
                                          c(list_coagulation$GN,
                                            list_erythrocyte$GN,
                                            list_platelet$GN),])
result_cv_coa$type <- "coagulation"
result_cv_ery$type <- "erythrocyte"
result_cv_pla$type <- "platelet"
result_cv_other$type <- "other protein"
result_cv <- rbind(result_cv_coa,result_cv_ery,result_cv_pla,result_cv_other)
result_cv$type <- factor(result_cv$type,
                         levels = c("coagulation",
                                    "erythrocyte",
                                    "platelet",
                                    "other protein"))

ggplot(result_cv,aes(x = type , y = CV, fill = type )) +
  geom_violin() +
  geom_boxplot(fill = "white",width = 0.2) +
  scale_fill_manual(values = c("coagulation" = "#E64B35FF",
                               "erythrocyte" = "#F39B7FFF",
                               "platelet" = "#7E6148FF",
                               "other protein" = "#3C5488FF")) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13))

# VN图测试----
# 生成示例数据
group1 <- letters[1:6]  # a-j
group2 <- letters[5:15]  # e-o

# 基本使用
result <- venn_plot(group1, group2)

# 自定义参数版本
venn_plot(
  set1 = group1,
  set2 = group2,
  categories = c("pre", "post"),
  title = "基因表达重叠分析",
  colors = c("#4daf4a", "#984ea3"),
  alpha = 0.6,
  print.mode = "raw"
)

# 查看交集信息
cat("总唯一元素数量:", result$total_unique, "\n")
cat("重叠元素数量:", result$overlap_count, "\n")
cat("具体重叠元素:", paste(result$overlap_elements, collapse = ", "))

