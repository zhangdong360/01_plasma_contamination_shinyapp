data_check <- function(data, data_group, cutoff = 0.9, 
                       group1 = NULL, group2 = NULL,
                       custom_erythrocyte = NULL,
                       custom_coagulation = NULL,
                       custom_platelet = NULL){
  # 函数说明 ----
  # data为输入数据矩阵，要求为样本为列，蛋白为行，矩阵中无非数值数据
  # data_group为样本分组，要求为：id列为样本名，和colnames(data)相同，group为输入样本的生物学分组
  # cutoff为相关性阈值，默认0.9
  
  # 加载必要的库
  library(Hmisc)
  library(readxl)
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(dplyr)
  library(stringr)
  
  # 检查输入数据的有效性
  if (!all(data_group$id %in% colnames(data))) {
    stop("Error: data_group中的id列与data的列名不匹配。")
  }
  
  # 定义分析标记物的函数
  # Function to filter markers and perform analysis
  analyze_markers <- function(data = data_ggplot, marker_list, title, filename_marker) {
    # 删除生物学分组中有差异的marker
    target_proteins <- marker_list$GN
    pattern <- paste0("\\b(", paste(target_proteins, collapse = "|"), ")\\b")
    data_filtered <- data_ggplot %>% filter(str_detect(key, pattern))
    
    # 绘制并保存箱线图
    p <- ggplot(data_filtered, aes(x = key, y = log2(value + 1), fill = group)) + 
      geom_boxplot() +
      theme_classic() +
      stat_compare_means(aes(group = group), label = "p.signif", method = "wilcox.test", size = 3) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black", size = 10),
            axis.text.y = element_text(hjust = 1, colour = "black", size = 10),
            plot.title = element_text(hjust = 0.5)) +
      labs(title = title)
    
    # 提取统计结果
    ggplot_build_obj <- ggplot_build(p)
    stat_results <- ggplot_build_obj$data[[2]]
    stat_results_with_keys <- data_filtered %>% distinct(key) %>% bind_cols(stat_results)
    keys_with_high_pvalue <- stat_results_with_keys %>% filter(p > 0.05) %>% pull(key)
    
    return(list(data = data_filtered, keys_with_high_pvalue = keys_with_high_pvalue,plot = p))
  }
  
  # 定义相关性分析函数
  # Correlation analysis function
  correlate_markers <- function(data_filtered) {
    data_spread <- spread(data_filtered, key, value)
    numeric_data <- data_spread %>% dplyr::select(-id, -group)
    corr_results <- Hmisc::rcorr(as.matrix(numeric_data), type = "pearson")
    return(list(r = corr_results$r, P = corr_results$P))
  }
  
  # 定义基于p值和相关性过滤标记物的函数
  # Filter markers based on p-value and correlation
  filter_markers <- function(keys_with_high_pvalue, corr_matrix_r) {
    # 提取相关性高的标记物
    high_corr_keys <- colnames(corr_matrix_r)[apply(corr_matrix_r, 2, function(x) any(abs(x) > cutoff & abs(x) < 1))]
    filtered_keys <- intersect(keys_with_high_pvalue, high_corr_keys)
    result <- list(high_corr_keys = high_corr_keys,
                   filtered_keys = filtered_keys)
    return(result)
  }
  
  # 定义绘制水平图的函数
  # Plot and save level plots for filtered markers
  plot_level <- function(data_filtered, filtered_keys, title) {
    if (length(filtered_keys) == 0) {
      warning(paste("No markers found for", title))
      return()
    }
    data_filtered <- data_filtered %>% filter(key %in% filtered_keys)
    p <- ggplot(data_filtered, aes(x = id, y = log2(value + 1), fill = id)) + 
      geom_boxplot() +
      theme_classic() +
      stat_compare_means(aes(group = id), label = "p.format", method = "anova", size = 3) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black", size = 10),
            axis.text.y = element_text(hjust = 1, colour = "black", size = 10),
            plot.title = element_text(hjust = 0.5),
            legend.position="none") +
      labs(title = paste(title, "level"))
    return(p)
  }
  
  # 收集所有统计信息
  analyze_markers_full <- function(data_ggplot, marker_list,analyze_markers_result,correlation_result) {
    target_proteins <- marker_list$GN
    data_filtered <- data_ggplot %>% 
      filter(key %in% target_proteins)
    ggplot_build_obj <- ggplot_build(analyze_markers_result$plot)
    analyze_markers_result$keys_with_high_pvalue
    stat_results <- ggplot_build_obj$data[[2]]
    stat_results <- subset(stat_results,select = c(p))
    # 计算统计信息
    stats <- data.frame(key = target_proteins,
                        exists = target_proteins%in%unique(data_filtered$key),
                        DE = ifelse(target_proteins%in%analyze_markers_result$keys_with_high_pvalue,FALSE,TRUE),
                        correlation = ifelse(target_proteins%in%correlation_result$high_corr_keys,FALSE,TRUE))
    stats[stats$exists%in% FALSE,"DE"] <- NA
    stats[stats$exists%in% FALSE,"correlation"] <- NA
    return(stats)
  }

  # Mann panel data input ----
  library(readxl)
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(dplyr)
  library(stringr)
  
  result <- list()
  result$rawdata <- data
  result$group <- data_group
  list_erythrocyte <- if (!is.null(custom_erythrocyte)) {
    data.frame(GN = custom_erythrocyte, ContaminationType = "erythrocyte")
  } else { 
    data.frame(GN = c("HBA1","HBB","CA1","HBD","PRDX2","CA2",
                      "CAT","SLC4A1","BLVRB","SPTA1","SPTB",
                      "ANK1","GAPDH","SOD1","PRDX6","BPGM",
                      "ACTB","ACTG1","SELENBP1","EPB41",
                      "LDHB","EIF5A","EIF5A2","EIF5AL1",
                      "PNP","ALDOA","TPI1","NIF3L1",
                      "UBC","UBB","RPS27A","UBA52",
                      "UBBP4","HBE1","HSPA8","PGK1","ALAD"),
               ContaminationType = "erythrocyte") }
  list_platelet <- if (!is.null(custom_platelet)) {
    data.frame(GN = custom_platelet, ContaminationType = "platelet")
  } else {
    data.frame(GN = c("FLNA","TLN1","MYH9","ACTB","VCL",
                      "ACTN1","TPM4","THBS1","TUBB1",
                      "YWHAZ","GSN","TUBA1B","ITGA2B",
                      "F13A1","PFN1","ITGB3","TAGLN2",
                      "FERMT3","RAP1B","PLEK","PPBP",
                      "GAPDH","MMRN1","MYL6","CFL1",
                      "PARVB","SDPR","TUBB4B","TMSB4X","PKM"),
               ContaminationType = "platelet")}
  list_coagulation <- if (!is.null(custom_coagulation)) {
    data.frame(GN = custom_coagulation, ContaminationType = "coagulation")
  } else {
    data.frame(GN = c("FGB","FGG","FGA","F13A1",
                      "SERPINC1","PPBP","F2",
                      "GP1BA","ECM1","CLU",
                      "DSP","WDR1","ATRN",
                      "GP5","SERPINA5","C1RL",
                      "MAN1A1","F11","KNG1",
                      "BCHE","GPLD1","THBS1",
                      "MASP1","TENM4","HSPG2",
                      "APOC3","CDH5","LTF","CST3",
                      "PROC","PF4","PF4V1"),
               ContaminationType = "coagulation")}
  # 检查缺失基因的函数
  check_missing_genes <- function(data, gene_list) {
    missing_genes <- gene_list[!gene_list %in% rownames(data)]
    return(missing_genes)
  }
  
  # 检查缺失的红细胞基因
  missing_erythrocyte_genes <- check_missing_genes(data, list_erythrocyte$GN)
  # 检查缺失的血小板基因
  missing_platelet_genes <- check_missing_genes(data, list_platelet$GN)
  # 检查缺失的凝血基因
  missing_coagulation_genes <- check_missing_genes(data, list_coagulation$GN)
  ## determine panel ----
  # 转换数据格式
  data <- as.data.frame(t(data))
  data$id <- rownames(data)
  data_merge <- merge(data_group, data, by = "id")
  data <- subset(data,select = -c(id))
  # 转换为长格式
  data_ggplot <- tidyr::gather(data_merge,
                               key = "key", 
                               value = "value", 
                               -c("id", "group"))
  
  ## Perform analysis for erythrocyte, platelet, and coagulation markers ----
  # 分析红细胞、血小板和凝血标记物
  erythrocyte_results <- analyze_markers(data_ggplot, list_erythrocyte,
                                         "Erythrocyte marker", 
                                         "/boxplot_Erythrocyte_markers")
  platelet_results <- analyze_markers(data_ggplot, list_platelet, 
                                      "Platelet marker",
                                      "/boxplot_Platelet_markers")
  coagulation_results <- analyze_markers(data_ggplot, list_coagulation, 
                                         "Coagulation marker", 
                                         "/boxplot_Coagulation_markers")
  
  ## Perform correlation analysis ----
  # 相关性分析
  corr_matrix_erythrocyte <- correlate_markers(erythrocyte_results$data)
  corr_matrix_platelet <- correlate_markers(platelet_results$data)
  corr_matrix_coagulation <- correlate_markers(coagulation_results$data)
  
  result_correlation <- list(
    erythrocyte = corr_matrix_erythrocyte,
    coagulation = corr_matrix_coagulation,
    platelet = corr_matrix_platelet
  )
  
  ## Filter markers based on p-value and correlation ----
  # 过滤标记物
  filtered_keys_erythrocyte <- filter_markers(erythrocyte_results$keys_with_high_pvalue,
                                              corr_matrix_erythrocyte$r)
  filtered_keys_platelet <- filter_markers(platelet_results$keys_with_high_pvalue, 
                                           corr_matrix_platelet$r)
  filtered_keys_coagulation <- filter_markers(coagulation_results$keys_with_high_pvalue,
                                              corr_matrix_coagulation$r)
  
  
  marker_list <- list(
    erythrocyte = filtered_keys_erythrocyte$filtered_keys,
    coagulation = filtered_keys_coagulation$filtered_keys,
    platelet = filtered_keys_platelet$filtered_keys
  )
  
  ## Plot and save level plots and result for filtered markers ----
  # 绘制并保存水平图
  Erythrocyte_pannel <- plot_level(erythrocyte_results$data, 
                                   marker_list$erythrocyte, 
                                   "Erythrocyte pannel")
  Platelet_pannel <-  plot_level(platelet_results$data,
                                 marker_list$platelet,
                                 "Platelet pannel")
  Coagulation_pannel <- plot_level(coagulation_results$data, 
                                   marker_list$coagulation, 
                                   "Coagulation pannel")
  stats_erythrocyte <- analyze_markers_full(data_ggplot = data_ggplot,
                                            marker_list = list_erythrocyte, 
                                            analyze_markers_result = erythrocyte_results,
                                            correlation_result = filtered_keys_erythrocyte)
  stats_platelet <- analyze_markers_full(data_ggplot = data_ggplot,
                                         marker_list = list_platelet, 
                                         analyze_markers_result = platelet_results,
                                         correlation_result = filtered_keys_platelet)
  stats_coagulation <- analyze_markers_full(data_ggplot = data_ggplot,
                                            marker_list = list_coagulation, 
                                            analyze_markers_result = coagulation_results,
                                            correlation_result = filtered_keys_coagulation)
  
  # 保存结果
  result <- list(
    rawdata = t(data),
    group = data_group,
    data = list(
      erythrocyte = spread(erythrocyte_results$data, key, value),
      coagulation = spread(coagulation_results$data, key, value),
      platelet = spread(platelet_results$data, key, value)
    ),
    plot_marker = list(erythrocyte = Erythrocyte_pannel,
                coagulation = Coagulation_pannel,
                platelet = Platelet_pannel),
    plot_contamination = list(erythrocyte = erythrocyte_results$plot,
                              coagulation = coagulation_results$plot,
                              platelet = platelet_results$plot),
    gene = list(
      erythrocyte = erythrocyte_results$keys_with_high_pvalue,
      coagulation = coagulation_results$keys_with_high_pvalue,
      platelet = platelet_results$keys_with_high_pvalue
    ),
    missing_genes = list(missing_erythrocyte_genes = missing_erythrocyte_genes,
                         missing_coagulation_genes = missing_coagulation_genes,
                         missing_platelet_genes = missing_platelet_genes),
    correlation = result_correlation,
    marker_list = marker_list,
    marker_stats = list(
      stats_erythrocyte = stats_erythrocyte,
      stats_platelet = stats_platelet,
      stats_coagulation = stats_coagulation)
  )
  return(result)
}
