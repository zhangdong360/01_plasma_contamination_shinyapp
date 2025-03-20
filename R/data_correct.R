data_correct <- function(data, 
                         type = "all",constraint = 1.2,
                         erythrocyte_marker,constraint_erythrocyte = 0.95,
                         coagulation_marker,constraint_coagulation = 0.95,
                         platelet_marker,constraint_platelet = 0.95) {
  # 函数说明 ----
  # data为输入数据，type为矫正类型，默认为"all"，表示矫正所有污染类型
  # type可以为erythrocyte，coagulation，platelet中任意几种
  # 强制参数在(0,1)范围内（包含边界检查）
  constraint_erythrocyte_original <- max(0, min(1, constraint_erythrocyte))
  constraint_platelet_original <- max(0, min(1, constraint_platelet))
  constraint_coagulation_original <- max(0, min(1, constraint_coagulation)) 
  # 当参数被修正时发出警告
  if(constraint_erythrocyte != constraint_erythrocyte_original) 
    warning("constraint_erythrocyte clamped to [0,1]")
  if(constraint_platelet != constraint_platelet_original) 
    warning("constraint_platelet clamped to [0,1]")
  if(constraint_coagulation != constraint_coagulation_original) 
    warning("constraint_coagulation clamped to [0,1]")
  library(MASS)
  # 检查必要的marker参数是否为空
  if ("all" %in% type) {
    if (missing(erythrocyte_marker) || is.null(erythrocyte_marker) || length(erythrocyte_marker) == 0) {
      stop("当type为'all'时，参数erythrocyte_marker不能为空")
    }
    if (missing(coagulation_marker) || is.null(coagulation_marker) || length(coagulation_marker) == 0) {
      stop("当type为'all'时，参数coagulation_marker不能为空")
    }
    if (missing(platelet_marker) || is.null(platelet_marker) || length(platelet_marker) == 0) {
      stop("当type为'all'时，参数platelet_marker不能为空")
    }
  } else {
    for (t in type) {
      if (t == "erythrocyte") {
        if (missing(erythrocyte_marker) || is.null(erythrocyte_marker) || length(erythrocyte_marker) == 0) {
          stop("当type包含'erythrocyte'时，参数erythrocyte_marker不能为空")
        }
      } else if (t == "coagulation") {
        if (missing(coagulation_marker) || is.null(coagulation_marker) || length(coagulation_marker) == 0) {
          stop("当type包含'coagulation'时，参数coagulation_marker不能为空")
        }
      } else if (t == "platelet") {
        if (missing(platelet_marker) || is.null(platelet_marker) || length(platelet_marker) == 0) {
          stop("当type包含'platelet'时，参数platelet_marker不能为空")
        }
      }
    }
  }
  # 定义计算均值的函数
  
  mean_2 <- function(data) {
    for_mean <- function(data) {
      # 标准化数据
      data <- log2(data + 1)
      M <- as.data.frame(t(data))
      a <- as.matrix(M)
      a <- M - rowMeans(M)
      mean <- sapply(a, mean, na.rm = TRUE)
      return(mean)
    }
    # 增加检测data$data$erythrocyte是否为NULL
    # erythrocyte_df
    erythrocyte_df <- data.frame(data$data$erythrocyte)
    rownames(erythrocyte_df) <- erythrocyte_df$id
    # platelet_df
    platelet_df <- data.frame(data$data$platelet)
    rownames(platelet_df) <- platelet_df$id
    # coagulation_df
    coagulation_df <- data.frame(data$data$coagulation)
    rownames(coagulation_df) <- coagulation_df$id
    
    rownames(data$data$erythrocyte) <- data$data$erythrocyte$id
    rownames(data$data$coagulation) <- data$data$coagulation$id
    rownames(data$data$platelet) <- data$data$platelet$id
    
    # 计算污染水平
    list_erythrocyte <- for_mean(erythrocyte_df[, colnames(erythrocyte_df) %in% erythrocyte_marker])
    list_platelet <- for_mean(platelet_df[, colnames(platelet_df) %in% platelet_marker])
    list_coagulation <- for_mean(coagulation_df[, colnames(coagulation_df) %in% coagulation_marker])
    
    smpl2 <- data.frame(
      erythrocyte = list_erythrocyte,
      Platelet = list_platelet,
      coagulation = list_coagulation
    )
    return(smpl2)
  }
  
  smpl2 <- mean_2(data)
  rawdata <- as.data.frame(data$rawdata)
  ndata1 <- matrix(nrow = nrow(rawdata), ncol = ncol(rawdata))
  rownames(ndata1) <- rownames(rawdata)
  colnames(ndata1) <- colnames(rawdata)
  b <- matrix(nrow = nrow(rawdata), ncol = ncol(smpl2))
  rownames(b) <- rownames(rawdata)
  colnames(b) <- colnames(smpl2)
  # 计算约束系数 ----
  result_cor <- test_result <- Hmisc::rcorr(as.matrix(t(data$rawdata)), type = "pearson")
  result_cor <- as.data.frame(result_cor$r)
  
  result_cor_erythrocyte <- abs(result_cor[,colnames(result_cor)%in%c(data$marker_list$erythrocyte)])
  result_cor_coagulation <- abs(result_cor[,colnames(result_cor)%in%c(data$marker_list$coagulation)])
  result_cor_platelet <- abs(result_cor[,colnames(result_cor)%in%c(data$marker_list$platelet)])
  result_cor_erythrocyte$avg <- NA
  result_cor_erythrocyte <- as.matrix(result_cor_erythrocyte)
  for (i in 1:dim(result_cor_erythrocyte)[1]) {
    result_cor_erythrocyte <- result_cor_erythrocyte[,order(result_cor_erythrocyte[i,])]
    start_col_erythrocyte <- min(max(round(dim(result_cor_erythrocyte)[2]*constraint_erythrocyte_original),1),dim(result_cor_erythrocyte)[2]-2)
    result_cor_erythrocyte[i,"avg"] <- mean(result_cor_erythrocyte[i,start_col_erythrocyte:dim(result_cor_erythrocyte)[2]-1])
  }
  
  result_cor_coagulation$avg <- NA
  result_cor_coagulation <- as.matrix(result_cor_coagulation)
  for (i in 1:dim(result_cor_coagulation)[1]) {
    result_cor_coagulation <- result_cor_coagulation[,order(result_cor_coagulation[i,])]
    start_col_coagulation <- min(max(round(dim(result_cor_coagulation)[2]*constraint_coagulation_original),1),
                                 dim(result_cor_coagulation)[2]-2)
    result_cor_coagulation[i,"avg"] <- mean(result_cor_coagulation[i,start_col_coagulation:dim(result_cor_coagulation)[2]-1])
  }
  
  result_cor_platelet$avg <- NA
  result_cor_platelet <- as.matrix(result_cor_platelet)
  for (i in 1:dim(result_cor_platelet)[1]) {
    result_cor_platelet <- result_cor_platelet[,order(result_cor_platelet[i,])]
    start_col_platelet <- min(max(round(dim(result_cor_platelet)[2]*constraint_platelet_original),1),dim(result_cor_platelet)[2]-2)
    result_cor_platelet[i,"avg"] <- mean(result_cor_platelet[i,start_col_platelet:dim(result_cor_platelet)[2]-1])
  }
  
  if (type == "all") {
    for (i in 1:nrow(ndata1)) {
      y <- rawdata[i, ] # expression across samples
      y <- t(y)
      x_eryth <- smpl2$erythrocyte
      x_plate <- smpl2$Platelet
      x_coagu <- smpl2$coagulation
      Constraint_factor_eryth <- result_cor_erythrocyte[colnames(y),"avg"]
      Constraint_factor_plate <- result_cor_platelet[colnames(y),"avg"]
      Constraint_factor_coagu <- result_cor_coagulation[colnames(y),"avg"]
      a <- summary(rlm(log2(y + 1) ~ x_eryth + x_plate + x_coagu, maxit = 30))
      for (j in 2:nrow(a$coefficients)) {
        b[i, j - 1] = a$coefficients[j, 1]
      }
      # new residuals:
      ndata1[i, ] <- log2(y + 1) - (a[["coefficients"]]["(Intercept)","Value"] + constraint * (
        x_eryth *a[["coefficients"]]["x_eryth","Value"]* Constraint_factor_eryth + 
          x_plate * a[["coefficients"]]["x_plate","Value"]* Constraint_factor_plate +
          x_coagu *a[["coefficients"]]["x_coagu","Value"]* Constraint_factor_coagu) ) + 
        a[["coefficients"]]["(Intercept)","Value"]
    }
  }
  if (type == "erythrocyte") {
    for (i in 1:nrow(ndata1)) {
      y <- rawdata[i, ] # expression across samples
      y <- t(y)
      x_eryth <- smpl2$erythrocyte
      a <- summary(rlm(log2(y + 1) ~ x_eryth, maxit = 30))
      for (j in 2:nrow(a$coefficients)) {
        b[i, j - 1] = a$coefficients[j, 1]
      }
      # new residuals:
      ndata1[i, ] <- log2(y + 1) - (a[["coefficients"]]["(Intercept)","Value"] + 
                                      constraint * x_eryth *a[["coefficients"]]["x_eryth","Value"]) + 
        a[["coefficients"]]["(Intercept)","Value"]
    }
  }
  if (type == "platelet") {
    for (i in 1:nrow(ndata1)) {
      y <- rawdata[i, ] # expression across samples
      y <- t(y)
      x_plate <- smpl2$Platelet
      # x_coagu <- smpl2$coagulation
      a <- summary(rlm(log2(y + 1) ~ x_plate, maxit = 30))
      for (j in 2:nrow(a$coefficients)) {
        b[i, j - 1] = a$coefficients[j, 1]
      }
      # new residuals:
      ndata1[i, ] <- log2(y + 1) - (a[["coefficients"]]["(Intercept)","Value"] + 
                                      constraint * x_plate *a[["coefficients"]]["x_plate","Value"]) + 
        a[["coefficients"]]["(Intercept)","Value"]
    }
  }
  if (type == "coagulation") {
    for (i in 1:nrow(ndata1)) {
      y <- rawdata[i, ] # expression across samples
      y <- t(y)
      x_coagu <- smpl2$coagulation
      a <- summary(rlm(log2(y + 1) ~ x_coagu, maxit = 30))
      for (j in 2:nrow(a$coefficients)) {
        b[i, j - 1] = a$coefficients[j, 1]
      }
      # new residuals:
      ndata1[i, ] <- log2(y + 1) - (a[["coefficients"]]["(Intercept)","Value"] + 
                                      constraint * x_coagu *a[["coefficients"]]["x_coagu","Value"]) + 
        a[["coefficients"]]["(Intercept)","Value"]
    }
  }
  result <- 2 ^ ndata1 - 1
  data$correct_data <- result
  data$contamination_level <- smpl2
  return(data)
}
