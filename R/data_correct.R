data_correct <- function(data, 
                         type = "all",constraint = 1.2,
                         erythrocyte_marker,
                         coagulation_marker,
                         platelet_marker) {
  # 函数说明 ----
  # data为输入数据，type为矫正类型，默认为"all"，表示矫正所有污染类型
  # type可以为erythrocyte，coagulation，platelet中任意几种
  
  library(MASS)
  # 定义计算均值的函数
  
  mean_2 <- function(data) {
    for_mean <- function(data) {
      # 标准化数据
      data <- log2(data + 1)
      M <- as.data.frame(t(data))
      a <- M
      for (i in 1:nrow(M)) {
        for (j in 1:ncol(M)) {
          a[i, j] <- M[i, j] - mean(as.numeric(M[i, ]), na.rm = TRUE)
        }
      }
      mean <- sapply(a, mean, na.rm = TRUE)
      return(mean)
    }
    # 增加检测data$data$erythrocyte是否为NULL
    rownames(data$data$erythrocyte) <- data$data$erythrocyte$id
    rownames(data$data$coagulation) <- data$data$coagulation$id
    rownames(data$data$platelet) <- data$data$platelet$id
    
    # 计算污染水平
    list1 <- for_mean(data$data$erythrocyte[, colnames(data$data$erythrocyte) %in% erythrocyte_marker])
    list2 <- for_mean(data$data$platelet[, colnames(data$data$platelet) %in% platelet_marker])
    list3 <- for_mean(data$data$coagulation[, colnames(data$data$coagulation) %in% coagulation_marker])
    
    smpl2 <- data.frame(
      erythrocyte = list1,
      Platelet = list2,
      coagulation = list3
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
    result_cor_erythrocyte[i,"avg"] <- mean(result_cor_erythrocyte[i,round(dim(result_cor_erythrocyte)[2]*percentage_erythrocyte):dim(result_cor_erythrocyte)[2]-1])
  }
  
  result_cor_coagulation$avg <- NA
  result_cor_coagulation <- as.matrix(result_cor_coagulation)
  for (i in 1:dim(result_cor_coagulation)[1]) {
    result_cor_coagulation <- result_cor_coagulation[,order(result_cor_coagulation[i,])]
    result_cor_coagulation[i,"avg"] <- mean(result_cor_coagulation[i,round(dim(result_cor_coagulation)[2]*percentage_coagulation):dim(result_cor_coagulation)[2]-1])
  }
  
  result_cor_platelet$avg <- NA
  result_cor_platelet <- as.matrix(result_cor_platelet)
  for (i in 1:dim(result_cor_platelet)[1]) {
    result_cor_platelet <- result_cor_platelet[,order(result_cor_platelet[i,])]
    result_cor_platelet[i,"avg"] <- mean(result_cor_platelet[i,round(dim(result_cor_platelet)[2]*percentage_platelet):dim(result_cor_platelet)[2]-1])
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
      library(MASS)
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
      library(MASS)
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
      library(MASS)
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
