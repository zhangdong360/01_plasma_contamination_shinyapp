data_correct <- function(data, type = "all") {
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
    list1 <- for_mean(data$data$erythrocyte[, colnames(data$data$erythrocyte) %in% data$gene$erythrocyte])
    list2 <- for_mean(data$data$platelet[, colnames(data$data$platelet) %in% data$gene$platelet])
    list3 <- for_mean(data$data$coagulation[, colnames(data$data$coagulation) %in% data$gene$coagulation])
    
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
  
  if (type == "all") {
    for (i in 1:nrow(ndata1)) {
      y <- rawdata[i, ] # expression across samples
      y <- t(y)
      x_eryth <- smpl2$erythrocyte
      x_plate <- smpl2$Platelet
      x_coagu <- smpl2$coagulation
      a <- summary(rlm(log2(y + 1) ~ x_eryth + x_plate + x_coagu, maxit = 30))
      for (j in 2:nrow(a$coefficients)) {
        b[i, j - 1] = a$coefficients[j, 1]
      }
      # new residuals:
      ndata1[i, ] <- log2(y + 1) - (a[["coefficients"]]["(Intercept)","Value"] + 
                                      x_eryth *a[["coefficients"]]["x_eryth","Value"] + 
                                      x_plate * a[["coefficients"]]["x_plate","Value"] +
                                      x_coagu *a[["coefficients"]]["x_coagu","Value"] ) + 
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
                                      x_eryth *a[["coefficients"]]["x_eryth","Value"]) + 
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
                                      x_plate *a[["coefficients"]]["x_plate","Value"]) + 
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
                                      x_coagu *a[["coefficients"]]["x_coagu","Value"]) + 
        a[["coefficients"]]["(Intercept)","Value"]
    }
  }
  result <- 2 ^ ndata1 - 1
  data$correct_data <- result
  data$contamination_level <- smpl2
  return(data)
}
