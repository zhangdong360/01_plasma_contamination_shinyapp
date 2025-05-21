data_correct <- function(data, 
                         type = c("coagulation","erythrocyte","platelet"),constraint = 1,
                         erythrocyte_marker,constraint_erythrocyte = 0.95,
                         coagulation_marker,constraint_coagulation = 0.95,
                         platelet_marker,constraint_platelet = 0.95) {
  # 函数说明 ----
  # data为输入数据，type为矫正类型，默认为c("coagulation","erythrocyte","platelet")，表示矫正所有污染类型
  # type可以为coagulation, erythrocyte, platelet中任意几种
  # constraint_*:参数为计算约束因子的最小分位数，默认为相关系数从小到大排序，计算后95%-100%相关系数的平均数
  # 强制constraint_*参数在(0,1)范围内（包含边界检查）
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
  # 定义计算均值的函数
  mean_2 <- function(data) {
    for_mean <- function(data) {
      # 标准化数据（不再需要转置）
      data <- log2(data + 1)
      # 直接使用原始基因 × 样本的结构
      M <- as.data.frame(data)
      # 计算行（基因）均值并中心化
      a <- M - rowMeans(M)
      # 对列（样本）取均值
      mean <- apply(a,2, mean, na.rm = TRUE)
      return(mean)
    }
    # 计算污染水平
    # list_erythrocyte <- for_mean(data$rawdata[rownames(data$rawdata) %in% erythrocyte_marker,])
    # list_platelet <- for_mean(data$rawdata[rownames(data$rawdata) %in% platelet_marker,])
    # list_coagulation <- for_mean(data$rawdata[rownames(data$rawdata) %in% coagulation_marker,])
    
    list_erythrocyte <-  if (any(rownames(data$rawdata) %in% erythrocyte_marker)) {
      for_mean(data$rawdata[rownames(data$rawdata) %in% erythrocyte_marker,,drop = FALSE])
    } else {
      rep(NA,dim(data$rawdata)[2])
    }
    list_platelet <- if (any(rownames(data$rawdata) %in% platelet_marker)) {
      for_mean(data$rawdata[rownames(data$rawdata) %in% platelet_marker,, drop = FALSE])
    } else {
      rep(NA,dim(data$rawdata)[2])
    }
    list_coagulation <- if (any(rownames(data$rawdata) %in% coagulation_marker)) {
      for_mean(data$rawdata[rownames(data$rawdata) %in% coagulation_marker,,drop = FALSE])
    } else {
      rep(NA,dim(data$rawdata)[2])
    }
    
    components  <- data.frame(
      erythrocyte = list_erythrocyte,
      Platelet = list_platelet,
      coagulation = list_coagulation
    )
    return(components )
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
  # 初始化约束因子为 NULL
  result_cor_erythrocyte <- result_cor_platelet <- result_cor_coagulation <- NULL
  # 动态计算各类型约束因子
  ## ery ----
  if (!is.null(erythrocyte_marker) && length(erythrocyte_marker) > 0) {
    result_cor_erythrocyte <- abs(result_cor[, colnames(result_cor) %in% erythrocyte_marker, drop = FALSE])
    result_cor_erythrocyte$avg <- NA
    result_cor_erythrocyte <- as.matrix(result_cor_erythrocyte)
    for (i in 1:dim(result_cor_erythrocyte)[1]) {
      result_cor_erythrocyte <- result_cor_erythrocyte[,order(result_cor_erythrocyte[i,])]
      if (dim(result_cor_erythrocyte)[2] < 3) {
        result_cor_erythrocyte[i,"avg"] <- mean(result_cor_erythrocyte[i,],na.rm = T)
      }else if (dim(result_cor_erythrocyte)[2] > 2) {
        start_col_erythrocyte <- min(max(round(dim(result_cor_erythrocyte)[2]*constraint_erythrocyte_original),1),
                                     dim(result_cor_erythrocyte)[2]-2)
        result_cor_erythrocyte[i,"avg"] <- mean(result_cor_erythrocyte[i,start_col_erythrocyte:dim(result_cor_erythrocyte)[2]-1])
        }
      }
    }
  ## cog ----
  if (!is.null(coagulation_marker) && length(coagulation_marker) > 0) {
    result_cor_coagulation <- abs(result_cor[, colnames(result_cor) %in% coagulation_marker, drop = FALSE])
    result_cor_coagulation$avg <- NA
    result_cor_coagulation <- as.matrix(result_cor_coagulation)
    for (i in 1:dim(result_cor_coagulation)[1]) {
      result_cor_coagulation <- result_cor_coagulation[,order(result_cor_coagulation[i,])]
      if (dim(result_cor_coagulation)[2] < 3) {
        result_cor_coagulation[i,"avg"] <- mean(result_cor_coagulation[i,],na.rm = T)
      }else if (dim(result_cor_coagulation)[2] > 2){
        start_col_coagulation <- min(max(round(dim(result_cor_coagulation)[2]*constraint_coagulation_original),1),
                                     dim(result_cor_coagulation)[2]-2)
        result_cor_coagulation[i,"avg"] <- mean(result_cor_coagulation[i,start_col_coagulation:dim(result_cor_coagulation)[2]-1])
      }
    }
  }
  
  ## platelet ----
  if (!is.null(platelet_marker) && length(platelet_marker) > 0) {
    result_cor_platelet <- abs(result_cor[, colnames(result_cor) %in% platelet_marker, drop = FALSE])
    result_cor_platelet$avg <- NA
    result_cor_platelet <- as.matrix(result_cor_platelet)
    for (i in 1:dim(result_cor_platelet)[1]) {
      result_cor_platelet <- result_cor_platelet[,order(result_cor_platelet[i,])]
      if (dim(result_cor_platelet)[2] < 3){
        result_cor_platelet[i,"avg"] <- mean(result_cor_platelet[i,],na.rm = T)
      }else if (dim(result_cor_platelet)[2] > 2){
        start_col_platelet <- min(max(round(dim(result_cor_platelet)[2]*constraint_platelet_original),1),
                                  dim(result_cor_platelet)[2]-2)
        result_cor_platelet[i,"avg"] <- mean(result_cor_platelet[i,start_col_platelet:dim(result_cor_platelet)[2]-1])
      }
    }
  }
  
  # correct ----
  if (setequal(type, c("coagulation","erythrocyte","platelet"))) {
    ## all ----
    for (i in 1:nrow(ndata1)) {
      y <- rawdata[i, ] # expression across samples
      y <- t(y)
      smpl2 <- smpl2[rownames(y),]
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
  } else if (setequal(type, c("coagulation","platelet"))) {
    ## coagulation platelet ----
    for (i in 1:nrow(ndata1)) {
      y <- rawdata[i, ] # expression across samples
      y <- t(y)
      smpl2 <- smpl2[rownames(y),]
      x_plate <- smpl2$Platelet
      x_coagu <- smpl2$coagulation
      Constraint_factor_plate <- result_cor_platelet[colnames(y),"avg"]
      Constraint_factor_coagu <- result_cor_coagulation[colnames(y),"avg"]
      a <- summary(rlm(log2(y + 1) ~ x_plate + x_coagu, maxit = 30))
      for (j in 2:nrow(a$coefficients)) {
        b[i, j - 1] = a$coefficients[j, 1]
      }
      # new residuals:
      ndata1[i, ] <- log2(y + 1) - (a[["coefficients"]]["(Intercept)","Value"] + constraint * (
        x_plate * a[["coefficients"]]["x_plate","Value"]* Constraint_factor_plate +
          x_coagu *a[["coefficients"]]["x_coagu","Value"]* Constraint_factor_coagu) ) + 
        a[["coefficients"]]["(Intercept)","Value"]
    }
  } else if (setequal(type, c("coagulation","erythrocyte")))  {
    ## coagulation erythrocyte ----
    for (i in 1:nrow(ndata1)) {
      y <- rawdata[i, ] # expression across samples
      y <- t(y)
      smpl2 <- smpl2[rownames(y),]
      x_eryth <- smpl2$erythrocyte
      x_coagu <- smpl2$coagulation
      Constraint_factor_eryth <- result_cor_erythrocyte[colnames(y),"avg"]
      Constraint_factor_coagu <- result_cor_coagulation[colnames(y),"avg"]
      a <- summary(rlm(log2(y + 1) ~ x_eryth + x_coagu, maxit = 30))
      for (j in 2:nrow(a$coefficients)) {
        b[i, j - 1] = a$coefficients[j, 1]
      }
      # new residuals:
      ndata1[i, ] <- log2(y + 1) - (a[["coefficients"]]["(Intercept)","Value"] + constraint * (
        x_eryth *a[["coefficients"]]["x_eryth","Value"]* Constraint_factor_eryth +
          x_coagu *a[["coefficients"]]["x_coagu","Value"]* Constraint_factor_coagu) ) + 
        a[["coefficients"]]["(Intercept)","Value"]
    }
  } else if (setequal(type, c("erythrocyte","platelet"))){
    ## erythrocyte platelet ----
    for (i in 1:nrow(ndata1)) {
      y <- rawdata[i, ] # expression across samples
      y <- t(y)
      smpl2 <- smpl2[rownames(y),]
      x_eryth <- smpl2$erythrocyte
      x_plate <- smpl2$Platelet
      Constraint_factor_eryth <- result_cor_erythrocyte[colnames(y),"avg"]
      Constraint_factor_plate <- result_cor_platelet[colnames(y),"avg"]
      a <- summary(rlm(log2(y + 1) ~ x_eryth + x_plate, maxit = 30))
      for (j in 2:nrow(a$coefficients)) {
        b[i, j - 1] = a$coefficients[j, 1]
      }
      # new residuals:
      ndata1[i, ] <- log2(y + 1) - (a[["coefficients"]]["(Intercept)","Value"] + constraint * (
        x_eryth *a[["coefficients"]]["x_eryth","Value"]* Constraint_factor_eryth + 
          x_plate * a[["coefficients"]]["x_plate","Value"]* Constraint_factor_plate) ) + 
        a[["coefficients"]]["(Intercept)","Value"]
    }
  } else if (identical(type, "erythrocyte")) {
    # 红细胞单独校正 ----
    for (i in 1:nrow(ndata1)) {
      y <- rawdata[i, ] # expression across samples
      y <- t(y)
      smpl2 <- smpl2[rownames(y),]
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
  } else if (identical(type, "platelet")) {
    # 血小板单独校正 ----
    for (i in 1:nrow(ndata1)) {
      y <- rawdata[i, ] # expression across samples
      y <- t(y)
      smpl2 <- smpl2[rownames(y),]
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
  } else if (identical(type, "coagulation"))  {
    # 凝血单独校正 ----
    for (i in 1:nrow(ndata1)) {
      y <- rawdata[i, ] # expression across samples
      y <- t(y)
      smpl2 <- smpl2[rownames(y),]
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
  } else {
    stop("Unsupported correction type combination: ", paste(type, collapse = ", "))
  }
  result <- 2 ^ ndata1 - 1
  data$correct_data <- result
  data$contamination_level <- smpl2
  return(data)
}
