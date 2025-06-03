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
  # 向量化计算约束因子
  compute_constraint_factors <- function(cor_matrix, marker_set, constraint_value) {
    if (is.null(marker_set) || length(marker_set) == 0) return(NULL)
    
    # 提取当前标记基因的相关系数子集
    sub_cor <- abs(cor_matrix[, marker_set, drop = FALSE])
    
    # 计算每个基因的约束因子
    constraint_factors <- apply(sub_cor, 1, function(x) {
      n <- length(x)
      if (n < 3) return(mean(x, na.rm = TRUE))
      
      # 计算起始位置
      start_idx <- max(1, min(round(n * constraint_value), n - 2))
      # 排序并取平均值
      sorted_vals <- sort(x)
      mean(sorted_vals[start_idx:(n - 1)], na.rm = TRUE)
    })
    
    return(constraint_factors)
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
  result_cor <- cor(as.matrix(t(data$rawdata)), method = "pearson")
  result_cor <- as.data.frame(result_cor)
  # 初始化约束因子为 NULL
  result_cor_erythrocyte <- result_cor_platelet <- result_cor_coagulation <- NULL
  # 计算约束因子
  result_cor_erythrocyte <- if (length(erythrocyte_marker) > 0) {
    compute_constraint_factors(cor_matrix = result_cor, 
                               erythrocyte_marker, 
                               constraint_erythrocyte_original)
  } else NULL
  
  result_cor_platelet <- if (length(platelet_marker) > 0) {
    compute_constraint_factors(cor_matrix = result_cor, 
                               platelet_marker, 
                               constraint_platelet_original)
  } else NULL
  
  result_cor_coagulation <- if (length(coagulation_marker) > 0) {
    compute_constraint_factors(cor_matrix = result_cor, 
                               coagulation_marker,
                               constraint_coagulation_original)
  } else NULL
  
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
      # 使用预计算的约束因子向量
      Constraint_factor_eryth <- if (!is.null(constraint_erythrocyte)) constraint_erythrocyte[colnames(y)] else 1
      Constraint_factor_plate <- if (!is.null(constraint_platelet)) constraint_platelet[colnames(y)] else 1
      Constraint_factor_coagu <- if (!is.null(constraint_coagulation)) constraint_coagulation[colnames(y)] else 1
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
      Constraint_factor_plate <- if (!is.null(constraint_platelet)) constraint_platelet[colnames(y)] else 1
      Constraint_factor_coagu <- if (!is.null(constraint_coagulation)) constraint_coagulation[colnames(y)] else 1
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
      Constraint_factor_eryth <- if (!is.null(constraint_erythrocyte)) constraint_erythrocyte[colnames(y)] else 1
      Constraint_factor_coagu <- if (!is.null(constraint_coagulation)) constraint_coagulation[colnames(y)] else 1
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
      Constraint_factor_eryth <- if (!is.null(constraint_erythrocyte)) constraint_erythrocyte[colnames(y)] else 1
      Constraint_factor_plate <- if (!is.null(constraint_platelet)) constraint_platelet[colnames(y)] else 1
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
      Constraint_factor_eryth <- if (!is.null(constraint_erythrocyte)) constraint_erythrocyte[colnames(y)] else 1
      a <- summary(rlm(log2(y + 1) ~ x_eryth, maxit = 30))
      for (j in 2:nrow(a$coefficients)) {
        b[i, j - 1] = a$coefficients[j, 1]
      }
      # new residuals:
      ndata1[i, ] <- log2(y + 1) - (a[["coefficients"]]["(Intercept)","Value"] + 
                                      constraint * x_eryth * a[["coefficients"]]["x_eryth","Value"] * Constraint_factor_eryth) + 
        a[["coefficients"]]["(Intercept)","Value"]
    }
  } else if (identical(type, "platelet")) {
    # 血小板单独校正 ----
    for (i in 1:nrow(ndata1)) {
      y <- rawdata[i, ] # expression across samples
      y <- t(y)
      smpl2 <- smpl2[rownames(y),]
      x_plate <- smpl2$Platelet
      Constraint_factor_plate <- if (!is.null(constraint_platelet)) constraint_platelet[colnames(y)] else 1
      # x_coagu <- smpl2$coagulation
      a <- summary(rlm(log2(y + 1) ~ x_plate, maxit = 30))
      for (j in 2:nrow(a$coefficients)) {
        b[i, j - 1] = a$coefficients[j, 1]
      }
      # new residuals:
      ndata1[i, ] <- log2(y + 1) - (a[["coefficients"]]["(Intercept)","Value"] + 
                                      constraint * x_plate * a[["coefficients"]]["x_plate","Value"] * Constraint_factor_plate) + 
        a[["coefficients"]]["(Intercept)","Value"]
    }
  } else if (identical(type, "coagulation"))  {
    # 凝血单独校正 ----
    for (i in 1:nrow(ndata1)) {
      y <- rawdata[i, ] # expression across samples
      y <- t(y)
      smpl2 <- smpl2[rownames(y),]
      x_coagu <- smpl2$coagulation
      Constraint_factor_coagu <- if (!is.null(constraint_coagulation)) constraint_coagulation[colnames(y)] else 1
      a <- summary(rlm(log2(y + 1) ~ x_coagu, maxit = 30))
      for (j in 2:nrow(a$coefficients)) {
        b[i, j - 1] = a$coefficients[j, 1]
      }
      # new residuals:
      ndata1[i, ] <- log2(y + 1) - (a[["coefficients"]]["(Intercept)","Value"] + 
                                      constraint * x_coagu * a[["coefficients"]]["x_coagu","Value"] * Constraint_factor_coagu) + 
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
