#############################################
# 逐像元 Granger 因果分析
#############################################

#############################################
# 参数区（只改这里）
#############################################

X_path <- "E:/download/edge_download/NPP/CO2"  
Y_path <- "E:/download/edge_download/NPP/NPP"   
output_path <- "E:/download/edge_download/NPP/结果-lag1"      

# 输出结果前缀
X_prefix <- "CO2"
Y_prefix <- "NPP"

# 滞后系数
lag      <- 1
p_cutoff <- 0.05

# 并行运行线程数
cores <- 12


#############################################
# 1. 加载包
#############################################
library(terra)
library(tseries)
library(lmtest)

# 检查目录是否存在
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

#############################################
# 2. 读取数据
#############################################
X_files <- list.files(X_path, pattern="tif$", full.names=TRUE)
Y_files <- list.files(Y_path, pattern="tif$", full.names=TRUE)
print(X_files)
print(Y_files)
stopifnot(length(X_files) == length(Y_files))

X <- rast(X_files)
Y <- rast(Y_files)

n_total_years <- nlyr(X)


#############################################
# 3. 像元函数（修正版本）
#############################################
pixel_granger <- function(v,n_year,lag){
  
  # 拆分
  x_ts <- v[1:n_year]
  y_ts <- v[(n_year+1):(2*n_year)]
  
  # NA 检查
  if (any(is.na(x_ts)) || any(is.na(y_ts))) {
    return(c(NA, NA))
  }
  
  # ---- 安全 ADF 函数 ----
  safe_adf <- function(ts){
    p <- try(tseries::adf.test(ts)$p.value, silent = FALSE)
    if (inherits(p, "try-error")) return(NA_real_)
    p
  }
  
  # ---- 在函数里调用 ----
  px <- safe_adf(x_ts)
  py <- safe_adf(y_ts)
  
  
  # 平稳性 & 差分
  if (!is.na(px) && px > 0.05) x_ts <- diff(x_ts)
  if (!is.na(py) && py > 0.05) y_ts <- diff(y_ts)
  
  n <- min(length(x_ts), length(y_ts))
  if (n < 10) return(c(NA, NA))
  
  x <- x_ts[1:n]
  y <- y_ts[1:n]
  df <- data.frame(x, y)
  
  # ============= 安全取 p 值函数（解决>1问题） =============
  safe_p <- function(res){
    if (inherits(res,"try-error")) return(NA)
    pvals <- res$`Pr(>F)`
    if (length(pvals) >= 2) return(pvals[2])
    NA
  }
  
  # X -> Y
  g1 <- try(lmtest::grangertest(y ~ x, data=df, order=lag), silent=FALSE)
  p1 <- safe_p(g1)
  
  # Y -> X
  g2 <- try(lmtest::grangertest(x ~ y, data=df, order=lag), silent=FALSE)
  p2 <- safe_p(g2)
  
  return(c(p1, p2))
}


#############################################
# 4. 逐像元执行（关键修复：用 lapp()）
#############################################
both <- c(X, Y)

result <- app(both, pixel_granger,n_year=n_total_years,lag=lag,cores=cores)

names(result) <- c(
  paste0("p_", X_prefix, "_to_", Y_prefix),
  paste0("p_", Y_prefix, "_to_", X_prefix)
)


#############################################
# 5. 因果方向分类
#############################################
r1 <- result[[1]] < p_cutoff
r2 <- result[[2]] < p_cutoff

direction <- r1*1 + r2*2

# 0: 无
# 1: X→Y
# 2: Y→X
# 3: 双向


#############################################
# 6. 输出
#############################################

name1 <- paste0("Granger_p_", X_prefix, "_to_", Y_prefix, "_lag",lag,".tif")
name2 <- paste0("Granger_p_", Y_prefix, "_to_", X_prefix, "_lag",lag, ".tif")

writeRaster(result[[1]],
            file.path(output_path, name1),
            overwrite=TRUE)

writeRaster(result[[2]],
            file.path(output_path, name2),
            overwrite=TRUE)

direction_name <- paste0("Granger_direction_", X_prefix, "_", Y_prefix,"_lag",lag,  ".tif")
writeRaster(direction,
            file.path(output_path, direction_name),
            overwrite=TRUE)


message("==================================")
message("Granger 因果分析完成）")
message("输出：")
message(paste0(output_path, "/",name1))
message(paste0(output_path, "/",name2))
message(paste0(output_path, "/",direction_name))
message("==================================")
