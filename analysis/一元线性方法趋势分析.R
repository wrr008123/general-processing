library(terra)

analyze_trend_terra <- function(input_dir, output_dir, output_prefix) {
  
  # 自动创建输出目录（如果不存在）
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 获取匹配前缀的所有tif文件
  tif_files <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)
  
  # 按文件名排序，确保顺序一致
  tif_files <- sort(tif_files)
  
  if (length(tif_files) < 2) {
    stop("❌ 至少需要两个 .tif 文件用于趋势分析。")
  }
  
  # 读取栅格栈
  s <- rast(tif_files)
  n <- nlyr(s)  # 图层数（时间点数）
  
  # 定义线性回归函数，使用时间序列索引 (1, 2, ..., n)
  linear_reg <- function(x) {
    if (all(is.na(x)) || sum(!is.na(x)) < 2) return(c(NA, NA))
    time_index <- which(!is.na(x))
    y <- x[time_index]
    t <- time_index
    fit <- lm(y ~ t)
    slope <- coef(fit)[2]
    p_val <- summary(fit)$coefficients[2, 4]
    return(c(slope, p_val))
  }
  
  # 应用回归分析（返回2层：斜率和p值）
  result <- app(s, fun = linear_reg)
  names(result) <- c("slope", "p_value")
  
  # 输出路径拼接
  slope_path <- file.path(output_dir, paste0(output_prefix, "_slope.tif"))
  pval_path  <- file.path(output_dir, paste0(output_prefix, "_p.tif"))
  sig_path   <- file.path(output_dir, paste0(output_prefix, "_sig0.5_slope.tif"))
  
  # 保存趋势和p值图
  writeRaster(result[[1]], slope_path, overwrite = TRUE)
  writeRaster(result[[2]], pval_path, overwrite = TRUE)
  
  # 显著性筛选
  slope_sig <- mask(result[[1]], result[[2]] <= 0.05, maskvalues = FALSE)
  writeRaster(slope_sig, sig_path, overwrite = TRUE)
  
  cat("✅ 趋势分析完成，结果已保存到：\n",
      slope_path, "\n", pval_path, "\n", sig_path, "\n")
}

#-------------使用--------------------
analyze_trend_terra(
  input_dir = "E:/China/p-e/",
  output_dir = "E:/China/p-e/",
  output_prefix = "pe"
)

