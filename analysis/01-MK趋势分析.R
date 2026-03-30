library(terra)
library(trend)

input_dir <- "E:/Administrator/Documents/wechat_data/xwechat_files/wxid_bsczvtwojnmv22_bd93/msg/file/2026-03/01/01"
out_dir <- 'E:/test'
prefix <- 'spei_daily'

# 并行运行线程数
cores <- 16

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

flnames <- list.files(path = file.path(input_dir), pattern = '.tif$')
fl <- file.path(input_dir, flnames)
print(fl)

firs <- rast(fl)
n_total_years <- length(fl)  # 总年份数


cat("总年份数:", n_total_years, "\n")


# Sen+MK计算 - 通用版本
fun_sen <- function(x,n_year) {
  valid_data <- na.omit(x)
  
  # 动态判断数据完整性 这里可能是需要修改的
  if (length(valid_data) < n_year) {
    return(c(NA, NA, NA))
  }
  
  MK_estimate <- trend::sens.slope(valid_data, conf.level = 0.95)
  return(c(MK_estimate$statistic, MK_estimate$estimates, MK_estimate$p.value))
}

firs_sen <- app(firs, fun_sen, n_year=n_total_years,cores=cores)
names(firs_sen) <- c("Z", "slope", "p-value")

# 保存结果
writeRaster(firs_sen[[1]], filename = file.path(out_dir, paste0(prefix,"_Z.tif")), overwrite = TRUE)
writeRaster(firs_sen[[2]], filename = file.path(out_dir, paste0(prefix,"_slope.tif")), overwrite = TRUE)
writeRaster(firs_sen[[3]], filename = file.path(out_dir, paste0(prefix,"_pvalue.tif")), overwrite = TRUE)
