library(terra)
library(trend)

base_dir <- "G:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/像元敏感性分析/滑动窗口敏感性分析CI_10y_1comp_lag1/KRD_sens-abs"
out_dir <- 'G:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/像元敏感性分析/滑动窗口敏感性分析CI_10y_1comp_lag1/趋势分析'
prefix <- 'CI_KRD_sens'

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

flnames <- list.files(path = file.path(base_dir), pattern = '.tif$')
fl <- file.path(base_dir, flnames)
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

firs_sen <- app(firs, fun_sen, n_year=n_total_years)
names(firs_sen) <- c("Z", "slope", "p-value")

# 保存结果
writeRaster(firs_sen[[1]], filename = file.path(out_dir, paste0(prefix,"_Z.tif")), overwrite = TRUE)
writeRaster(firs_sen[[2]], filename = file.path(out_dir, paste0(prefix,"_slope.tif")), overwrite = TRUE)
writeRaster(firs_sen[[3]], filename = file.path(out_dir, paste0(prefix,"_pvalue.tif")), overwrite = TRUE)
