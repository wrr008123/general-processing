library(terra)
library(trend)
#输入一个文件夹内的单波段TIFF数据，在这里是历年的NDVI年最大值

output_name_prefix <- 'ET'
base_dir <- "F:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/03_基础数据_标准化/ET"
out_dir <- 'F:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/趋势分析/数据标准化后趋势结果'

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

flnames <- list.files(path=file.path(base_dir),pattern='.tif$')
fl<-file.path(base_dir,flnames)
print(fl)

firs<-rast(fl)

#Sen+MK计算
fun_sen<-function(x){
  if(length(na.omit(x))<20) return(c(NA, NA, NA))#删除数据不连续含有NA的像元
  MK_estimate <- trend::sens.slope(ts(na.omit(x), frequency = 1), conf.level = 0.95) #Sen斜率估计
  slope <- MK_estimate$estimates
  MK_test <- MK_estimate$p.value
  Zs <- MK_estimate$statistic
  return(c(Zs, slope, MK_test))
}

firs_sen = app(firs, fun_sen, cores=4 )
names(firs_sen) <- c("Z", "slope", "p-value")


# 分别写出三个结果
writeRaster(firs_sen[[1]], filename = file.path(out_dir, paste0(output_name_prefix, "_Z.tif")), overwrite=TRUE)
writeRaster(firs_sen[[2]], filename = file.path(out_dir, paste0(output_name_prefix,"_slope.tif")), overwrite=TRUE)
writeRaster(firs_sen[[3]], filename = file.path(out_dir, paste0(output_name_prefix,"_pvalue.tif")), overwrite=TRUE)