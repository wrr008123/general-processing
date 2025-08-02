library(terra)

# 定义函数
align_rasters <- function(template_path, input_dir, output_dir, method = "bilinear") {
  # 读取模板栅格
  template <- rast(template_path)
  
  # 获取输入路径下所有.tif文件
  raster_files <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)
  
  # 创建输出目录（如不存在）
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 遍历每个栅格文件
  for (raster_file in raster_files) {
    cat("Processing:", raster_file, "\n")
    
    # 读取原始栅格
    r <- rast(raster_file)
    
    # 对齐到模板（resample）
    r_aligned <- resample(r, template, method = method)
    
    # 输出文件名
    output_filename <- file.path(output_dir, basename(raster_file))
    
    # 写入对齐后的栅格
    writeRaster(r_aligned, output_filename, overwrite = TRUE)
  }
  
  cat("All rasters aligned and saved to:", output_dir, "\n")
}



align_rasters(
  template_path = 'F:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/基础数据_对齐/SFR_时间插值/karst_SFR_2000.tif',
  input_dir = "F:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/基础数据_对齐/SFR_多年均值",
  output_dir = 'F:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/基础数据_对齐/SFR_多年均值对齐',
  method = "bilinear"  # 或 "near"（最近邻），适用于分类数据
)

