library(terra)

# ===========================
# 函数：仅投影（不对齐分辨率/范围）
# ===========================
batch_project_crs_only <- function(template_tif, input_dir, output_dir, method = "bilinear") {
  
  # 1. 读取模板，提取坐标系
  template <- rast(template_tif)
  target_crs <- crs(template)
  cat("========== 模板TIF坐标系 ==========\n")
  cat(crs(template, proj = TRUE), "\n\n")  # 简洁格式
  
  
  # 2. 获取所有tif
  tif_files <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)
  
  # 3. 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 4. 循环处理
  for (i in seq_along(tif_files)) {
    
    cat("正在处理：", tif_files[i], "\n")
    
    
    r <- rast(tif_files[i])
    cat(crs(r, proj = TRUE), "\n\n")  # 简洁格式
    
    
    # ===========================
    # 核心：只按CRS投影
    # ===========================
    r_proj <- project(r, target_crs, method = method)
    
    # 输出路径
    out_path <- file.path(output_dir, basename(tif_files[i]))
    
    writeRaster(r_proj, out_path, overwrite = TRUE)
  }
  
  cat("全部完成！\n")
}

# ===========================
# 使用示例
# ===========================
template_tif <- "E:/0_paper1/SFR_China/SFR_align/SFR_2000.tif"
input_dir   <- "E:/0_data/MaxNDVI_ypr"
output_dir  <- "E:/0_data/MaxNDVI_ypr_TY"

batch_project_crs_only(template_tif=template_tif, input_dir=input_dir, output_dir=output_dir)
