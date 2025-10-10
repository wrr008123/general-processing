# 批量裁剪栅格

library(terra)

# 1. 设置路径
tif_folder <- "E:/0_data/Annual_CO2_concentration_1km_2003-2023"   # 例如 "D:/data/tifs"
shp_path <- "E:/0_data/Map/杨-mask84-中国边界/mask_84.shp"       # 例如 "D:/data/clip_area.shp"
output_folder <- "E:/0_data_process/China_CO2_2003_2023" # 例如 "D:/data/clipped"

# 确保输出文件夹存在
if(!dir.exists(output_folder)) dir.create(output_folder)

# 2. 读取shp
clip_shp <- vect(shp_path)

# 3. 批量读取tif文件
tif_files <- list.files(tif_folder, pattern = "\\.tif$", full.names = TRUE)

# 4. 循环裁剪并保存
for (tif in tif_files) {
  # 读取栅格
  r <- rast(tif)
  
  # 裁剪
  r_clipped <- crop(r, clip_shp) |> mask(clip_shp)
  
  # 生成输出文件名
  out_name <- file.path(output_folder, basename(tif))
  
  # 保存结果
  writeRaster(r_clipped, out_name, overwrite = TRUE)
  
  cat("已裁剪并保存:", out_name, "\n")
}

cat("所有tif文件裁剪完成！\n")
