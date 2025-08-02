library(terra)
library(zoo)

# 设置输入输出路径
input_dir <- "F:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/基础数据_对齐/SFR"
output_dir <- "F:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/基础数据_对齐/SFR_时间插值"


input_path <- input_dir
output_path <- output_dir

# 檢查輸出路徑是否存在，如果不存在則建立它
if (!dir.exists(output_path)) {
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  cat("Output directory created at:", output_path, "\n")
}

# 獲取輸入資料夾中所有的 .tif 檔案
# 使用 sort() 確保檔案按名稱順序排列，這代表了時間順序
cat("Searching for TIFF files in:", input_path, "\n")
tif_files <- sort(list.files(path = input_path, 
                             pattern = "\\.tif$", 
                             full.names = TRUE))

print(tif_files)

if (length(tif_files) == 0) {
  stop("No .tif files found in the specified input directory.")
} else {
  cat("Found", length(tif_files), "files to process.\n")
}


# --- 4. 讀取TIF檔案並堆疊為一個 SpatRaster 物件 ---

# terra::rast() 可以將多個檔案直接讀取為一個多圖層的 SpatRaster 物件
# 每個圖層代表一個時間點
cat("Loading and stacking raster files...\n")
raster_stack <- rast(tif_files)


# --- 5. 定義插值函數 ---

# 這個函數將會被應用到每一個像元的時間序列上
# x 是一個代表單一像元在所有時間點上的值的向量 (e.g., c(12, 15, NA, 20))
interpolation_function <- function(x) {
  # 使用 zoo::na.approx 進行線性插值
  # rule = 2 的作用是：如果序列的開頭或結尾是NA，則使用最近的有效值來填充
  # 這對於填補序列開始或結束的缺失很有用
  interpolated_values <- zoo::na.approx(x, rule = 2, na.rm = FALSE)
  return(interpolated_values)
}


# --- 6. 執行像元級插值 ---

# terra::app() 會將指定的函數應用於 SpatRaster 的每一個像元上
cat("Performing pixel-wise temporal interpolation... (This may take a while)\n")
interpolated_stack <- app(raster_stack, fun = interpolation_function)


# --- 7. 儲存插值後的結果 ---

# 產生輸出的檔案路徑與名稱，與輸入檔名保持一致
output_filenames <- file.path(output_path, basename(tif_files))

# 將插值後的多層 SpatRaster 物件寫入到多個單獨的檔案中
# terra::writeRaster 可以自動處理這個過程
cat("Saving interpolated rasters to:", output_path, "\n")
writeRaster(interpolated_stack, 
            filename = output_filenames, 
            overwrite = TRUE) # 如果檔案已存在，則覆蓋

cat("Processing complete!\n")
