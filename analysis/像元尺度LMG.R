library(terra)
library(relaimpo)


# 设置工作目录
setwd("F:/project/wrr/data/aligned/03-对齐数据")

# 读取数据
et <- rast("ET_aligned.tif")
ndvi <- rast("NDVI_aligned.tif")
p <- rast("P_aligned.tif")
tmp <- rast("TMP_aligned.tif")
sfr <- rast("SFR_aligned.tif")

# 校验
if (!compareGeom(et, ndvi, p, tmp, sfr)) stop("Input rasters are not aligned!")
if (nlyr(et) != 21) stop("Each raster must have 21 bands!")

# 合并栅格（105 个图层: 5 个变量 × 21 年）
all_vars <- c(sfr, et, ndvi, p, tmp)

# 定义函数：单像素回归分析，返回各变量相对贡献率
pixel_lmg <- function(values) {
  # 拆分每种变量的 21 年数据
  sfr_vals <- values[1:21]
  et_vals <- values[22:42]
  ndvi_vals <- values[43:63]
  p_vals <- values[64:84]
  tmp_vals <- values[85:105]
  
  # 有缺失就返回 NA
  if (any(is.na(values))) return(rep(NA, 4))
  
  # 构建数据框并拟合模型
  df <- data.frame(
    SFR = sfr_vals,
    ET = et_vals,
    NDVI = ndvi_vals,
    P = p_vals,
    TMP = tmp_vals
  )
  lm_model <- tryCatch(lm(SFR ~ ET + NDVI + P + TMP, data = df), error = function(e) return(NULL))
  if (is.null(lm_model)) return(rep(NA, 4))
  
  relimp <- calc.relimp(lm_model, type = "lmg", rela = TRUE)
  c(relimp$lmg["ET"], relimp$lmg["NDVI"], relimp$lmg["P"], relimp$lmg["TMP"])
}

# 应用函数到整个栅格栈，输出为 4 层栅格（每层一个变量的贡献率）
contrib_stack <- app(all_vars, pixel_lmg)

# 设置图层名称
names(contrib_stack) <- c("ET_contrib", "NDVI_contrib", "P_contrib", "TMP_contrib")

# 保存结果
writeRaster(contrib_stack, 
            filename = c("ET_contribution.tif", "NDVI_contribution.tif", "P_contribution.tif", "TMP_contribution.tif"), 
            overwrite = TRUE)

cat("Finished: All contribution maps saved.\n")
