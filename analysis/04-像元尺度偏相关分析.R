library(terra)
library(ppcor)

#' @title 使用 terra 包进行像元尺度的偏相关分析
#'
#' @description 对多个栅格时间序列数据，在像元级别计算指定变量间的偏相关性。
#'
#' @param pattern_list list 一个列表，包含各个变量的栅格文件路径和匹配模式。
#'   格式: `list(var_name = list(path = "...", pattern = "..."))`
#' @param dependent_var character 因变量的名称 (必须是 pattern_list 中的一个键)。
#' @param independent_var character 自变量的名称 (必须是 pattern_list 中的一个键)。
#' @param control_vars character 向量，包含一个或多个控制变量的名称 (必须是 pattern_list 中的键)。
#' @param output_dir character 输出相关系数和p值栅格文件的文件夹。
#' @param cores integer 用于并行计算的核心数。默认为 1。
#'
#' @return 返回一个包含两层（相关系数和p-value）的 SpatRaster 对象。
#'
calc_partial_corr_raster_terra <- function(pattern_list,
                                           dependent_var,
                                           independent_var,
                                           control_vars,
                                           output_dir,
                                           cores = 1) {
  
  # ----------- 1. 加载和准备栅格数据 -----------
  
  cat("Step 1: Loading raster files...\n")
  
  # 将所有变量名整合到一个有序列表中，确保因变量和自变量在前
  all_vars <- c(dependent_var, independent_var, control_vars)
  
  raster_stacks <- list()
  n_layers <- NULL # 用于验证所有变量的时间序列长度是否一致
  
  for (var in all_vars) {
    if (!var %in% names(pattern_list)) {
      stop(paste("Variable '", var, "' not found in pattern_list."))
    }
    
    # 获取文件列表并创建 SpatRaster
    files <- list.files(
      path = pattern_list[[var]]$path,
      pattern = pattern_list[[var]]$pattern,
      full.names = TRUE
    )
    
    
    if (length(files) == 0) {
      stop(paste("No files found for variable '", var, "' with pattern '", pattern_list[[var]]$pattern, "' in path '", pattern_list[[var]]$path, "'"))
    }
    
    raster_stacks[[var]] <- rast(files)
    cat(paste("  - Loaded", length(files), "rasters for variable:", var, "\n"))
    
    for (f in files) {
      cat(basename(f), "\n")
    }
    
    # 验证图层数（时间序列长度）是否一致
    if (is.null(n_layers)) {
      n_layers <- nlyr(raster_stacks[[var]])
    } else if (nlyr(raster_stacks[[var]]) != n_layers) {
      stop("Mismatch in the number of layers (time steps) between variables.")
    }
  }
  
  # 将所有 SpatRaster 合并为一个多层对象
  # 顺序为：因变量, 自变量, 控制变量1, 控制变量2, ...
  combined_stack <- do.call(c, c(
    list(raster_stacks[[dependent_var]]),
    list(raster_stacks[[independent_var]]),
    raster_stacks[control_vars]
  ))
  
  cat("\nStep 2: All rasters stacked. Total layers:", nlyr(combined_stack), "\n")
  
  
  # ----------- 2. 定义像元级偏相关函数 -----------
  
  pcor_pixel_func <- function(x) {
    # x 是一个向量，包含了单个像元上所有变量在所有时间点的值
    
    # 如果像元值全部为 NA，则直接返回 NA
    if (all(is.na(x))) {
      return(c(NA, NA)) # 分别对应 coefficient 和 p-value
    }
    
    # 将一维向量 x 重塑为 [时间序列长度 x 变量数] 的矩阵
    # terra 的 app 函数按像元应用函数，x 的顺序与 combined_stack 的图层顺序一致
    num_vars <- length(all_vars)
    pixel_matrix <- matrix(x, ncol = num_vars, byrow = FALSE)
    
    # 转换为数据框并命名，方便索引
    df <- as.data.frame(pixel_matrix)
    colnames(df) <- all_vars
    
    # 检查数据有效性（例如，方差是否为零）
    # 如果因变量、自变量或任何控制变量没有变化，则无法计算相关性
    if (any(sapply(df, function(col) length(unique(na.omit(col))) < 2))) {
      return(c(NA, NA))
    }
    
    # 使用 tryCatch 来处理计算中可能出现的错误
    result <- tryCatch({
      # 执行偏相关检验
      pcor_res <- ppcor::pcor.test(
        x = df[[dependent_var]],
        y = df[[independent_var]],
        z = df[control_vars],
        method = "pearson" # 可按需更改为 "spearman"
      )
      
      # 返回相关系数和 p-value
      c(pcor_res$estimate, pcor_res$p.value)
      
    }, error = function(e) {
      # 如果计算出错，同样返回 NA
      c(NA, NA)
    })
    
    return(result)
  }
  
  
  # ----------- 3. 应用函数并保存结果 -----------
  
  cat(paste("\nStep 3: Starting pixel-wise partial correlation calculation on", cores, "core(s)...\n"))
  
  # 使用 terra::app 在每个像元上应用 pcor_pixel_func 函数
  # app 函数会自动处理并行化
  start_time <- Sys.time()
  
  pcor_result_raster <- app(
    combined_stack,
    fun = pcor_pixel_func,
    cores = cores
  )
  
  end_time <- Sys.time()
  cat("Calculation finished in:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes.\n")
  
  # 分离结果图层
  pcor_coef <- pcor_result_raster[[1]]
  pcor_pval <- pcor_result_raster[[2]]
  
  # 为图层命名
  names(pcor_coef) <- paste0("pcor_coef_", dependent_var, "_vs_", independent_var)
  names(pcor_pval) <- paste0("pcor_pval_", dependent_var, "_vs_", independent_var)
  
  cat("\nStep 4: Writing output rasters...\n")
  
  
  # 如果目录不存在，则分别创建
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  output_coef_file = file.path(output_dir,paste0('pcorr_',dependent_var,'_',independent_var,'_coef.tif'))
  output_pval_file = file.path(output_dir,paste0('pcorr_',dependent_var,'_',independent_var,'_pval.tif'))
  
  # 保存输出文件
  writeRaster(pcor_coef, output_coef_file, overwrite = TRUE)
  writeRaster(pcor_pval, output_pval_file, overwrite = TRUE)
  
  cat(paste("  - Coefficient raster saved to:", output_coef_file, "\n"))
  cat(paste("  - P-value raster saved to:", output_pval_file, "\n"))
  
  # 返回包含两个图层的 SpatRaster 对象
  return(invisible(c(pcor_coef, pcor_pval)))
}


# -----------使用------------
# 可以接着添加多个变量
patterns <- list(
  sfr = list(path = "F:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/基础数据_对齐/SFR_时间插值", pattern = "karst_SFR_\\d+\\.tif$"),
  krd = list(path = "F:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/03_基础数据_标准化/KRD", pattern = "KRD_\\d+\\.tif$"),
  p = list(path = "F:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/基础数据_对齐/P_总量", pattern = "P_\\d+\\.tif$"),
  tmp = list(path = "F:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/基础数据_对齐/TMP", pattern = "TMP_\\d+\\.tif$"),
  et = list(path = "F:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/基础数据_对齐/ET", pattern = "ET_\\d+\\.tif$")
)

# 因变量
dependent_var = "sfr"
# 自变量
independent_var = "et"
# 控制变量
control_vars = c('tmp','krd','p')

# 还得是一个核心最快 cores=1最快
result <- calc_partial_corr_raster_terra(
  pattern_list = patterns,
  dependent_var = dependent_var,
  independent_var = independent_var,
  control_vars = control_vars,
  output_dir = "F:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/偏相关分析_new0729",
  cores = 1
)
