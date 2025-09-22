library(terra)
library(relaimpo)

# 这里不用改，使用的时候改下面【使用部分】
compute_lmg_rasters <- function(
    years = 2000:2020,
    response = list(dir = "./resp", pattern = "resp_%d.tif"),
    predictors = list(
      list(name = "NDVI",  dir = "./pred1", pattern = "NDVI_%d.tif"),
      list(name = "Temp",  dir = "./pred2", pattern = "Temp_%d.tif")
    ),
    output_dir = "lmg_outputs",
    file_suffix = ".tif",
    relative = TRUE,   # 新增参数：是否输出百分比贡献
    overwrite = FALSE,
    cores = 1
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 辅助函数：根据目录和模板解析文件列表
  resolve_files <- function(dir, pattern) {
    vapply(years, function(y) file.path(dir, sprintf(pattern, y)), FUN.VALUE = "")
  }
  
  # 读取响应栅格
  resp_files <- resolve_files(response$dir, response$pattern)
  if(!all(file.exists(resp_files))) {
    stop("响应文件缺失：", paste(resp_files[!file.exists(resp_files)], collapse=", "))
  }
  message("响应变量文件：")
  print(resp_files)
  resp_stack <- rast(resp_files)
  
  # 读取预测因子栅格
  pred_stacks <- list()
  pred_names <- character()
  for(i in seq_along(predictors)) {
    p <- predictors[[i]]
    if(is.null(p$name)) stop(sprintf("第 %d 个预测因子缺少 name 字段", i))
    pred_names[i] <- p$name
    fns <- resolve_files(p$dir, p$pattern)
    if(!all(file.exists(fns))) {
      stop(sprintf("预测因子 %s 文件缺失：%s", p$name,
                   paste(fns[!file.exists(fns)], collapse=", ")))
    }
    message(sprintf("预测因子 [%s] 文件：", p$name))
    print(fns)
    pred_stacks[[i]] <- rast(fns)
    if(!compareGeom(resp_stack, pred_stacks[[i]], stopOnError = FALSE)) {
      stop(sprintf("预测因子 %s 的几何与响应不一致，请先重采样/对齐。", p$name))
    }
  }
  
  npreds <- length(pred_stacks)
  ny <- length(years)
  
  # 合并所有栅格
  all_layers <- c(resp_stack, do.call(c, pred_stacks))
  
  pixel_fun <- function(vals) {
    if(all(is.na(vals))) return(rep(NA_real_, npreds))
    resp_v <- vals[1:ny]
    preds_v <- matrix(vals[(ny+1):length(vals)], ncol = npreds, byrow = FALSE)
    ok_idx <- which(!is.na(resp_v) & apply(!is.na(preds_v), 1, all))
    if(length(ok_idx) < (npreds + 1)) return(rep(NA_real_, npreds))
    df <- data.frame(y = resp_v[ok_idx])
    for(j in seq_len(npreds)) df[[paste0("x", j)]] <- preds_v[ok_idx, j]
    tryCatch({
      lm0 <- lm(y ~ ., data = df)
      re <- calc.relimp(lm0, type = "lmg", rela = FALSE)
      out <- re$lmg
      v <- sapply(paste0("x", seq_len(npreds)), function(nm) if(nm %in% names(out)) out[nm] else NA_real_)
      
      if(relative) {
        s <- sum(v, na.rm = TRUE)
        if(!is.na(s) && s > 0) v <- v / s  # 转为比例
      }
      v
    }, error = function(e) rep(NA_real_, npreds))
  }
  
  if(cores > 1) terraOptions(nthreads = cores)
  
  out <- app(all_layers, fun = pixel_fun, cores = cores)
  names(out) <- paste0("LMG_", pred_names)
  
  # 输出结果，每个因子单独一个 tif
  out_paths <- character(npreds)
  for(i in seq_len(npreds)) {
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", pred_names[i])  # 防止中文或特殊字符导致文件名问题
    out_file <- file.path(output_dir, paste0("LMG_contribution_", safe_name, file_suffix))
    writeRaster(out[[i]], filename = out_file, overwrite = overwrite)
    out_paths[i] <- normalizePath(out_file)
  }
  
  message("输出完成：", paste(out_paths, collapse=", "))
  invisible(out_paths)
}

#---------------注意-------------------
#
# 很多人第一次用 LMG 时的常见疑惑，就是为什么结果加起来不等于1
# LMG贡献之和等于模型的R2,是R2的分解,而不是 1,
# 可以设置relative参数为true，即可输出相对百分比（总和为1）
# %d 为年份的占位符
#
#---------------使用-------------------
paths <- compute_lmg_rasters(
  years = 2000:2020,
  response = list(dir = "F:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/基础数据_对齐/SFR_时间插值", pattern = "karst_SFR_%d.tif"),
  predictors = list(
    list(name = "T", dir = "F:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/基础数据_对齐/TMP", pattern = "TMP_%d.tif"),
    list(name = "P", dir = "F:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/基础数据_对齐/P_总量", pattern = "P_%d.tif"),
    list(name = "ET", dir = "F:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/基础数据_对齐/ET", pattern = "ET_%d.tif"),
    list(name = "KRDI", dir = "F:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/03_基础数据_标准化/KRD", pattern = "KRD_%d.tif")
  ),
  output_dir = "F:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/03_基础数据_标准化/06_LMG重要性",
  relative = TRUE,
  cores = 1,
  overwrite = TRUE
)
