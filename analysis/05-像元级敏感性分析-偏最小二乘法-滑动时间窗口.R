library(terra)
library(pls)

run_pls_sensitivity <- function(
    dependent,
    independents,
    years,
    output_dir,
    ncomp = 1,
    scale = FALSE
) {
  
  #--------------------------------------------
  # 基本信息
  n_year <- length(years)
  all_names <- c(
    dependent$name,
    sapply(independents, function(x) x$name)
  )
  n_coef <- length(all_names) - 1
  
  #--------------------------------------------
  # 年份匹配读取函数
  read_yearly_stack <- function(dir, years) {
    files <- list.files(dir, pattern = "\\.tif$", full.names = TRUE)
    
    files_year <- sapply(years, function(y) {
      f <- files[grepl(as.character(y), basename(files))]
      
      message(sprintf(
        "Year %d -> %s",
        y, basename(f)
      ))
      
      if (length(f) != 1) return(NA)
      f
    })
    
    
    
    if (any(is.na(files_year))) {
      stop(paste("年份匹配失败：", dir))
    }
    
    rast(files_year)
  }
  
  #--------------------------------------------
  # 读取因变量
  y_stack <- read_yearly_stack(dependent$path, years)
  
  # 读取自变量
  x_stacks <- lapply(independents, function(v) {
    read_yearly_stack(v$path, years)
  })
  
  # 合并（y 在前）
  final_stack <- c(y_stack, rast(x_stacks))
  names(final_stack)
  
  #--------------------------------------------
  # 像元级 PLS 函数
  pls_pixel <- function(v, n_year, all_names, ncomp, scale) {
    
    n_coef <- length(all_names) - 1
    
    # 缺失值检查
    if (any(is.na(v))) {
      return(rep(NA_real_, n_coef))
    }
    
    # 拆分为时间 × 变量
    df <- matrix(v, nrow = n_year, byrow = FALSE)
    colnames(df) <- all_names
    
    Y <- as.matrix(df[, 1])
    X <- as.matrix(df[, -1])
    
    fit <- try(
      plsr(
        Y ~ X,
        ncomp = ncomp,
        validation = "none",
        scale = scale
      ),
      silent = TRUE
    )
    
    if (inherits(fit, "try-error")) {
      return(rep(NA_real_, n_coef))
    }
    
    as.vector(coef(fit, ncomp = ncomp))
  }
  
  #--------------------------------------------
  # 运行像元计算
  result_sens <- app(
    final_stack,
    pls_pixel,
    n_year   = n_year,
    all_names = all_names,
    ncomp    = ncomp,
    scale    = scale
  )
  
  names(result_sens) <- all_names[-1]
  
  #--------------------------------------------
  # 输出结果
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  start_year <- min(years)
  end_year   <- max(years)
  
  for (i in seq_len(n_coef)) {
    out_file <- file.path(
      output_dir,
      paste0(
        dependent$name, "_",
        all_names[-1][i], "_sens_",
        start_year, "_", end_year, ".tif"
      )
    )
    
    writeRaster(
      result_sens[[i]],
      filename = out_file,
      overwrite = TRUE
    )
  }
  
  return(result_sens)
}


run_pls_sensitivity_sliding <- function(
    dependent,
    independents,
    all_years,
    window,
    output_dir,
    ncomp = 1,
    scale = FALSE
) {
  
  start_years <- all_years[1:(length(all_years) - window + 1)]
  
  results <- list()
  
  for (sy in start_years) {
    
    ey <- sy + window - 1
    years_i <- sy:ey
    
    message("Running window: ", sy, "-", ey)
    
    res_i <- run_pls_sensitivity(
      dependent    = dependent,
      independents = independents,
      years        = years_i,
      output_dir   = output_dir,
      ncomp        = ncomp,
      scale        = scale
    )
    
    results[[paste0(sy, "_", ey)]] <- res_i
  }
  
  return(results)
}


# -------使用-----------------------------------------------
# 因变量：名称 + 路径
dependent <- list(
  name = "SFRHC",
  path = "G:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/02_基础数据_标准化/SFR_CI"
)

# 自变量（可增减多个）
independents <- list(
  list(name="KRD",   path="G:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/02_基础数据_标准化/KRD"),
  list(name="P", path="G:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/02_基础数据_标准化/P_总量"),
  list(name="T", path="G:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/02_基础数据_标准化/TMP"),
  list(name="ET", path="G:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/数据_v3/02_基础数据_标准化/ET")
)
# 数据年份（前后都包括的）
years <- 2000:2020
# 滑动时间窗口
window <- 10
# 输出路径
output_dir <- "G:/phd/03-第一篇论文/碳酸盐岩风化成土对石漠化恢复的响应_数据资料/像元敏感性分析/滑动窗口敏感性分析CI_10y_1comp"


# 运行滑动窗口 PLS
all_results <- run_pls_sensitivity_sliding(
  dependent    = dependent,
  independents = independents,
  all_years    = years,
  window       = window,
  output_dir   = output_dir,
  ncomp        = 1,
  scale        = FALSE
)


