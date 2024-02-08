fileRegression <- function(
  dataFilePath,
  outFilePath,
  targetVariable,
  predictors,
  fixedEffects,
  regionName,
  timeName
) {
  data <- read_csv(dataFilePath)
  result <- runRegression(
    data,
    targetVariable,
    predictors,
    fixedEffects,
    regionName,
    timeName
  )
  writeRDS(result, file = outFilePath)
}


runRegression <- function(
  data,
  targetVariable,
  predictors,
  fixedEffects,
  regionName,
  timeName
) {

  frml <- as.formula(
    paste0(
      targetVariable,
      " ~ ",
      paste0(c(predictors, fixedEffects), collapse = " + ")))

  pt <- proc.time()
  cat("running linear regression... ")
  fit <- lm(frml, data = data)
  cat("duration:", (proc.time()-pt)[3], "s\n")

  pt <- proc.time()
  cat("calculating variances... ")
  vcovList <-
    list(
      vcovPL = sandwich::vcovPL(fit, cluster = as.formula(paste("~", regionName, "+", timeName)))
    )
  cat("duration:", (proc.time()-pt)[3], "s\n")

  pt <- proc.time()
  cat("calculating influence... ")
  influence <- influence(fit)
  cat("duration:", (proc.time()-pt)[3], "s\n")

  coeffs <- coef(fit)[predictors]
  coefVcovList <- lapply(vcovList, \(v) v[predictors, predictors])

  result <- list(
    coeffs = coeffs, # p
    coefVcovList = coefVcovList, # list(p x p)
    leverage = infl$hat, # n
    looSd = infl$sigma, # n
    influence = infl$coefficients[, predictors], # n x p
  )

  return(result)
}
