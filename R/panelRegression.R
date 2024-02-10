#' @export
runRegressions <- function(
  dataFilePath,
  outFilePath,
  variableTable,
  regionName,
  timeName,
  timeRange = NULL,
  regionRegex = NULL,
  nBatches = 1,
  batchIndex = 1
) {

  data <- read_csv(dataFilePath)

  batch <- cerUtility::splitAndGetOneBatch(
    "regressionVariables",
    seq_len(nrow(variableTable)),
    nBatches,
    batchIndex)

  regressionResults <- runRegressionTable(
    data,
    variableTable[batch, ],
    regionName,
    timeName,
    timeRange,
    regionRegex
  )

  batchOutFilePath <-
    cerUtility::removeFileNameEnding(outFilePath) |>
    paste0("_", batchIndex, ".RDS")
  cerUtility::makeDirsIfNecessary(batchOutFilePath)
  saveRDS(
    lst(
      regressionResults,
      dataFilePath,
      outFilePath,
      variableTable,
      regionName,
      timeName,
      timeRange,
      regionRegex,
      nBatches,
      batchIndex),
    file = batchOutFilePath)
}


runRegressionTable <- function(
  data,
  variableTable,
  regionName,
  timeName,
  timeRange = NULL,
  regionRegex = NULL
) {

  if (hasValueString(regionRegex)) {
    data <- data |> filter(str_detect(!!sym(regionName), regionRegex))
  }
  if (length(timeRange) == 2) {
    data <- data |> filter(!!sym(timeName) >= timeRange[1], !!sym(timeName) <= timeRange[2])
  }

  resultList <- lapply(seq_len(nrow(variableTable)), \(i) {
    info <- lapply(variableTable, `[[`, i)
    result <- runRegressionCore(
      data,
      info$targetVariable,
      info$predictors,
      info$fixedEffects,
      regionName,
      timeName
    )
    result$name <- info$name
    return(result)
  })
  names(resultList) <- variableTable$name
  return(resultList)
}


runRegressionCore <- function(
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

  coeffs <- coef(fit)[predictors]
  nonNaPredictors <- predictors[!is.na(coeffs)]

  pt <- proc.time()
  cat("calculating variances... ")
  coefVcovList <- list(
      vcovPL = sandwich::vcovPL(
        fit,
        order.by = data[[timeName]],
        lag = 2)[nonNaPredictors, nonNaPredictors],
      vcovCLtime = sandwich::vcovCL(
        fit,
        cluster = data[[timeName]])[nonNaPredictors, nonNaPredictors],
      vcovCLregion = sandwich::vcovCL(
        fit,
        cluster = data[[regionName]])[nonNaPredictors, nonNaPredictors],
      vcovCLcountry = sandwich::vcovCL(
        fit,
        cluster = str_sub(data[[regionName]], 1, 3))[nonNaPredictors, nonNaPredictors]
  )
  cat("duration:", (proc.time()-pt)[3], "s\n")

  pt <- proc.time()
  cat("calculating influence... ")
  infl <- influence(fit)
  cat("duration:", (proc.time()-pt)[3], "s\n")

  result <- list(
    predictors = predictors, # p+
    nonNaPredictors = nonNaPredictors, # p
    coeffs = coeffs, # p
    coefVcovList = coefVcovList, # list(p x p)
    leverage = infl$hat, # n
    looSd = infl$sigma, # n
    influence = infl$coefficients[, nonNaPredictors] # n x p
  )

  return(result)
}

