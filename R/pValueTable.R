getPValueTable <- function(data, predictors, target, timeTrendOrder, cumulate) {

  if (length(predictors) == 0) {
    return(NULL)
  }

  frml <- formula(paste(
      target, "~",
      paste(predictors, collapse = " + "),
      "- 1"))
  fit <- lm(frml, data = data)

  coefFun <- buildCoefFun(cumulate = cumulate)
  coeffs <- coefFun(fit, predictors)
  correctionFactor <- getFeCorrectionFactor(data, predictors, timeTrendOrder)

  vcovList <- list(
    homoscedastic = sandwich::vcovHC(fit, type = "const"),
    heteroscedastic = sandwich::vcovHC(fit, type = "HC1"),
    clusterRegion = sandwich::vcovCL(fit, cluster = data$ID),
    clusterYear = sandwich::vcovCL(fit, cluster = data$theYear),
    DriscollKraay0 = sandwich::vcovPL(fit, order.by = data$theYear, lag = 0, aggregate = TRUE),
    DriscollKraay1 = sandwich::vcovPL(fit, order.by = data$theYear, lag = 1, aggregate = TRUE),
    DriscollKraay2 = sandwich::vcovPL(fit, order.by = data$theYear, lag = 2, aggregate = TRUE),
    DriscollKraay3 = sandwich::vcovPL(fit, order.by = data$theYear, lag = 3, aggregate = TRUE)
  )
  clusterList <- list(
    clusterRegion = data$ID,
    clusterYear = data$theYear)
  numberOfClusters <- sapply(clusterList, \(x) x |> unique() |> length())

  se <- vapply(
    vcovList,
    \(vcov) cumulateSe(vcov, correctionFactor, cumulate),
    double(length(coeffs)))

  pValues <- 2 * (1 - pnorm(abs(coeffs), sd=se))
  pValues <- matrix(pValues, nrow = length(vcovList), byrow = TRUE, dimnames = list(names(vcovList), names(coeffs)))

  pValuesTable <-
    pValues |>
    as_tibble(rownames="type") |>
    left_join(tibble(type = names(numberOfClusters), numberOfClusters = numberOfClusters), join_by(type)) |>
    select(type, numberOfClusters, everything()) |>
    rename("#clusters" = numberOfClusters)
  pValuesTableWithCoef <-
    bind_rows(
      c(list(type = "coefficient", "#clusters" = NA), as.list(coeffs)),
      pValuesTable)
  pValuesTableColorValues <-
    bind_rows(
      list(type = "coefficient"),
      pValuesTable)
  pvTable <-
    pValuesTableWithCoef |>
    makeTable("p-values") |>
    colorizeTable(3:ncol(pValuesTableWithCoef), pValuesTableColorValues)
  return(pvTable)
}
