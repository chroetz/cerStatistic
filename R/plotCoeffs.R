saveCoeffPlots <- function(
  rdsFilePath,
  outFilePath,
  name = NULL,
  vcovName = NULL,
  sigLevel = 0.95,
  cumulative = FALSE,
  width=1920,
  height=1080,
  dpi=300
) {
  plt <- createCoeffPlots(rdsFilePath, name, vcovName, sigLevel, cumulative)
  cerUtility::makeDirsIfNecessary(outFilePath)
  ggsave(outFilePath, plt, width=width, height=height, units="px", dpi=dpi)
}


createCoeffPlots <- function(
  rdsFilePath,
  name = NULL,
  vcovNames = NULL,
  sigLevel = 0.95,
  cumulative = FALSE
) {

  rds <- readRDS(rdsFilePath)
  if (hasValueString(name)) {
    reg <- rds$regressionResults[[name]]
  } else {
    reg <- rds$regressionResults[[1]]
    name <- reg$predictors[1]
  }
  vcovList <- if (hasValueString(vcovNames)) reg$coefVcovList[vcovNames] else reg$coefVcovList
  prefix <- cerUtility::longestCommonPrefix(reg$predictors)
  label <- c(reg$predictors[1], str_sub(reg$predictors[-1], nchar(prefix) + 1))
  lag <- reg$predictors |> str_extract("_lag[0-9]+$") |> str_remove("_lag") |> as.integer()
  lag[is.na(lag)] <- 0

  a <- 1-(1-sigLevel)/2
  q <- qnorm(a)

  if (cumulative) {
    coefVarList <-
      lapply(vcovList, \(vcov) vapply(seq_along(reg$coeffs), \(i) sum(vcov[1:i, 1:i]), numeric(1)))
    names(coefVarList) <- names(vcovList)
  } else {
    coefVarList <- lapply(vcovList, diag)
    names(coefVarList) <- names(vcovList)
  }

  pltDataCoef <- tibble(
    label = label,
    lag = lag,
    coef = if (cumulative) cumsum(reg$coeffs) else reg$coeffs
  )
  pltDataVar <-
    pltDataCoef |>
    bind_cols(coefVarList) |>
    pivot_longer(cols = -c(label, lag, coef), names_to = "vcov", values_to = "var") |>
    mutate(
      confiRadius = q * sqrt(var),
      lower = coef - confiRadius,
      upper = coef + confiRadius,
      signif = upper < 0)
  plt <-
    pltDataCoef |>
    ggplot(aes(x = lag, y = coef))
    plt <-
    plt +
    geom_ribbon(data = pltDataVar, aes(ymin = lower, ymax = upper, fill = vcov, color = vcov), alpha = 0.1)
  plt <- plt +
    geom_hline(yintercept=0) +
    geom_line(linewidth=1.0) +
    geom_point(size=2) +
    ggtitle(
      sprintf("%s, %g%% confidence", name, sigLevel*100),
      subtitle = sprintf("n = %d, rank = %d", rds$fit$residuals |> length(), rds$fit$rank)) +
    ylab(if (cumulative) "cumulative coefficient" else "coefficient")

  return(plt)
}
