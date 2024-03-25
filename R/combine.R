#' @export
combineVariables <- function(
  outFilePath,
  lags,
  regionRegex = NULL,
  timeRange = NULL,
  targetRegionName = "GID_1",
  targetTimeName = "year",
  variableDescriptorList,
  requireNotNA = NULL,
  saveMerge = FALSE,
  ...
) {

  data <- mergeData(
    regionRegex,
    timeRange,
    targetRegionName,
    targetTimeName,
    variableDescriptorList
  )

  if (saveMerge) {
    write_csv(
      data,
      paste0(cerUtility::removeFileNameEnding(outFilePath), "_merge.csv"))
  }

  dataLaggedClean <- addDiffAndLag(
    data,
    targetRegionName,
    targetTimeName,
    lags,
    requireNotNA,
    ...)

  write_csv(dataLaggedClean, outFilePath)
}


#' @export
mergeData <- function(
    regionRegex,
    timeRange,
    targetRegionName,
    targetTimeName,
    variableDescriptorList
  ) {
  dataList <- lapply(
    variableDescriptorList,
    loadData,
    targetRegionName = targetRegionName,
    targetTimeName = targetTimeName,
    regionRegex = regionRegex,
    timeRange = timeRange)
  data <- Reduce(
    function(x, y) full_join(x, y, by = c(targetRegionName, targetTimeName)),
    dataList)
  return(data)
}


#' @export
addDiffAndLag <- function(
    data,
    targetRegionName,
    targetTimeName,
    lags,
    requireNotNA,
    ...
) {
  numericVars <- setdiff(names(data)[sapply(data, is.numeric)], c(targetRegionName, targetTimeName))

  for (variable in numericVars) {
    data[[paste0(variable, "_2")]] <- data[[variable]]^2
  }

  dataLag1 <- data
  for (variable in numericVars) {
    dataLag1 <- addLagged(dataLag1, data, 1, variable, targetRegionName, targetTimeName)
  }

  data <-
    dataLag1 |>
    mutate(...) |>
    select(-ends_with("_lag1"))

  numericVars <- setdiff(names(data)[sapply(data, is.numeric)], c(targetRegionName, targetTimeName))

  dataLag1 <- data
  for (variable in numericVars) {
    dataLag1 <- addLagged(dataLag1, data, 1, variable, targetRegionName, targetTimeName)
  }

  data <- dataLag1
  for (vn in numericVars) {
    data <- data |> mutate(!!paste0(vn, "_d") := !!sym(vn) - !!sym(paste0(vn, "_lag1")))
  }
  data <-
    data |>
    select(-ends_with("_lag1"))

  numericVars <- setdiff(names(data)[sapply(data, is.numeric)], c(targetRegionName, targetTimeName))

  dataLagged <- data
  for (variable in numericVars) {
    for (lag in lags) {
      dataLagged <- addLagged(dataLagged, data, lag, variable, targetRegionName, targetTimeName)
    }
  }

  data <-
    dataLagged |>
    drop_na(all_of(requireNotNA))

  return(data)
}


loadData <- function(
  info,
  targetRegionName, targetTimeName,
  regionRegex = NULL, timeRange = NULL
) {
  data <-
    read_csv(info$filePath) |>
    select(all_of(c(info$regionName, info$timeName, info$variableNames))) |>
    rename(
      !!targetRegionName := !!sym(info$regionName),
      !!targetTimeName := !!sym(info$timeName))
  if (hasValueString(regionRegex)) {
    data <- data |> filter(str_detect(!!sym(info$regionName), regionRegex))
  }
  if (length(timeRange) == 2) {
    data <- data |> filter(!!sym(info$timeName) >= timeRange[1], !!sym(info$timeName) <= timeRange[2])
  }
  return(data)
}


addLagged <- function(x, y, lag, variable, regionName, timeName) {
  out <-
    x |>
    mutate(laggedTime = !!sym(timeName) - lag) |>
    left_join(
      y |>
        select(!!sym(regionName), !!sym(timeName), all_of(variable)) |>
        rename_with(\(x) paste0(x, "_lag", lag), all_of(variable)),
      join_by(!!sym(regionName), laggedTime == !!sym(timeName))) |>
    select(-laggedTime)
  return(out)
}
