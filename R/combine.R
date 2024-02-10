combineVariables <- function(
  outFilePath,
  lags,
  regionRegex = NULL,
  timeRange = NULL,
  targetRegionName = "GID_1",
  targetTimeName = "year",
  variableDescriptorList,
  requireNotNA = NULL,
  ...
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

  dataLag1 <- data
  for (variable in setdiff(names(data), c(targetRegionName, targetTimeName))) {
    dataLag1 <- addLagged(dataLag1, data, 1, variable)
  }

  dataDelta <-
    dataLag1 |>
    mutate(...) |>
    select(-ends_with("_lag1"))

  dataLagged <- dataDelta
  for (variable in setdiff(names(dataDelta), c(targetRegionName, targetTimeName))) {
    for (lag in lags) {
      dataLagged <- addLagged(dataLagged, dataDelta, lag, variable)
    }
  }

  dataLaggedClean <-
    dataLagged |>
    filter(IS_CLEAN) |>
    drop_na(all_of(requireNotNA))

  write_csv(dataLaggedClean, outFilePath)
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


addLagged <- function(x, y, lag, variable) {
    out <-
      x |>
      mutate(laggedYear = year - lag) |>
      left_join(
        y |>
          select(GID_1, year, all_of(variable)) |>
          rename_with(\(x) paste0(x, "_lag", lag), all_of(variable)),
        join_by(GID_1, laggedYear == year)) |>
      select(-laggedYear)
    return(out)
  }
