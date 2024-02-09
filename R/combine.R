combineVariables <- function(
  lags = 1:10,
  regionRegex = NULL,
  timeRange = c(1950, 2020),
  variablePath = "~/cerROOT/variable",
  targetRegionName = "GID_1",
  targetTimeName = "year",
  variableDescriptorList = list(
    dose = list(
      filePath = file.path(variablePath, "DOSE/DOSE_V2.csv"),
      regionName = "GID_1",
      timeName = "year",
      variableNames = c(grppc = "grp_pc_usd_2015", "StructChange")),
    storm = list(
      filePath = file.path(variablePath, "storm/storm.csv"),
      regionName = "GID_1",
      timeName = "year",
      variableNames = c(storm = "thresh34kt")),
    drought = list(
      filePath = file.path(variablePath, "drought/drought.csv"),
      regionName = "region",
      timeName = "year",
      variableNames = c(drought = "percentile50")),
    flood = list(
      filePath = file.path(variablePath, "flood/flood_flopros.csv"),
      regionName = "GID_1",
      timeName = "year",
      variableNames = c(flood = "shareOfPeopleAffectedByFlddphGt0")),
    tas = list(
      filePath = file.path(variablePath, "gswp3-w5e5/tasDoseYear_pop_1950-2019.csv"),
      regionName = "iso",
      timeName = "year",
      variableNames = c(tas = "tasMean_pop")),
    pr = list(
      filePath = file.path(variablePath, "gswp3-w5e5/prDoseYear_pop_1950-2019.csv"),
      regionName = "iso",
      timeName = "year",
      variableNames = c(pr = "prMean_pop")),
    pop = list(
      filePath = file.path(variablePath, "HYDE/HYDE3p3_population_DoseArcmin5.csv"),
      regionName = "GID_1",
      timeName = "year",
      variableNames = c(pop = "population")))
) {

  dataList <- lapply(
    infos,
    loadData,
    targetRegionName = targetRegionName,
    targetTimeName = targetTimeName,
    regionRegex = regionRegex,
    timeRange = timeRange)
  data <- Reduce(
    function(x, y) full_join(x, y, by = c(targetRegionName, targetTimeName)),
    dataList)
  prevData <-
    data |>
    rename_with(\(x) paste0("prev_", x), where(is.numeric))
  dataDelta <-
    data |>
    mutate(prev_year = year - 1) |>
    left_join(prevData, by=c(targetRegionName, "prev_year")) |>
    mutate( # TODO: generalize
      log2growth = log2(grppc/prev_grppc),
      dstorm = storm - prev_storm,
      ddrought = drought - prev_drought,
      dflood = flood - prev_flood,
      dtas = tas - prev_tas,
      dpr = pr - prev_pr,
      dpop = pop - prev_pop) |>
    select(-starts_with("prev_"))
  write_csv(dataDelta, file.path(variablePath, "dataDelta.csv"))

  dataLagged <- dataDelta
  for (variable in setdiff(names(dataDelta), c(targetRegionName, targetTimeName))) {
    for (lag in lags) {
      dataLagged <- addLagged(dataLagged, dataDelta, lag, variable)
    }
  }
  write_csv(dataLagged, file.path(variablePath, "dataLagged.csv"))

  dataLaggedClean <-
    dataLagged |>
    filter(StructChange == 0) |>
    drop_na(log2growth)
  write_csv(dataLaggedClean, file.path(variablePath, "dataLaggedClean.csv"))
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
          rename_with(\(x) paste0(x, "_", lag), all_of(variable)),
        join_by(GID_1, laggedYear == year)) |>
      select(-laggedYear)
    return(out)
  }
