#' @export
renderModelSelection <- function(
  outFileName,
  outDirPath = getwd(),
  outFormat = "HTML",
  envir = new.env(),
  quiet = FALSE,
  ...
) {
  yamlParams <- list(...)
  paths <- str_subset(names(yamlParams), "Path$")
  for (path in paths) {
    yamlParams[[path]] <- normalizePath(yamlParams[[path]])
  }

  rmdSourceFilePath <- system.file("rmarkdown/ModelSelection.Rmd", package = "cerStatistic")

  if (!dir.exists(outDirPath)) dir.create(outDirPath, recursive = TRUE)

  outFormat <- tolower(outFormat)[[1]]
  if (outFormat == "pdf") {
    outFormat <- "pdf_document"
  } else if (outFormat == "html") {
    outFormat <- "html_document"
  } else if (outFormat == "rmd") {
    return(.modelSelectionRmd(yamlParams, rmdSourceFilePath, outDirPath, outFileName))
  } else {
    stop("Unknown format: ", outFormat)
  }

  yamlParams <- expressionsToObject(yamlParams)

  rmarkdown::render(
    rmdSourceFilePath,
    intermediates_dir = tempdir(),
    output_dir = outDirPath,
    output_file = outFileName,
    output_format = outFormat,
    params = yamlParams,
    envir = envir,
    quiet = quiet)
}


.modelSelectionRmd <- function(yamlParams, rmdSourceFilePath, outDirPath, outFileName) {
  linesMain <- readLines(rmdSourceFilePath)
  delimiters <- grep("^(---|\\.\\.\\.)\\s*$", linesMain)
  headerMain <- linesMain[(delimiters[1]):(delimiters[2])]
  yml <- yaml::yaml.load(
    headerMain,
    handlers = list(r = function(x) ymlthis::yml_params_code(!!rlang::parse_expr(x))))
  yamlParams <- nonBaseElementsToExpressions(yamlParams)
  yamlParams <- expressionsToYmlCode(yamlParams)
  baseYaml <- ymlthis::as_yml(yml)
  newYamlParams <- baseYaml$params
  newYamlParams[names(yamlParams)] <- yamlParams
  newYaml <- ymlthis::yml_replace(
    baseYaml,
    params = newYamlParams,
    date = format(Sys.Date()))

  rmdDstFilePath <- file.path(outDirPath, paste0(outFileName, ".Rmd"))
  file.copy(rmdSourceFilePath, rmdDstFilePath)
  ymlthis::use_rmarkdown(
    newYaml,
    path = rmdDstFilePath,
    template = rmdSourceFilePath,
    include_yaml = FALSE,
    overwrite = TRUE,
    quiet = TRUE)
}


expressionsToYmlCode <- function(params) {
  if (length(params) == 0) return(params)
  selExpr <- sapply(params, \(x) typeof(x) %in% c("symbol", "language"))
  params[selExpr] <- lapply(
    params[selExpr],
    \(x) ymlthis::yml_params_code(!!x)
  )
  return(params)
}

expressionsToObject <- function(params) {
  if (length(params) == 0) return(params)
  selExpr <- sapply(params, \(x) typeof(x) %in% c("symbol", "language"))
  params[selExpr] <- lapply(
    params[selExpr],
    rlang::eval_bare
  )
  return(params)
}

nonBaseElementsToExpressions <- function(params) {
  if (length(params) == 0) return(params)
  selExpr <- sapply(params, \(x) typeof(x) %in% c("symbol", "language"))
  selBase <- sapply(params, \(x) (is.null(names(x)) && is.atomic(x) && length(x) == 1))
  sel <- !selExpr & !selBase
  params[sel] <- lapply(params[sel], objectToExpression)
  return(params)
}

objectToExpression <- function(x) {
  dput(x) |>
    capture.output() |>
    paste(collapse=" ") |>
    rlang::parse_expr()
}
