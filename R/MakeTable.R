#' @export
makeTable <- function(data, caption, digits = 4, leftColNums = NULL) {
  align <- ifelse(
    seq_along(data) %in% leftColNums,
    "l",
    "r"
  )
  data |>
    kableExtra::kbl(
      caption = caption,
      digits = digits,
      booktabs = TRUE,
      align = align,
      escape = FALSE,
      table.attr = paste0(
        'class="table table-condensed table-striped" ',
        'style = "',
        'white-space: nowrap;',
        'width: auto !important;',
        'margin-left:',
        'auto; margin-right: auto;"'))
}

valueToColor <- function(value) {
  colors <- case_when(
    is.na(value) ~ rgb(0.6,0.6,0.6),
    value > 0.1 ~ rgb(1,0.5,0.5),
    value <= 0.1 & value > 0.05 ~ rgb(1,0.75,0.5),
    value <= 0.05 & value > 0.01 ~ rgb(1.0,1,0.5),
    value <= 0.01 ~ rgb(0.5,1,0.5))
  colors[is.na(colors)] <- rgb(0.6,0.6,0.6)
  colors
}

colorizeTable <- function(k, colorCols, valueTbl) {
  for (colId in colorCols) {
    k <- kableExtra::column_spec(
      k,
      colId,
      background = valueToColor(valueTbl[[colId]]))
  }
  k
}

valueToColor2 <- function(value) {
  colorsFun <- colorRamp(c("#B0FFFF", "#00FFFF", "#FFA0FF", "#BF78BF", "#7F507F"))
  colStr <- apply(
    colorsFun(value^(0.7)),
    1,
    function(x) {
      if (any(is.na(x))) rgb(0.6,0.6,0.6) else rgb(x[1]/256, x[2]/256, x[3]/256)
    })
  colStr[value == 0] <- "#FFFFFF"
  return(colStr)
}

colorizeTable2 <- function(k, colorCols, valueTbl) {
  for (colId in colorCols) {
    k <- kableExtra::column_spec(
      k,
      colId,
      background = valueToColor2(valueTbl[[colId]]))
  }
  k
}
