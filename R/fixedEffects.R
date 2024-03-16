#' @export
removeFixedEffects <- function(data, fixedEffects, exclude=NULL) {

  fixedEffectVariable <-
    fixedEffects |>
    str_extract_all("[\\w_]+") |>
    unlist() |>
    unique()
  numericVariables <-
    data |>
    select(where(is.numeric)) |>
    names()
  variables <- setdiff(numericVariables, c(fixedEffectVariable, exclude))

  Z <- data |>
    select(all_of(variables)) |>
    as.matrix()
  frmlRemoveFe <- formula(paste("Z ~", paste(fixedEffects, collapse = " + ")))
  fitLmFe <- lm(frmlRemoveFe, data = data)
  Z_noFe <- Z - fitLmFe$fitted.values
  dataNoFe <-
    bind_cols(
      data |> select(-all_of(variables)),
      as_tibble(Z_noFe))

  return(dataNoFe)
}
