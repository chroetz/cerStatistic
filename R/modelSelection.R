fitLmWithSummary <- function(data, predictors, target, seFun, coefFun, timeTrendOrder) {
  if (length(predictors) == 0) {
    res <- cbind(coef=NA, se=NA, t=NA, pVal=NA)
    fit <- lm(paste(target, "~ 0"), data = data)
  } else {
    frml <- formula(paste(
      target, "~",
      paste(predictors, collapse = " + "),
      "- 1"))
    fit <- lm(frml, data = data)
    se <- seFun(fit, data, predictors, timeTrendOrder)
    coef <- coefFun(fit, predictors)
    t <- coef / se
    pVal <- 2 * pt(abs(t), df = fit$df.residual-timeTrendOrder, lower.tail = FALSE)
    res <- cbind(coef=coef, se, t, pVal)
  }
  aic <- AIC(fit)
  bic <- BIC(fit)
  return(lst(res, aic, bic, fit))
}

#' @export
buildSeFun <- function(clusterFun, cumulate = NULL) {
  force(clusterFun)
  force(cumulate)
  seFun <- function(fit, data, predictors, timeTrendOrder) {
    vcov <- sandwich::vcovCL(fit, cluster = clusterFun(data))[predictors, predictors, drop=FALSE]
    correctionFactor <- getFeCorrectionFactor(data, predictors, timeTrendOrder)
    se <- sqrt(diag(vcov)*correctionFactor)
    if (length(cumulate) > 0) { # TODO: validate implementation
      lagNum <- str_extract(predictors, "lag\\d+") |> str_sub(4) |> as.integer()
      baseName <- str_remove(predictors, "_lag\\d+")
      lagNum[is.na(lagNum) | !(baseName %in% cumulate)] <- 0
      for (j in which(lagNum > 0)) {
        lagIdxs <- j:(j-lagNum[j]) # predictors must be ordered by lag
        variance <- sum(vcov[lagIdxs, lagIdxs])
        se[j] <- sqrt(variance*correctionFactor)
      }
    }
    return(se)
  }
  return(seFun)
}

#' @export
buildCoefFun <- function(cumulate = NULL) {
  force(cumulate)
  coefFun <- function(fit, predictors) {
    coef <- coef(fit)[predictors]
    if (length(cumulate) > 0) { # TODO: validate implementation
      lagNum <- str_extract(predictors, "lag\\d+") |> str_sub(4) |> as.integer()
      baseName <- str_remove(predictors, "_lag\\d+")
      lagNum[is.na(lagNum) | !(baseName %in% cumulate)] <- 0
      for (j in which(lagNum > 0)) {
        coef[j] <- coef[j] + coef[j-1] # predictors must be ordered by lag
        names(coef)[j] <- paste0("\U2211", names(coef)[j])
      }
    }
    return(coef)
  }
  return(coefFun)
}

#' @export
buildClusterFun <- function(clusterVariableNames) {
  force(clusterVariableNames)
  clusterFun <- function(data) {
    cluster <- rep("_", nrow(data))
    for (nm in clusterVariableNames) {
      cluster <- paste(cluster, data[[nm]], sep = "_")
    }
    return(cluster)
  }
  return(clusterFun)
}

getFeCorrectionFactor <- function(data, predictors, timeTrendOrder) {
  n <- nrow(data)
  p1 <- length(predictors)
  nFe <- getNumFe(data, timeTrendOrder)
  p2 <- nFe + length(predictors)
  correctionFactor <- (n-p1)/(n-p2)
  return(correctionFactor)
}

getNumFe <- function(data, timeTrendOrder) {
  nYears <- length(unique(data$theYear))
  nIDs <- length(unique(data$ID))
  nFe <- nYears + (timeTrendOrder+1)*(nIDs-1)
  return(nFe)
}

#' @export
fitAndPredictLm. <- function(xTrain, yTrain, xValidation) {
  fit <- .lm.fit(xTrain, yTrain)
  as.vector(xValidation %*% fit$coefficients)
}

#' @export
fitAndPredictLm <- function(xTrain, yTrain, xValidation) {
  coeff <- solve.default(crossprod(xTrain), crossprod(xTrain, yTrain))
  as.vector(xValidation %*% coeff)
}

cvPredPowOne <- function(x, y, cluster, nFolds, fitAndPredictFun, lossFun) {
  if (length(x) == 0) {
    return(rep(0, nFolds))
  }
  allClusters <- unique(cluster)
  clustersFoldId <- sample(rep(1:nFolds, length.out = length(allClusters)))
  sapply(seq_len(nFolds), \(i) {
    validationClusters <- allClusters[clustersFoldId == i]
    isValidation <- cluster %in% validationClusters
    yValidation <- y[isValidation]
    yTrain <- y[!isValidation]
    xValidation <- x[isValidation, , drop=FALSE]
    xTrain <- x[!isValidation, , drop=FALSE]
    esti <- fitAndPredictFun(xTrain, yTrain, xValidation)
    lossFun(yValidation, esti)
  })
}

kFoldCvOne <- function(x, y, cluster, nFolds, fitAndPredictFun, lambda = 1) {
  allClusters <- unique(cluster)
  clustersFoldId <- sample(rep(1:nFolds, length.out = length(allClusters)))
  cvResult <- lapply(seq_len(nFolds), \(i) {
    validationClusters <- allClusters[clustersFoldId == i]
    isValidation <- cluster %in% validationClusters
    yValidation <- y[isValidation]
    yTrain <- y[!isValidation]
    if (length(x) == 0) {
      loss <- vapply(lambda, \(lmbd) yValidation^2, FUN.VALUE=numeric(length(yValidation)))
      referenceLoss <- loss
    } else {
      xValidation <- x[isValidation, , drop=FALSE]
      xTrain <- x[!isValidation, , drop=FALSE]
      esti <- fitAndPredictFun(xTrain, yTrain, xValidation)
      loss <- vapply(lambda, \(lmbd) (yValidation - lmbd*esti)^2, FUN.VALUE=numeric(length(yValidation)))
      referenceLoss <- yValidation^2
    }
    tibble(
      loss = loss,
      referenceLoss = referenceLoss,
      cluster = cluster[isValidation],
      foldNr = i)
  }) |>
    bind_rows()
  return(cvResult)
}

scaleToUnit <- function(x) {
  (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE))
}

colVars <- function(x) {
  colMeans((x - rep(colMeans(x), each = nrow(x)))^2)
}

kFoldCv <- function(
    data,
    predictors,
    target,
    nReps,
    cluster,
    nFolds,
    fitAndPredictFun,
    lambda
) {
  x <- data |> select(all_of(predictors)) |> as.matrix()
  y <- data |> pull(target)

  cvLossTable <- lapply(seq_len(nReps), \(repNr) {
    lossTable <- kFoldCvOne(
      x, y, cluster,
      nFolds, fitAndPredictFun, lambda)
    clusteredLossTable <-
      lossTable |>
      summarise(
        loss = matrix(colMeans(loss - referenceLoss), nrow=1),
        referenceLoss = mean(referenceLoss),
        foldNr = foldNr[1],
        .by = cluster)
    varIn <-
      clusteredLossTable |>
      summarise(
        var = matrix(colVars(loss), nrow = 1),
        .by = foldNr
      ) |>
      pull(var) |>
      colMeans()
    tibble(
      repNr = repNr,
      lambda = lambda,
      mean = colMeans(clusteredLossTable$loss),
      reference = mean(clusteredLossTable$referenceLoss),
      varOut = colVars(clusteredLossTable$loss),
      varIn = varIn)
  }) |> bind_rows()

  return(cvLossTable)
}


getModelListSelectionValues <- function(
    data,
    predictorList,
    target,
    nReps,
    nFolds,
    seFun,
    coefFun,
    timeTrendOrder,
    fitAndPredictFun,
    cluster,
    lambda
) {

  nCluster <- length(unique(cluster))
  cvResults <- lapply(seq_along(predictorList), \(i) {
    preds <- predictorList[[i]]
    cvRes <- kFoldCv(
      data,
      preds,
      target,
      nReps,
      cluster,
      nFolds,
      fitAndPredictFun,
      lambda
    ) |>
    summarize(
      CV = mean(mean),
      CV_sigma = sqrt(mean(varIn) / nCluster),
      .by = lambda
    ) |>
    mutate(
      CV0 = CV,
      CV1 = CV + CV_sigma,
      CV2 = CV + 2*CV_sigma,
      CV3 = CV + 3*CV_sigma)
    tibble(
      predictors = list(preds),
      CV = cvRes |> filter(lambda == 1) |> pull(CV) |> first(),
      CVsd = cvRes |> filter(lambda == 1) |> pull(CV_sigma) |> first(),
      CV0 = cvRes |> filter(CV0 == min(CV0)) |> pull(CV) |> first(),
      CV0sd = cvRes |> filter(CV0 == min(CV0)) |> pull(CV_sigma) |> first(),
      CV0lambda = cvRes |> filter(CV0 == min(CV0)) |> pull(lambda) |> first(),
      CV1 = cvRes |> filter(CV1 == min(CV1)) |> pull(CV) |> first(),
      CV1sd = cvRes |> filter(CV1 == min(CV1)) |> pull(CV_sigma) |> first(),
      CV1lambda = cvRes |> filter(CV1 == min(CV1)) |> pull(lambda) |> first(),
      CV2 = cvRes |> filter(CV2 == min(CV2)) |> pull(CV) |> first(),
      CV2sd = cvRes |> filter(CV2 == min(CV2)) |> pull(CV_sigma) |> first(),
      CV2lambda = cvRes |> filter(CV2 == min(CV2)) |> pull(lambda) |> first(),
      CV3 = cvRes |> filter(CV3 == min(CV3)) |> pull(CV) |> first(),
      CV3sd = cvRes |> filter(CV3 == min(CV3)) |> pull(CV_sigma) |> first(),
      CV3lambda = cvRes |> filter(CV3 == min(CV3)) |> pull(lambda) |> first())
  }) |>
    bind_rows()
  scale <- if (min(cvResults$CV0) >= 0) diff(range(cvResults$CV)) else -min(cvResults$CV0)
  cvResults <- cvResults |>
    mutate(across(matches("^CV\\d*(sd)?$"), \(x) x/scale))

  predictorsMax <- sapply(predictorList, length) |> max()

  rowsPval <- lapply(seq_along(predictorList), \(i) {
    preds <- predictorList[[i]]
    res <- fitLmWithSummary(data, preds, target, seFun, coefFun, timeTrendOrder)
    tbl <-
      rep(NA_character_, predictorsMax) |>
      setNames(paste0("name_", seq_len(predictorsMax))) |>
      as.list() |>
      as_tibble()
    for (j in seq_len(nrow(res$res))) {
      nm <- rownames(res$res)[j]
      tbl[1,paste0("name_", j)] <- if (is.null(nm)) NA else nm
      tbl[1,paste0("pValue_", j)] <- res$res[j, "pVal"]
      tbl[1,paste0("coef_", j)] <- res$res[j, "coef"]
    }
    ic <- list(aic = res$aic, bic = res$bic)
    bind_cols(tbl, as_tibble(ic))
  })

  pValTable <- bind_rows(rowsPval)

  res <- bind_cols(
    pValTable,
    cvResults
  ) |>
    mutate( # first model is reference and has value 0 for information criteria
      aic = aic - aic[1],
      bic = bic - bic[1])

  return(res)
}

#' @export
getModelListSelectionTable <- function(
    data,
    predictorList,
    target,
    nReps,
    nFolds,
    seFun,
    coefFun,
    timeTrendOrder,
    fitAndPredictFun,
    cluster,
    lambda,
    csvPath
) {

  res <-
    getModelListSelectionValues(
      data,
      predictorList,
      target,
      nReps,
      nFolds,
      seFun,
      coefFun,
      timeTrendOrder,
      fitAndPredictFun,
      cluster,
      lambda
    ) |>
    mutate(AIC = scaleToUnit(aic), BIC = scaleToUnit(bic))

  if (!dir.exists(dirname(csvPath))) dir.create(dirname(csvPath))
  write_csv(res, csvPath)

  resShow <- res
  predictorsMax <- sapply(predictorList, length) |> max()
  for (j in seq_len(predictorsMax)) {
    resShow[[paste0("X", j)]] <- sprintf(
      "%s<br>%.3g<br>(p: %.4f)",
      res[[paste0("name_", j)]],
      res[[paste0("coef_", j)]],
      res[[paste0("pValue_", j)]])
  }
  resShow$CV <- sprintf("%.3g<br>\U00B1%.3g", res$CV, res$CVsd)
  for (j in 0:3) {
    resShow[[paste0("CV", j)]] <- sprintf(
      "%.3g<br>\U00B1%.3g<br>(s: %.2f)",
      res[[paste0("CV", j)]],
      res[[paste0("CV", j, "sd")]],
      res[[paste0("CV", j, "lambda")]])
  }
  resShow <-
    resShow |>
    select(starts_with("X"), AIC, BIC, CV, all_of(paste0("CV", 0:3)))
  resValue <-
    res |>
    select(starts_with("pValue_"), AIC, BIC, CV, all_of(paste0("CV", 0:3))) |>
    mutate(across(starts_with("CV"), scaleToUnit))

  kblTable <- resShow |>
    makeTable("Model Selection") |>
    colorizeTable(seq_len(predictorsMax), resValue) |>
    colorizeTable2(predictorsMax + 1:7, resValue)
  return(kblTable)
}
