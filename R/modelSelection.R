fitLmWithSummary <- function(data, predictors, target, seFun, timeTrendOrder) {
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
    t <- coef(fit) / se
    pVal <- 2 * pt(abs(t), df = fit$df.residual-timeTrendOrder, lower.tail = FALSE)
    res <- cbind(coef=coef(fit), se, t, pVal)
  }
  aic <- AIC(fit)
  bic <- BIC(fit)
  return(lst(res, aic, bic, fit))
}

#' @export
buildSeFun <- function(clusterFun) {
  force(clusterFun)
  seFun <- function(fit, data, predictors, timeTrendOrder) {
    vcov <- sandwich::vcovCL(fit, cluster = clusterFun(data))[predictors, predictors, drop=FALSE]
    correctionFactor <- getFeCorrectionFactor(data, predictors, timeTrendOrder)
    se <- sqrt(diag(vcov)*correctionFactor)
    return(se)
  }
  return(seFun)
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

scaleToUnit <- function(x) {
  (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE))
}

kFoldCv <- function(
    data,
    predictors,
    target,
    nReps,
    nFolds,
    fitAndPredictFun,
    lossFun,
    clusterVariableName
) {
  x <- data |> select(all_of(predictors)) |> as.matrix()
  y <- data |> pull(target)
  cv <- replicate(
    nReps,
    cvPredPowOne(
      x, y,
      cluster=data[[clusterVariableName]],
      nFolds=nFolds,
      fitAndPredictFun=fitAndPredictFun,
      lossFun=lossFun
    ),
    simplify = FALSE)
  cvFoldMean <- sapply(cv, mean)
  res <- list()
  res$cvMean <- mean(cvFoldMean)
  res$cvMeanSd <- sd(cvFoldMean) / sqrt(length(cvFoldMean))
  return(as_tibble(res))
}

modelSelectionRowsToResults <- function(rows) {
  rows |>
    bind_rows() |>
    mutate(
      CV = cvMean,#scaleToUnit(cvMean),
      CVsd = cvMeanSd)#/(max(cvMean, na.rm=TRUE)-min(cvMean, na.rm=TRUE)))
}


getModelListSelectionValues <- function(
    data,
    predictorList,
    target,
    nReps,
    nFolds,
    seFun,
    timeTrendOrder,
    fitAndPredictFun,
    lossFun,
    clusterVariableName
) {

  rows <- lapply(seq_along(predictorList), \(i) {
    preds <- predictorList[[i]]
    cvRes <- kFoldCv(
      data,
      preds,
      target,
      nReps,
      nFolds,
      fitAndPredictFun,
      lossFun,
      clusterVariableName
    )
    bind_cols(tibble(predictors = list(preds)), cvRes)
  })

  modelSelectionResults <- modelSelectionRowsToResults(rows)

  predictorsMax <- sapply(predictorList, length) |> max()

  rowsPval <- lapply(seq_along(predictorList), \(i) {
    preds <- predictorList[[i]]
    res <- fitLmWithSummary(data, preds, target, seFun, timeTrendOrder)
    tbl <-
      rep(NA_character_, predictorsMax) |>
      setNames(paste0("name_", seq_len(predictorsMax))) |>
      as.list() |>
      as_tibble()
    for (j in seq_len(nrow(res$res))) {
      nm <- rownames(res$res)[j]
      tbl[1,paste0("name_", j)] <- if (is.null(nm)) NA else nm
      tbl[1,paste0("pValue_", j)] <- res$res[j, "pVal"]
      tbl[1,paste0("tStatistic_", j)] <- res$res[j, "t"]
      coef <- res$res[j, "coef"]
      tbl[1,paste0("coef_", j)] <- coef
      if (is.na(coef)) {
        tbl[1,paste0("sign_", j)] <- "?"
      } else {
        tbl[1,paste0("sign_", j)] <- if (coef > 0) "+" else "-"
      }
    }
    ic <- list(aic = res$aic, bic = res$bic)
    bind_cols(tbl, as_tibble(ic))
  })

  pValTable <- bind_rows(rowsPval)

  res <- bind_cols(
    pValTable,
    modelSelectionResults |> select(CV, CVsd)) |>
    mutate( # first model is reference and has value 0
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
    timeTrendOrder,
    fitAndPredictFun,
    lossFun,
    clusterVariableName,
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
      timeTrendOrder,
      fitAndPredictFun,
      lossFun,
      clusterVariableName
    ) |>
    mutate(AIC = scaleToUnit(aic), BIC = scaleToUnit(bic))

  predictorsMax <- sapply(predictorList, length) |> max()
  for (j in seq_len(predictorsMax)) {
    res[[paste0("X", j)]] <- sprintf(
      "%s<br>%.1g (%.4f)",
      res[[paste0("name_", j)]],
      res[[paste0("coef_", j)]],
      res[[paste0("pValue_", j)]])
  }
  resShow <-
    res |>
    select(starts_with("X"), AIC, BIC, CV, CVsd)
  if (!dir.exists(dirname(csvPath))) dir.create(dirname(csvPath))
  write_csv(resShow, csvPath)
  resValue <-
    res |>
    select(starts_with("pValue_"))
  kblTable <- resShow |>
    makeTable("Model Selection") |>
    colorizeTable(seq_len(predictorsMax), resValue) |>
    colorizeTable2(predictorsMax + 1:3, resShow)
  return(kblTable)
}


getModelListSelectionHalfVariationValues <- function(
    data,
    predictorList,
    target,
    nRepsInner,
    nFolds,
    seFun,
    timeTrendOrder,
    fitAndPredictFun,
    lossFun,
    clusterVariableName,
    nRepsOuter
) {

  cluster <- data[[clusterVariableName]]
  allClusters <- unique(cluster)

  reps <- replicate(
    nRepsOuter,
    {
      halfClusters <- sample(allClusters, length(allClusters) %/% 2)
      isHalf <- cluster %in% halfClusters
      data1 <- data[isHalf, ]
      data2 <- data[!isHalf, ]
      res1 <- getModelListSelectionValues(
        data1,
        predictorList,
        target,
        nRepsInner,
        nFolds,
        seFun,
        timeTrendOrder,
        fitAndPredictFun,
        lossFun,
        clusterVariableName
      )
      res2 <- getModelListSelectionValues(
        data2,
        predictorList,
        target,
        nRepsInner,
        nFolds,
        seFun,
        timeTrendOrder,
        fitAndPredictFun,
        lossFun,
        clusterVariableName
      )
      return(list(res1, res2))
    },
    simplify = FALSE
  )

  return(reps)
}


#' @export
getModelListSelectionHalfVariationResultTables <- function(
    data,
    predictorList,
    target,
    nRepsInner,
    nFolds,
    seFun,
    timeTrendOrder,
    fitAndPredictFun,
    lossFun,
    clusterVariableName,
    nRepsOuter
) {
  pairsList <- getModelListSelectionHalfVariationValues(
    data,
    predictorList,
    target,
    nReps = nRepsInner,
    nFolds = nFolds,
    seFun = seFun,
    timeTrendOrder = timeTrendOrder,
    fitAndPredictFun = fitAndPredictFun,
    lossFun = lossFun,
    clusterVariableName = clusterVariableName,
    nRepsOuter = nRepsOuter
  )
  vars <- lapply(pairsList , \(pairs) {
    x1 <- pairs[[1]]
    x2 <- pairs[[2]]
    nms <- names(x1)
    xVar <- lapply(
      nms,
      \(nm) if (is.numeric(x1[[nm]])) (x1[[nm]] - x2[[nm]])^2/2 else x1[[nm]] == x2[[nm]])
    names(xVar) <- nms
    as_tibble(xVar)
  })
  means <- lapply(pairsList , \(pairs) {
    x1 <- pairs[[1]]
    x2 <- pairs[[2]]
    nms <- names(x1)
    xMean <- lapply(
      nms,
      \(nm) if (is.numeric(x1[[nm]])) (x1[[nm]] + x2[[nm]])/2 else x1[[nm]])
    names(xMean) <- nms
    as_tibble(xMean)
  })
  var1 <- vars[[1]]
  nms <- names(var1)

  meanOfVars <- lapply(nms, \(nm) sapply(vars, \(x) x[[nm]]) |> rowMeans())
  names(meanOfVars) <- nms
  meanOfVars <- as_tibble(meanOfVars)

  varOfVars <- lapply(nms, \(nm) sapply(vars, \(x) x[[nm]]) |> apply(1, var))
  names(varOfVars) <- nms
  varOfVars <- as_tibble(varOfVars)

  mean1 <- means[[1]]
  nms <- names(mean1)

  meanOfMeans <- lapply(nms, \(nm) {
    if (is.numeric(mean1[[nm]])) {
      sapply(means, \(x) x[[nm]]) |> rowMeans()
    } else {
      means[[1]][[nm]]
    }
  })
  names(meanOfMeans) <- nms
  meanOfMeans <- as_tibble(meanOfMeans)

  varOfMeans <- lapply(nms, \(nm) {
    if (is.numeric(mean1[[nm]])) {
      sapply(means, \(x) x[[nm]]) |> apply(1, var)
    } else {
      sapply(means, \(x) x[[nm]]) |> apply(1, \(x) length(unique(x[!is.na(x)])))
    }
  })
  names(varOfMeans) <- nms
  varOfMeans <- as_tibble(varOfMeans)

  return(lst(meanOfMeans, varOfMeans, meanOfVars, varOfVars))
}
