#' @export
analyzeCorrelation <- function(data, predictors, target, variable) {

  x <- data |> select(all_of(predictors)) |> as.matrix()
  y <- data |> pull(target)
  fit <- .lm.fit(x, y)
  residuals <- (y - x %*% fit$coefficients) |> as.vector()
  data <- data |> mutate(residual = residuals)

  dataIdWide <-
    data |>
    select(ID, theYear, !!variable) |>
    pivot_wider(names_from=ID, values_from=!!variable) |>
    arrange(theYear)
  dataYearWide <-
    data |>
    select(ID, theYear, !!variable) |>
    arrange(theYear) |>
    pivot_wider(names_from=theYear, values_from=!!variable) |>
    arrange(ID)

  corrMatYear <- cor(
    dataYearWide |> select(-ID),
    use = "pairwise.complete.obs")

  corrMatId <- cor(
    dataIdWide |> select(-theYear),
    use = "pairwise.complete.obs")

  return(lst(corrMatId, corrMatYear))
}


#' @export
showCorrelations <- function(corrMatId, corrMatYear, hclustMethod, hclustCut) {

  cat("\n\n**Correlation Between Regions**\n\n")

  corrMatId[!is.finite(corrMatId)] <- 0
  hcluster <- hclust(as.dist(1-corrMatId), method = hclustMethod)
  corrMatIdOrdered <- corrMatId[hcluster$order, hcluster$order]
  showCorrImage(corrMatIdOrdered)

  clusters <- getClusters(hcluster, hclustCut)
  meanCorrelations <- sapply(clusters, \(cluster) {
    subMat <- corrMatId[cluster, cluster]
    mean(subMat[upper.tri(subMat, diag=FALSE)])
  })
  correlationClusterTable <- tibble(
    IDs =  sapply(clusters, paste, collapse = ","),
    n = sapply(clusters, length),
    meanCorrelation = meanCorrelations
  ) |>
    arrange(desc(meanCorrelation))
  correlationClusterTable |> makeTable("Correlation Cluster Table") |> cat("\n\n")

  cat("\n\n**Correlation Between Years**\n\n")

  showCorrImage(corrMatYear)

  n <- nrow(corrMatYear)
  autoCorData <- lapply(1:(n-1), autoCor, corrMatYear) |> bind_rows()

  cat("\n\n")
  plt <- autoCorData |>
    ggplot(aes(x=k, y=cor)) +
    geom_point() + geom_smooth() +
    geom_point(data = autoCorData |> summarize(cor = mean(cor), .by = k), color="red", stroke = 2, size=2, shape =4) +
    geom_hline(yintercept = 0, color = "green") +
    labs(title = "Autocorrelation of Residuals", x = "Lag", y = "Correlation")
  plot(plt)
  cat("\n\n")
}


getClusters <- function(hcluster, h) {
  tr <- cutree(hcluster, h = h)
  tbl <- table(tr) |> as.vector()
  clusterIds <- tr[tr %in% which(tbl > 1)]
  clusterIdUnique <- unique(clusterIds)
  clusters <- lapply(clusterIdUnique, \(id) names(clusterIds)[clusterIds == id])
  clusters
}


showCorrImage <- function(mat) {
  cat("\n\n")
  par(mar = c(0,0,0,0))
  image(
    z = mat,
    zlim = c(-1, 1),
    col = colorRampPalette(c("darkblue", "white", "darkred"))(255),
    useRaster = TRUE,
    asp = 1,
    main = NULL,
    axes = FALSE,
    frame.plot = FALSE)
  cat("\n\n")
}


autoCor <- function(k, mat) {
  n <- nrow(mat)
  v <- mat[cbind(1:(n-k), (1+k):n)]
  tibble(k = k, cor = v)
}

