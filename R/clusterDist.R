mclustDist <- function(X) {
  X.scaled <- apply(X, 2, function(x) {
    robustscaleQn(x)
  })
  X.clean <- X.scaled[apply(X.scaled, 1, function(x) all(x >= -7 & x <= 7)),]
  #m <- Qclust(X.clean, K = 3, q = 0.9, modelNames = 'VVV')
  m <- Mclust(X.clean, G = c(1,2,3), modelNames = 'VVV', verbose = FALSE)
  groups <- table(m$classification)
  pivot.group <- which(groups == max(groups))[1]
  eligible.groups <- setdiff(which(groups > (dim(X.clean)[1] / 10)), pivot.group)
  dist <- 0
  if (length(eligible.groups) > 0) {
    for (g in eligible.groups) {
      dist <- c(dist,
                mahalanobis(x = t(as.matrix(m$parameters$mean[, g])),
                            center = m$parameters$mean[, pivot.group],
                            cov = m$parameters$variance$sigma[,,pivot.group]))
    }
  }
  return(max(dist))
}