mclust.selection <- function(pattern.coef, pattern.error, pattern.error.bounds = c(-Inf, 1.5), pattern.coef.bounds = c(3, Inf)) {
  pattern.error.scaled <- robustscaleQn(x = pattern.error, mid = 'modal', nmods = 2, orderout = 1)
  pattern.coef.scaled <- robustScaleQn(x = pattern.coef, idx = pattern.error.scaled < pattern.error.bounds[2], mid = 'fixed', mid.value = 0)
  X <- cbind(pattern.coef.scaled,
             pattern.error.scaled)
  m <- Mclust(X, G = c(1, 2, 3), modelNames = 'VVV')
  if (m$G > 1) {
    center.norms <- apply(m$parameters$mean, 2, function(x) sum(x ^ 2))
    base_idx <- which(center.norms == min(center.norms))
    good_idx <- which(apply(m$parameters$mean[, -base_idx], 2, function(x) {
      x[1] > pattern.coef.bounds[1] & x[1] < pattern.coef.bounds[2]&
        x[2] > pattern.error.bounds[1] & x[2] < pattern.error.bounds[2]
    }))
    out <- ifelse(m$classification == base_idx, 'base',
                  ifelse(m$classification %in% good_idx, 'good', 'residuals'))
  } else {
    out <- rep('base', length(pattern.coef))
  }
  return(data.frame(out = out,
                    pattern.coef.scaled = pattern.coef.scaled,
                    pattern.error.scaled = pattern.error.scaled))
}

robust.selection <- function(pattern.coef, pattern.error, pattern.error.bounds = c(-Inf, 1.5), pattern.coef.bounds = c(3, Inf)) {
  pattern.error.scaled <- robustscaleQn(x = pattern.error, mid = 'modal', nmods = 2, orderout = 1)
  pattern.coef.scaled <- robustScaleQn(x = pattern.coef, idx = pattern.error.scaled < pattern.error.bounds[2], mid = 'fixed', mid.value = 0)
  dist <- qclustDist(X = cbind(pattern.coef.scaled, pattern.error.scaled))
  if (dist > 10) {
    out <- ifelse(pattern.error.scaled > pattern.error.bounds[2] | pattern.error.scaled < pattern.error.bounds[1], 'residuals',
                  ifelse(pattern.coef.scaled > pattern.coef.bounds[1] & pattern.coef.scaled < pattern.coef.bounds[2], 'good', 'base'))
  } else {
    out <- 'base'
  }
  return(data.frame(out = out,
                    pattern.coef.scaled = pattern.coef.scaled,
                    pattern.error.scaled = pattern.error.scaled))
}