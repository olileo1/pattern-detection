robustMoments <- function(x, alpha = 0.1) {
  cutoff <- length(x) * alpha
  tmp <- x[order(abs(x - median(x)))[1:floor(length(x)*(1-alpha))]]
  rob.mean <- mean(tmp)
  rob.sd <- mean((tmp - rob.mean) ^ 2)
  return(list(mean = rob.mean,
              sd = rob.sd))
}

robustScaling <- function(x, alpha = 0.1) {
  tmp <- robustMoments(x = x, alpha = alpha)
  return((x - tmp$mean) / (tmp$sd))
}

robustscaleMAD <- function(x) {
  med <- median(x)
  mad <- median(abs(x - med))
  return((x - med) / mad)
}

modals <- function(x, nmods = 2, orderout = 2) {
  m <- tclust(x = x, k = nmods)
  return(sort(m$centers)[orderout])
}

robustscaleQn <- function(x, idx = seq_along(x), mid = NULL, mid.value = 0, prob = 0.5, nmods = 2, orderout = 1) {
  mean <- switch(mid,
                 quantile = quantile(x[idx], probs = prob),
                 modal = modals(x[idx], nmods = nmods, orderout = orderout),
                 fixed = mid.value)
  sd <- s_Qn(x[idx])
  return((x - mean) / sd)
}
