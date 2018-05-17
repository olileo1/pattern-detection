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

robustscaleQn <- function(x, idx = NULL, prob = 0.5, mid = NULL) {
  if (!is.null(idx)) {
    sd <- s_Qn(x[idx])
    if (!is.null(mid)) {
      mean <- mid
    } else {
      mean <- quantile(x[idx], probs = prob)
    }
  } else {
    sd <- s_Qn(x)
    if (!is.null(mid)) {
      mean <- mid
    } else {
      mean <- quantile(x, probs = prob)
    }
    mean <- quantile(x, probs = prob)
  }
  return((x - mean) / sd)
}
