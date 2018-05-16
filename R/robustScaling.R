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

robustScalingMAD <- function(x) {
  med <- median(x)
  mad <- median(abs(x - med))
  return((x - med) / mad)
}

mid <- function(x, alpha.begin = 0.1, alpha.end = 0.1) {
  n <- length(x)
  start <- ceiling(alpha.begin * n)
  end <- floor((1 - alpha.end) * n)
  return(x[max(1, start) : min(n, end)])
}
