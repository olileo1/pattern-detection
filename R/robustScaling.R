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
