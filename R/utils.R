mid <- function(x, alpha.begin = 0.1, alpha.end = 0.1) {
  n <- length(x)
  start <- ceiling(alpha.begin * n)
  end <- floor((1 - alpha.end) * n)
  return(x[max(1, start) : min(n, end)])
}