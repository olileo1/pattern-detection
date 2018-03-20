windowFeatures <- function(y, windows = matrix(c(1, length(y)), ncol =2), feature = function(x) c(mean = mean(x), sd = sd(x)), return.data.table = TRUE) {
  out <- lapply(1:dim(windows)[1], function(i) {
    if (is.null(dim(y))) {
      n <- length(y)
      append(as.list(feature(y[max(1, windows[i, 1]):min(windows[i, 2], n)])),
             list(
               win.start = windows[i, 1],
               win.end = windows[i, 2]
             ))
    } else {
      n <- dim(y)[1]
      append(as.list(feature(y[max(1, windows[i, 1]):min(windows[i, 2], n),])),
             list(
               win.start = windows[i, 1],
               win.end = windows[i, 2]
             )) 
    }
  })
  if (return.data.table) {
    return(rbindlist(lapply(out, as.data.table)))
  } else {
    return(out)
  }
}

getSlidingWindows <- function(n, window.length, step.length) {
  if (window.length > n) return(matrix(c(1, n), ncol = 2))
  n.windows <- ceiling((n - window.length) / step.length)
  out <- cbind(1 + c(0:n.windows) * step.length, window.length + c(0:n.windows) * step.length)
  out[dim(out)[1], dim(out)[2]] <- n
  return(out)
}

slidingFeatures <- function(y, window.lengths, step.lengths,
                            feature = function(x) mean(x),
                            return.data.table = TRUE) {
  if (is.null(dim(y))) {
    n <- length(y)
  } else {
    n <- dim(y)[1]
  }
  out <- lapply(1:length(window.lengths), function(i) {
    windows <- get.timewins(n = n, window.length = window.lengths[i], step.length = step.lengths[i])
    windowFeatures(y = y, windows = windows, feature = feature, return.data.table = return.data.table)
  })
  if (return.data.table) {
    return(rbindlist(out))
  } else {
    return(out)
  }
}
