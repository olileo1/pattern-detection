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
    windows <- getSlidingWindows(n = n, window.length = window.lengths[i], step.length = step.lengths[i])
    windowFeatures(y = y, windows = windows, feature = feature, return.data.table = return.data.table)
  })
  if (return.data.table) {
    return(rbindlist(out))
  } else {
    return(out)
  }
}

patternSearch.gangnamstyle <- function(y,
                                       window.lengths = floor(c(0.1, 0.5) * length(y)),
                                       step.lengths = rep(floor(0.1 * length(y)), length(window.lengths)),
                                       pattern = list(x = c(0, 0.1, 1), y = c(0, 1, 0)),
                                       eligible.window.constraints = c('m1.max.sign.pat.coef > 0.1', 'm2.max.sign.pat.coef > 0.1',
                                                                       'm1.lrvar.scaled < 1', 'm1.lrvar.scaled < 1', 'm2.slope.coef.scaled < 2')) {
  featuretable <- slidingFeatures(y = y,
                                  window.lengths = window.lengths,
                                  step.lengths = step.lengths,
                                  feature = function(x) {patternFitting.fiesta(y = x, pat = pattern)},
                                  return.data.table = TRUE)
  set(featuretable, j = 'm1.max.sign.pat.coef.scaled', value = robustScalingMAD(featuretable[['m1.max.sign.pat.coef']]))
  set(featuretable, j = 'm2.max.sign.pat.coef.scaled', value = robustScalingMAD(featuretable[['m2.max.sign.pat.coef']]))
  set(featuretable, j = 'm1.lrvar.scaled', value = robustScalingMAD(featuretable[['m1.lrvar']]))
  set(featuretable, j = 'm2.lrvar.scaled', value = robustScalingMAD(featuretable[['m2.lrvar']]))
  set(featuretable, j = 'm2.slope.coef.scaled', value = robustScalingMAD(abs(featuretable[['m2.slope.coef']] / sqrt(featuretable[['m2.slope.coef.var']]))))
  set(featuretable, j = 'eligible.window', value = 1:nrow(featuretable) %in% featuretable[eval(parse(text = paste(eligible.window.constraints, collapse = ' & '))), which = TRUE])
  return(featuretable)
}

patternSearch.robust <- function(y,
                                 window.lengths = floor(c(0.1, 0.5) * length(y)),
                                 step.lengths = rep(floor(0.1 * length(y)), length(window.lengths)),
                                 pattern = list(x = c(0, 0.1, 1), y = c(0, 1, 0)),
                                 m1.pattern.threshold = 0,
                                 m2.pattern.threshold = 0.05,
                                 plot = NULL) {
  featuretable <- slidingFeatures(y = y,
                                  window.lengths = window.lengths,
                                  step.lengths = step.lengths,
                                  feature = function(x) {patternFitting.partyparty(y = x, pat = pattern)},
                                  return.data.table = TRUE)
  
  m1.outlier.model <- lmrob(featuretable[['m1.max.sign.pat.coef']] ~ 1)
  set(featuretable, j = 'm1.outliers', value = m1.outlier.model$rweights)
  m2.outlier.model <- lmrob(featuretable[['m2.max.sign.pat.coef']] ~ 1)
  set(featuretable, j = 'm2.outliers', value = m2.outlier.model$rweights)
  rmcd <- covMcd(cbind(featuretable[['m1.pat.coef']],featuretable[['m2.pat.coef']]))
  mdist <- mahalanobis(cbind(featuretable[['m1.pat.coef']],featuretable[['m2.pat.coef']]),
                       center = rmcd$center, cov = rmcd$cov  )
  set(featuretable, j = 'mdist', value = mdist)
  out <- featuretable[m1.max.sign.pat.coef > m1.pattern.threshold &
                        m2.max.sign.pat.coef > m2.pattern.threshold][order(m2.outliers)]
  if (nrow(out) == 0) {
    return(NULL)
  } else {
    if (!is.null(plot) && plot > 0) {
      p <- lapply(1:min(nrow(out), plot), function(i) {
        ggplot(data = data.frame(x = 1:length(y), y = y), aes(x = x, y = y)) +
          geom_rect(aes(xmin = out[i][['win.start']], xmax = out[i][['win.end']], ymin = -Inf, ymax = Inf)) +
          geom_line() +
          ggtitle(paste0('m2.pat.coef: ', round(out[i][['m2.pat.coef']], 3),
                         ' m2.outliers: ', round(out[i][['m2.outliers']], 3)))
      })
      do.call('grid.arrange', c(p, ncol = 2))
    }
    return(out)
  }
}