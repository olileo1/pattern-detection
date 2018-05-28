#' For a given timeseries y as vector and a set of window start and end points given as matrix this function applies
#' the provided feature function on each window of the timeseries.
#'
#' @param y univariate timeseries as vector
#' @param windows start and end points of the windows given as matrix (each row is a window)
#' @param feature function that takes a vector as Input and returns certain features as list
#' 
#' @return tibble where each row contains the features calculated on a certain window
windowFeatures <- function(y, windows = matrix(c(1, length(y)), ncol = 2), feature = function(x) c(mean = mean(x), sd = sd(x))) {
  out <- lapply(1:dim(windows)[1], function(i) {
    n <- length(y)
    append(as.list(feature(y[max(1, windows[i, 1]):min(windows[i, 2], n)])),
           list(
             win.start = windows[i, 1],
             win.end = windows[i, 2]
           ))
  })
  return(bind_rows(lapply(out, as_tibble)))
}

#' For a given length n, window.length and step.length this function create a matrix containing
#' the start and end points for the sliding windows
#'
#' @param n positive integer giving the length of the series
#' @param window.length integer giving the width of the sliding window
#' @param step.length integer giving the step.length with which the window gets slided through the series
#' 
#' @return matrix where each row contains start and end point of a window
getSlidingWindows <- function(n, window.length, step.length) {
  if (window.length > n) return(matrix(c(1, n), ncol = 2))
  n.windows <- ceiling((n - window.length) / step.length)
  out <- cbind(1 + c(0:n.windows) * step.length, window.length + c(0:n.windows) * step.length)
  if ((n - out[dim(out)[1], 1]) < (window.length / 2)) {
    out <- out[-(dim(out)[1]),]
  } else {
    out[dim(out)[1], 2] <- n
  }
  return(out)
}

#' For a given timeseries y, given window lengths window.lengths and corresponding step.lengths
#' this function applies the provided feature function to each set of sliding windows defined by 
#' window.lengths and step.lengths
#'
#' @param y univariat timeseries
#' @param window.lengths integer vector of length m containing the window lengths
#' @param step.length integer vector of length m containing the step length for each window length
#' @param feature function that takes a vector as Input and returns certain features as list
#' 
#' @return data.frame where each row contains the features calculated on a certain window
slidingFeatures <- function(y, window.lengths, step.lengths,
                            feature = function(x) mean(x)) {
  n <- length(y)
  out <- lapply(1:length(window.lengths), function(i) {
    if (n > window.lengths[i]) {
      windows <- getSlidingWindows(n = n, window.length = window.lengths[i], step.length = step.lengths[i])
      windowFeatures(y = y, windows = windows, feature = feature)
    }
  })
  return(bind_rows(out))
}

#' This function performs a patternSearch on a timeseries y. Window lengths and corresponding step lengths
#' have to be supported by the user, as well as the coordinates of the desired pattern. The normalization
#' method for the windows can be selected among a variety of alternatives as well the linear regression method.
#' A function measuring the extend of the residual error has to be supported.
#'
#' @param y univariat timeseries
#' @param window.lengths integer vector of length m containing the window lengths
#' @param step.lengths integer vector of length m containing the step length for each window length
#' @param pattern list object with x and y elements containing the coordinates for the pattern
#' @param window.normalization character defining the desired method to normalize the windows. Can be 'peak', 'begin', 'end', 'standard'
#' @param lmtype character defining what method to use for the linear regression. For now only 'robust' is possible
#' @param error.measure any function that takes a vector (residuals of the regression) as input and returns a scalar. 
#' 
#' @return data.frame where each row contains the features calculated on a certain window
patternSearch <- function(y,
                          window.lengths = floor(c(0.1, 0.3, 0.5) * length(y)),
                          step.lengths = rep(floor(c(0.1) * length(y)), length(window.lengths)),
                          pattern = list(x = c(0, 0.1, 1), y = c(0, 1, 0)),
                          window.normalization = 'peak',
                          lmtype = 'robust',
                          error.measure = log.lrvar.smooth.trim,
                          cluster.detection = 'mclust',
                          selection.criteria = 'robust.scale',
                          top.windows = 5
){
  lmfunc <- switch(lmtype,
                   robust = lmrobust.estimation)
  featuretable <- slidingFeatures(y = y,
                                  window.lengths = window.lengths,
                                  step.lengths = step.lengths,
                                  feature = function(x) {
                                    pattern.fitting.wayne(y = x, pat = pattern,
                                                         normalization = window.normalization,
                                                         lmfunc = lmfunc,
                                                         error.measure = error.measure)
                                  })
  
  featuretable <- featuretable %>% filter(!is.na(pattern.error) & !is.na(pattern.coef) & !is.na(base.error))
  
  # anomaly <- switch(cluster.detection,
  #                   none = TRUE,
  #                   mclust = mclustDist(X =  cbind(featuretable$pattern.coef,
  #                                                  featuretable$pattern.error)) > 10)
  
  if (TRUE) {
    featuretable[['pattern.error.scaled']] <- robustscaleQn(featuretable[['pattern.error']], 
                                                            idx = featuretable[['pattern.error']] < featuretable[['base.error']]*0.95)
    featuretable[['pattern.coef.scaled']] <- robustscaleQn(featuretable[['pattern.coef']],
                                                           mid = 0,
                                                           idx = featuretable[['pattern.error']] < featuretable[['base.error']]*0.95)
    eligible.windows <- featuretable %>%
      filter(pattern.error.scaled < 1.5) %>%
      arrange(-pattern.coef.scaled)
    if (nrow(eligible.windows) > 0) {
      out <- eligible.windows
      return(list(
        y = y,
        windowresults = out
        ))
    }
    else {return(NULL)}
  } else {
    return(NULL)
  }
}

#################################################################################################
#################################################################################################
#################################################################################################
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