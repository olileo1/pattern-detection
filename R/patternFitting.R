patternFitting.orcutt <- function(y,
                                  pat = rep(1, length(y))) {
  n <- length(y)
  y <- robustScaling(y)
  X <- cbind(
    intercept = rep(1, n),
    slope = 1:n,
    pattern = approx(pat, n = length(y))$y
  )
  m.rob <- lmrob.fit(x = as.matrix(X),
                     y = y,
                     control = lmrob.control(setting = 'KS2011', maxit.scale = 500))
  y.smoothed <- y
  y.smoothed[which(abs(m.rob$residuals) > qnorm(0.975, mean = 0, sd = m.rob$scale))] <-
    m.rob$fitted.values[which(abs(m.rob$residuals) > qnorm(0.975, mean = 0, sd = m.rob$scale))]
  m.co.base <- cochraneOrcuttFast(y = y.smoothed, X = X[, c('intercept', 'slope')])
  m.co.pattern <- cochraneOrcuttFast(y = y.smoothed, X = X[, c('intercept', 'pattern')])
  return(list(
    pattern.coef = m.co.pattern$coefficients['pattern'],
    pattern.p-value = m.co.pattern$p.value['pattern'],
    f.statistic = m.co.pattern$co.residuals.var / m.co.base$co.residuals.var,
    f.pvalue = pf(m.co.pattern$co.residuals.var / m.co.base$co.residuals.var, df1 = n - 2, df2 = n - 2, lower.tail = FALSE)
  ))
}

patternFitting.roblm <- function(y,
                                  pat = rep(1, length(y))) {
  n <- length(y)
  y <- robustScaling(y)
  X <- cbind(
    intercept = rep(1, n),
    slope = 1:n,
    pattern = approx(pat, n = length(y))$y
  )
  m.rob <- lmrob.fit(x = as.matrix(X),
                     y = y,
                     control = lmrob.control(setting = 'KS2011', maxit.scale = 500))
  y.smoothed <- y
  y.smoothed[which(abs(m.rob$residuals) > qnorm(0.975, mean = 0, sd = m.rob$scale))] <-
    m.rob$fitted.values[which(abs(m.rob$residuals) > qnorm(0.975, mean = 0, sd = m.rob$scale))]
  m.co.base <- cochraneOrcuttFast(y = y.smoothed, X = X[, c('intercept', 'slope')])
  m.co.pattern <- cochraneOrcuttFast(y = y.smoothed, X = X[, c('intercept', 'pattern')])
  if (m.co.pattern$coefficients['pattern'] > 0 &
      m.co.pattern$p.value['pattern'] < 0.01) {
    return(m.co.pattern$co.residuals.var / m.co.base$co.residuals.var)
  } else {
    return(Inf)
  }
  return(out)
}
