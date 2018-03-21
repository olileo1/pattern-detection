patternFitting.orcutt <- function(y,
                                  pat = rep(1, length(y)),
                                  max.iter = 100,
                                  yule = TRUE) {
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
  m.co.base <- cochraneOrcuttFast(y = y.smoothed, X = X[, c('intercept', 'slope')], max.iter = max.iter, yule = yule)
  m.co.pattern <- cochraneOrcuttFast(y = y.smoothed, X = X[, c('intercept', 'slope', 'pattern')], max.iter = max.iter, yule = yule)
  return(list(
    pattern.coef = m.co.pattern$coefficients['pattern'],
    pattern.p.value = m.co.pattern$p.value['pattern'],
    f.statistic = m.co.base$co.residuals.var / m.co.pattern$co.residuals.var,
    f.pvalue = pf(m.co.base$co.residuals.var / m.co.pattern$co.residuals.var, df1 = n - 2, df2 = n - 2, lower.tail = FALSE)
  ))
}

patternFitting.roblm <- function(y,
                                 pat = rep(1, length(y)),
                                 alpha = 0.1) {
  n <- length(y)
  y <- robustScaling(y)
  X <- cbind(
    intercept = rep(1, n),
    slope = 1:n,
    pattern = approx(pat, n = length(y))$y
  )
  m1 <- lmrob.fit(x = as.matrix(X[, c('intercept', 'slope')]),
                  y = y,
                  control = lmrob.control(setting = 'KS2011', maxit.scale = 500))
  m2 <- lmrob.fit(x = as.matrix(X[, c('intercept', 'pattern')]),
                  y = y,
                  control = lmrob.control(setting = 'KS2011', maxit.scale = 500))
  return(list(
    pattern.coef = m2$coefficients['pattern'],
    pattern.p.value = pt(abs(m2$coefficients['pattern']) / m2$cov['pattern', 'pattern'], df = m2$degree.freedom, lower.tail = FALSE) * 2,
    f.statistic = m1$scale / m2$scale,
    f.pvalue = pf(m1$scale / m2$scale, df1 = m1$df.residual, df2 = m2$df.residual, lower.tail = FALSE)
  ))
  return(out)
}
