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

patternFitting.roblm2 <- function(y,
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
                  control = lmrob.control(setting = 'KS2011', maxit.scale = 1000, max.it = 1000, refine.tol = 1e-6, rel.tol = 1e-6, solve.tol = 1e-6))
  m2 <- lmrob.fit(x = as.matrix(X[, c('intercept', 'pattern')]),
                  y = y,
                  control = lmrob.control(setting = 'KS2011', maxit.scale = 1000, max.it = 1000, refine.tol = 1e-6, rel.tol = 1e-6, solve.tol = 1e-6))
  m3 <- lmrob.fit(x = as.matrix(X[, c('intercept', 'slope', 'pattern')]),
                  y = y,
                  control = lmrob.control(setting = 'KS2011', maxit.scale = 1000, max.it = 1000, refine.tol = 1e-6, rel.tol = 1e-6, solve.tol = 1e-6))
  y.smoothed <- y
  y.smoothed[which(abs(m3$residuals) > qnorm(0.975, mean = 0, sd = m3$scale))] <-
    m3$fitted.values[which(abs(m3$residuals) > qnorm(0.975, mean = 0, sd = m3$scale))]
  res <- robustHarveyTest(y = y.smoothed,
                          X = X[, c('slope', 'pattern')])
  return(list(
    pattern_coef_m2 = m2$coefficients['pattern'],
    pattern_coef_m3 = m3$coefficients['pattern'],
    m1_scale = m1$scale,
    m2_scale = m2$scale,
    m3_scale = m3$scale,
    m1_lrvar = lrvar(m1$residuals, type = 'Newey-West'),
    m2_lrvar = lrvar(m2$residuals, type = 'Newey-West'),
    m3_lrvar = lrvar(m3$residuals, type = 'Newey-West'),
    harvey_z = res$z['pattern'],
    harvey_p = res$p.value['pattern']
  ))
  return(out)
}

patternFitting.roblm3 <- function(y,
                                  pat = rep(1, length(y)),
                                  alpha = 0.1) {
  n <- length(y)
  y <- y / y[1]
  X <- cbind(
    intercept = rep(1, n),
    slope = 1:n,
    pattern = approx(pat, n = length(y))$y
  )
  m1 <- lmrob.fit(x = as.matrix(X[, c('intercept', 'slope')]),
                  y = y,
                  control = lmrob.control(setting = 'KS2011', maxit.scale = 1000, max.it = 1000, refine.tol = 1e-6, rel.tol = 1e-6, solve.tol = 1e-6))
  m2 <- lmrob.fit(x = as.matrix(X[, c('intercept', 'pattern')]),
                  y = y,
                  control = lmrob.control(setting = 'KS2011', maxit.scale = 1000, max.it = 1000, refine.tol = 1e-6, rel.tol = 1e-6, solve.tol = 1e-6))
  m3 <- lmrob.fit(x = as.matrix(X[, c('intercept', 'slope', 'pattern')]),
                  y = y,
                  control = lmrob.control(setting = 'KS2011', maxit.scale = 1000, max.it = 1000, refine.tol = 1e-6, rel.tol = 1e-6, solve.tol = 1e-6))
  y.smoothed <- y
  y.smoothed[which(abs(m3$residuals) > qnorm(0.975, mean = 0, sd = m3$scale))] <-
    m3$fitted.values[which(abs(m3$residuals) > qnorm(0.975, mean = 0, sd = m3$scale))]
  res <- robustHarveyTest(y = y.smoothed,
                          X = X[, c('slope', 'pattern')])
  return(list(
    slope_coef_m1 = m1$coefficients['slope'],
    pattern_coef_m2 = m2$coefficients['pattern'],
    slope_coef_m3 = m3$coefficients['slope'],
    pattern_coef_m3 = m3$coefficients['pattern'],
    m1_scale = m1$scale,
    m2_scale = m2$scale,
    m3_scale = m3$scale,
    m1_lrvar = lrvar(m1$residuals, type = 'Newey-West'),
    m2_lrvar = lrvar(m2$residuals, type = 'Newey-West'),
    m3_lrvar = lrvar(m3$residuals, type = 'Newey-West'),
    harvey_slope.0 = res$beta.0['slope'],
    harvey_pattern.0 = res$beta.0['pattern'],
    harvey_slope.1 = res$beta.1['slope'],
    harvey_pattern.1 = res$beta.1['pattern'],
    harvey_s_slope.0 = res$s.0['slope'],
    harvey_s_pattern.0 = res$s.0['pattern'],
    harvey_s_slope.1 = res$s.1['slope'],
    harvey_s_pattern.1 = res$s.1['pattern'],
    harvey_lambda = res$lambda
  ))
  return(out)
}

patternFitting.roblm4 <- function(y,
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
                  control = lmrob.control(setting = 'KS2011', maxit.scale = 1000, max.it = 1000, refine.tol = 1e-6, rel.tol = 1e-6, solve.tol = 1e-6))
  m2 <- lmrob.fit(x = as.matrix(X[, c('intercept', 'pattern')]),
                  y = y,
                  control = lmrob.control(setting = 'KS2011', maxit.scale = 1000, max.it = 1000, refine.tol = 1e-6, rel.tol = 1e-6, solve.tol = 1e-6))
  m3 <- lmrob.fit(x = as.matrix(X[, c('intercept', 'slope', 'pattern')]),
                  y = y,
                  control = lmrob.control(setting = 'KS2011', maxit.scale = 1000, max.it = 1000, refine.tol = 1e-6, rel.tol = 1e-6, solve.tol = 1e-6))
  y.smoothed <- y
  y.smoothed[which(abs(m3$residuals) > qnorm(0.975, mean = 0, sd = m3$scale))] <-
    m3$fitted.values[which(abs(m3$residuals) > qnorm(0.975, mean = 0, sd = m3$scale))]
  res <- robustHarveyTest(y = y.smoothed,
                          X = X[, c('slope', 'pattern')])
  return(list(
    slope_coef_m1 = m1$coefficients['slope'],
    pattern_coef_m2 = m2$coefficients['pattern'],
    slope_coef_m3 = m3$coefficients['slope'],
    pattern_coef_m3 = m3$coefficients['pattern'],
    m1_scale = m1$scale,
    m2_scale = m2$scale,
    m3_scale = m3$scale,
    m1_lrvar = lrvar(m1$residuals, type = 'Newey-West'),
    m2_lrvar = lrvar(m2$residuals, type = 'Newey-West'),
    m3_lrvar = lrvar(m3$residuals, type = 'Newey-West'),
    harvey_slope.0 = res$beta.0['slope'],
    harvey_pattern.0 = res$beta.0['pattern'],
    harvey_slope.1 = res$beta.1['slope'],
    harvey_pattern.1 = res$beta.1['pattern'],
    harvey_s_slope.0 = res$s.0['slope'],
    harvey_s_pattern.0 = res$s.0['pattern'],
    harvey_s_slope.1 = res$s.1['slope'],
    harvey_s_pattern.1 = res$s.1['pattern'],
    harvey_lambda = res$lambda
  ))
  return(out)
}

convoluteTimeseries <- function(y, length.out) {
  smoothed <- sapply(split(y, ceiling(seq_along(y) * length.out / length(y))), function(x) median(x))
  smoothed <- (smoothed - min(smoothed)) / (max(smoothed) - min(smoothed))
  return(
    list(
      smoothed = smoothed
    )
  )
}

movingMedian <- function(y, n.back = 2, n.ahead = 2) {
  n <- length(y)
  y.smoothed <- sapply(1:n, function(i) {
    median(y[max(0, (i - n.back)):min((i + n.ahead), n)])
  })
  return(y.smoothed)
}

patternFitting.harvey <- function(y,
                                  pat = list(x = c(0, 0.1, 1), y = c(0, 1, 0))) {
  n <- length(y)
  y.smoothed <- movingMedian(y = y, n.back = 2, n.ahead = 2)
  y.smoothed <- y.smoothed / median(y.smoothed[1:floor(n / 20)])
  X <- cbind(
    intercept = rep(1, n),
    slope = 1:n,
    pattern = approx(pat, n = length(y))$y
  )
  harvey1 <- robustHarveyTest(y = y.smoothed,
                              X = as.matrix(X[, c('slope')]))
  harvey2 <- robustHarveyTest(y = y.smoothed,
                              X = as.matrix(X[, c('pattern')]))
  harvey3 <- robustHarveyTest(y = y.smoothed,
                              X = as.matrix(X[, c('pattern', 'slope')]))
  return(list(
    harvey1.var = harvey1$lambda * harvey1$var.0 + (1 - harvey1$lambda) * harvey1$var.1,
    harvey2.var = harvey2$lambda * harvey2$var.0 + (1 - harvey2$lambda) * harvey2$var.1,
    harvey3.var = harvey3$lambda * harvey3$var.0 + (1 - harvey3$lambda) * harvey3$var.1,
    harvey2.beta0 = harvey2$beta.0,
    harvey2.s0 = harvey2$s.0,
    harvey2.beta1 = harvey2$beta.1,
    harvey2.s1 = harvey2$s.1,
    harvey2.lambda = harvey2$lambda,
    harvey3.beta0 = harvey3$beta.0['pattern'],
    harvey3.s0 = harvey3$s.0['pattern'],
    harvey3.beta1 = harvey3$beta.1['pattern'],
    harvey3.s1 = harvey3$s.1['pattern'],
    harvey3.lambda = harvey3$lambda
  ))
}

patternFitting.roblm.new <- function(y,
                                     pat = list(x = c(0, 0.1, 1), y = c(0, 1, 0))) {
  n <- length(y)
  y <- y / median(y[1:floor(n / 20)])
  X <- cbind(
    intercept = rep(1, n),
    slope = 1:n,
    pattern = approx(pat, n = length(y))$y
  )
  m1 <- lmrob.fit(x = as.matrix(X[, c('intercept', 'slope')]),
                  y = y,
                  control = lmrob.control(setting = 'KS2011', maxit.scale = 1000, max.it = 1000, refine.tol = 1e-6, rel.tol = 1e-6, solve.tol = 1e-6))
  m2 <- lmrob.fit(x = as.matrix(X[, c('intercept', 'pattern')]),
                  y = y,
                  control = lmrob.control(setting = 'KS2011', maxit.scale = 1000, max.it = 1000, refine.tol = 1e-6, rel.tol = 1e-6, solve.tol = 1e-6))
  m3 <- lmrob.fit(x = as.matrix(X[, c('intercept', 'slope', 'pattern')]),
                  y = y,
                  control = lmrob.control(setting = 'KS2011', maxit.scale = 1000, max.it = 1000, refine.tol = 1e-6, rel.tol = 1e-6, solve.tol = 1e-6))
  return(list(
    m1.scale = m1$scale,
    m2.scale = m2$scale,
    m3.scale = m3$scale,
    m2.beta = m2$coefficients['pattern'],
    m2.s = m2$cov['pattern', 'pattern'],
    m3.beta = m3$coefficients['pattern'],
    m3.s = m3$cov['pattern', 'pattern']
  ))
}

patternFitting.partyparty <- function(y,
                                      pat = list(
                                        x = c(0, 0.1, 1),
                                        y = c(0, 0.1, 1)
                                        )
                                      ) {
  n <- length(y)
  y <- y / median(y[1:floor(n / 20)])
  X <- cbind(
    intercept = rep(1, n),
    slope = 1:n,
    pattern = approx(pat, n = length(y))$y
  )
  m.base <- lmrob.fit(x = as.matrix(X[, c('intercept', 'slope')]),
                      y = y,
                      control = lmrob.control(setting = 'KS2011', maxit.scale = 1000, max.it = 1000, k.max = 1000, refine.tol = 1e-6, rel.tol = 1e-6, solve.tol = 1e-6))
  m1 <- lmrob.fit(x = as.matrix(X[, c('intercept', 'pattern')]),
                  y = y,
                  control = lmrob.control(setting = 'KS2011', maxit.scale = 1000, max.it = 1000, k.max = 1000, refine.tol = 1e-6, rel.tol = 1e-6, solve.tol = 1e-6))
  m2 <- lmrob.fit(x = as.matrix(X[, c('intercept', 'slope', 'pattern')]),
                  y = y,
                  control = lmrob.control(setting = 'KS2011', maxit.scale = 1000, max.it = 1000, k.max = 1000, refine.tol = 1e-6, rel.tol = 1e-6, solve.tol = 1e-6))
  return(
    list(
      m.base.scale = m.base$scale,
      m.base.slope.coef = m.base$coefficients['slope'],
      m1.scale = m1$scale,
      m1.pat.coef = m1$coefficients['pattern'],
      m1.pat.coef.sd = m1$cov['pattern', 'pattern'],
      m1.max.sign.pat.coef = m1$coefficients['pattern'] -
        qt(0.99, df = m1$df.residual, lower.tail = TRUE) * sqrt(m1$cov['pattern', 'pattern']),
      m1.res = list(res = m1$residuals),
      m1.fit = list(fit = m1$fitted.values),
      m2.scale = m2$scale,
      m2.pat.coef = m2$coefficients['pattern'],
      m2.pat.coef.sd = m2$cov['pattern', 'pattern'],
      m2.max.sign.pat.coef = m2$coefficients['pattern'] -
        qt(0.99, df = m2$df.residual, lower.tail = TRUE) * sqrt(m2$cov['pattern', 'pattern']),
      m2.slope.coef = m2$coefficients['slope'],
      m2.res = list(res = m2$residuals),
      m2.fit = list(fit = m2$fitted.values)
    )
  )
}
