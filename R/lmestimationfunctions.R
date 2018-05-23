lmrobust.estimation <- function(y, X) {
  out <- tryCatch({
    m <- lmrob.fit(x = X,
                   y = y,
                   control = lmrob.control(setting = 'KS2011',
                                           maxit.scale = 10000,
                                           max.it = 10000,
                                           k.max = 10000,
                                           refine.tol = 1e-6,
                                           rel.tol = 1e-6,
                                           solve.tol = 1e-6))
    return(list(
      coefficients = m$coefficients,
      residuals = m$residuals
    ))
  },
  warning = function(w) {
    coefficients <- solve(t(X) %*% X, t(X) %*% y)[,1]
    return(list(
      coefficients = coefficients,
      residuals = (y - X %*% coefficients)[,1]
    ))
  },
  error = function(e) {
    coefficients <- solve(t(X) %*% X, t(X) %*% y)[,1]
    return(list(
      coefficients = coefficients,
      residuals = (y - X %*% coefficients)[,1]
    ))
  })
  return(out)
}
