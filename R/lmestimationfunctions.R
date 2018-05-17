lmrobust.estimation <- function(y, X) {
  tryCatch({
    m <- lmrob.fit(x = X,
                   y = y,
                   control = lmrob.control(setting = 'KS2011',
                                           maxit.scale = 10000,
                                           max.it = 10000,
                                           k.max = 10000,
                                           refine.tol = 1e-6,
                                           rel.tol = 1e-6,
                                           solve.tol = 1e-6))
    out$coefficients <- m$coefficients
    out$residuals <- m$residuals
  },
  warning = function(w) {
    coefficients <- solve(t(X) %*% X, t(X) %*% y)[,1]
    out$residuals <- (y - X %*% coefficients)[,1]
    out$coefficients <- coefficients
  })
  return(out)
}