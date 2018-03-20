cochraneOrcuttFast <- function(y, X, convergence = 8) {
  n <- length(y)
  e <- y - X %*% solve(t(X) %*% X, t(X) %*% y)
  rho <- sum(e[-1] * e[-n]) / sum(e[-1] ^ 2)
  rho.hist <- rho
  X.trans <- X[-1, ] - rho * X[-n, ]
  y.trans <- y[-1] - rho * y[-n]
  e <- y - X %*% solve(t(X.trans) %*% X.trans, t(X.trans) %*% y.trans)
  rho <- sum(e[-1] * e[-n]) / sum(e[-1] ^ 2)
  rho.hist[2] <- rho
  i <- 2
  while(round(rho.hist[i - 1], convergence) != round(rho.hist[i], convergence)) {
    X.trans <- X[-1, ] - rho * X[-n, ]
    y.trans <- y[-1] - rho * y[-n]
    beta <- solve(t(X.trans) %*% X.trans, t(X.trans) %*% y.trans)
    e <- y - X %*% beta
    rho <- sum(e[-1] * e[-n]) / sum(e[-1] ^ 2)
    i <- i + 1
    rho.hist[i] <- rho
  }
  out <- list()
  out$fitted.values <- (X %*% beta)[,1]
  out$residuals <- y - out$fitted.values
  out$co.fitted.values <- (X.trans %*% beta)[,1]
  out$co.residuals <- y.trans - out$co.fitted.values
  out$coefficients <- beta[,1]
  out$co.residuals.var <- sum(out$co.residuals ^ 2) / (n-2)
  out$rho <- rho.hist[(i - 1)]
  out$t.value <- (beta / sqrt(out$co.residuals.var * diag(solve(t(X.trans) %*% X.trans))))[,1]
  out$p.value <- pt(abs(out$t.value), df = n - 1 - length(beta), lower.tail = FALSE) * 2
  return(out)
}
