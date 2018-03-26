vacf.cov <- function(X, lag=floor(dim(as.matrix(X))[1]/4), demean=TRUE){
  X <- as.matrix(X)
  T <- dim(X)[1]
  lags <- c(0:min(lag,T-1))
  if(demean){
    X <- apply(X,2,function(x) x - sum(x)/length(x))
  }
  vacf <- lapply(lags, function(l) (t(X[(1+l):T,])%*%X[1:(T-l),])/T)
  vacf <- array(unlist(vacf), dim=c(nrow(vacf[[1]]),ncol(vacf[[1]]),length(vacf)))
  dimnames(vacf)[[3]] <- as.list(lags)
  return(vacf)
}

weigths.lr <- function(M,T,kernel='Bartlett'){
  lags <- switch(kernel,
                 Truncated = if(floor(M)==0){NULL}else{1:min(floor(M),(T-1))},
                 Bartlett = if(floor(M)==0){NULL}else{1:min(floor(M),(T-1))},
                 Bohmann = if(floor(M)==0){NULL}else{1:min(floor(M),(T-1))},
                 Daniell = 1:(T-1),
                 Quadratic = 1:(T-1))
  k <- switch(kernel,
              Truncated = lags^0,
              Bartlett = 1 - abs(lags/M),
              Bohmann = (1-abs(lags/M))*cos(pi*(lags/M))+sin(pi*abs(lags/M))/pi,
              Daniell = sin(pi*lags/M)/(pi*lags/M),
              Quadratic = (25/(12*pi^2*(lags/M)^2))*((sin((6*pi*(lags/M))/5)/((6*pi*(lags/M))/5)-cos((6*pi*(lags/M))/5)))
  )
  return(list(weigths=c(1,k), upper=length(lags)))
}

vlongrun.var <- function(Y,M=4*(dim(as.matrix(Y))[1]/100)^(2/9), kernel='Bartlett', demean=TRUE){
  Y <- as.matrix(Y)
  T <- dim(Y)[1]
  w <- weigths.lr(M,T,kernel=kernel)
  if(w$upper == 0){
    omega <- vacf.cov(Y,lag=0,demean=demean)[,,1]
    delta <- omega
    sigma <- omega
  }
  else{
    k <- w$weigths
    vacf <- vacf.cov(Y, lag=w$upper, demean=demean)
    dimnames(vacf) <- NULL
    omega <- vacf[,,1]
    delta <- vacf[,,1]
    for(j in 2:(w$upper+1)){
      omega <- omega + k[j]*(vacf[,,j]+t(vacf[,,j]))
      delta <- delta + k[j]*vacf[,,j]
    }
    sigma <- vacf[,,1]
  }
  return(list(Sigma=as.matrix(sigma), Omega=as.matrix(omega), Delta=as.matrix(delta)))
}

PP.Z <- function(y,trend=c(0,1),kernel='Bartlett',M=(4*(length(y)/100)^(2/9)),detrend='one-step',test='coefficient'){
  y <- as.vector(y)
  T <- length(y)
  y.t <- y[2:T]
  y.lag <- y[1:(T-1)]
  if(is.null(trend)){
    y.t.detrend <- y.t
    y.lag.detrend <- y.lag
  }
  else{
    D <- sapply(trend, function(s) c(1:T)^s)
    if(detrend=='two-step'){
      theta <- solve(t(D)%*%D,t(D)%*%y)
      y.t.detrend <- y.t-D[2:T,]%*%theta
      y.lag.detrend <- y.lag - D[1:(T-1),]%*%theta
    }
    else{
      y.t.detrend <- y.t - D[2:T,]%*%solve(t(D[2:T,])%*%D[2:T,],t(D[2:T,])%*%y.t)
      y.lag.detrend <- y.lag - D[2:T,]%*%solve(t(D[2:T,])%*%D[2:T,],t(D[2:T,])%*%y.lag)
    }
  }
  rho <- sum(y.t.detrend*y.lag.detrend)/sum(y.lag.detrend^2)
  u.hat <- y.t.detrend-rho*y.lag.detrend
  sigma <- sum(u.hat^2)/(T-1)
  t <- (rho-1)*sqrt(sum(y.lag.detrend^2)/sigma)
  if(M == 'AND91'){
    M <- And.HAC91(u.hat,kernel=kernel)
  }
  omega <- vlongrun.var(u.hat, M=M, kernel=kernel, demean=FALSE)$Omega[1,1]
  if(test == 'tstatistic'){
    Z <- sqrt(sigma/omega)*t-0.5*(omega-sigma)/sqrt(omega*sum(y.lag.detrend^2)/(T^2))
  }
  else{
    Z <- T*(rho-1)-0.5*(omega-sigma)*(T^2)/sum(y.lag.detrend^2)
  }
  return(list(Z=as.vector(Z),M=as.vector(M)))
}

robustHarveyTest <- function(y, X) {
  n <- length(y)
  # trend coefficient for the I(0) - case
  y <- y - mean(y)
  X <- apply(X, 2, function(x) x - mean(x))
  beta.0 <- solve(t(X) %*% X, t(X) %*% y)[,1]
  u.0 <- y - X %*% beta.0
  omega.0 <- lrvar(u.0, type = 'Newey-West')
  s.0 <- sqrt(omega.0 * n * diag(solve(t(X) %*% X)))
  y.diff <- y[-1] - y[-n]
  X.diff <- X[-1,] - X[-n,]
  beta.1 <- solve(t(X.diff) %*% X.diff, t(X.diff) %*% y.diff)[,1]
  u.1 <- y.diff - X.diff %*% beta.1
  omega.1 <- lrvar(u.1, type = 'Newey-West')
  s.1 <- sqrt(omega.1 * (n - 1) * diag(solve(t(X.diff) %*% X.diff)))
  U <- PP.Z(y = u.0)$Z
  S <- sum(cumsum(u.0) ^ 2) / ((n ^ 2) * omega.0)
  lambda <- exp(-((U/S)^2))
  z <- (1 - lambda) * (beta.0 / s.0) + lambda * (beta.1 / s.1)
  return(
    list(z = z,
         p.value = pnorm(abs(z), mean = 0, sd = 1, lower.tail = FALSE),
         lambda = lambda,
         beta.0 = beta.0,
         s.0 = s.0,
         beta.1 = beta.1,
         s.1 = s.1)
  )
}