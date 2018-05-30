# Simualtion to test if our method works for simulated data
## Function to create a pattern with specific length and height
create.pattern <- function(points = cbind(x = c(0, 0.5, 1), y = c(0, 1, 0)),
                           n = 100,
                           height = 100) {
  out <- approx(points, n = n)$y * height
  return(out)
}

## function to create the output of a polynom for given start value and end value and given polynomcoefficients
## the first two coefficients are irrelevant, since we need two degrees of freedom to fit the polynom to the start and end condition
create.poly <- function(start = 0,
                        end = 100,
                        n = 100,
                        p.coef = c(1,1,1)) {
  if (length(p.coef) < 2) return(rep(NA, n))
  if (length(p.coef) == 2) {
    s1 <- 0
    s2 <- 0
  } else {
    s1 <- sum(p.coef[-c(1, 2)] * (-1) ^ (2:(length(p.coef) - 1)))
    s2 <- sum(p.coef[-c(1, 2)])
  }
  coef.out <- c((start + end - s1 -s2) / 2, (-start + end + s1 -s2) / 2, p.coef[-c(1, 2)])
  x <- seq(-1, 1, length.out = n)
  X <- do.call('cbind', lapply(0:(length(p.coef) - 1), function(p) x^p))
  return((X %*% coef.out)[,1])
}

## function that creates a complete stream:
## Starting with a polynomial like function, then comes a pattern and then again a polynomial like function
create.stream <- function(pre.start,
                          pre.coef,
                          pre.n,
                          pattern.start,
                          pattern.points,
                          pattern.n,
                          pattern.height,
                          past.start,
                          past.coef,
                          past.n,
                          past.end) {
  pre <- create.poly(start = pre.start,
                     end = pattern.start,
                     n = pre.n,
                     p.coef = pre.coef)
  pattern.poly <- create.poly(start = pattern.start,
                              end = past.start,
                              n = pattern.n,
                              p.coef = c(0,0))
  pattern <- create.pattern(points = pattern.points,
                            n = pattern.n,
                            height = pattern.height)
  past <- create.poly(start = past.start,
                      end = past.end,
                      n = past.n,
                      p.coef = past.coef)
  return(c(pre, pattern + pattern.poly, past))
}

# function to create noise with spikes, it is also possible to create the noise as randomwalk
noise_generator <- function(n = 50, sd = 1, spike.prob = 0.1, ar1 = 0.8, season = NULL) {
  if (!is.null(season)) {
    noise <- sin(2 * pi * (1:n) / season) * sd * 0.8 + rnorm(n, mean = 0, sd = sd * 0.2)
  } else {
    innovation_sd <- sd * sqrt(1 - ar1 ^ 2)
    noise <- rnorm(1, mean = 0, sd = innovation_sd)
    for(i in 2:n) {
      noise[i] <- ar1 * noise[(i - 1)] + rnorm(1, mean = 0, sd = innovation_sd)
    }
  }
  spikes <- rbinom(n, size = 1, prob = spike.prob) * 6 * sd
  return(noise + spikes)
}

create.randompattern <- function(type = 'type0',
                                 magnitude = 1,
                                 length = 1000,
                                 plot = FALSE,
                                 rsd = 0.1,
                                 ar1 = 0.1,
                                 spike.prob = 0.01) {
  pattern.points <- switch(type,
                           type0 = cbind(x = c(0, 1), y = c(0, 0)),
                           type1 = cbind(x = c(0, 0.1, 1), y = c(0, 1, 0)),
                           type2 = cbind(x = c(0, 0.2, 1), y = c(0, 1, 0)),
                           type3 = cbind(x = c(0, 0.5, 1), y = c(0, 1, 0)),
                           type4 = cbind(x = c(0, 0.9, 1), y = c(0, 1, 0)))
  start.height <- runif(n = 1, min = 0.9 * magnitude, max = 1.1 * magnitude)
  # total length of the timeseries is random
  total.n <- runif(n = 1, min = 0.8 * length, max = 1.2 * length)
  splits <- runif(4, 0.1, 0.9)
  splits <- c(splits[1], sum(splits[2:3]), splits[4])
  splits <- splits / sum(splits)
  pre.growth <- 0
  pattern.growth <- 0
  past.growth <- 0
  pre.coef <- c(1, 1)
  pre.n <- floor(splits[1] * total.n)
  past.coef <- c(1, 1)
  past.n <- floor(splits[3] * total.n)
  pattern.n <- floor(splits[2] * total.n)
  pattern.height.ratio <- runif(1, 0.1, 0.4)
  pattern.height <- start.height * pattern.height.ratio
  
  x <- create.stream(pre.start = start.height,
                     pre.coef = pre.coef,
                     pre.n = pre.n,
                     pattern.start = start.height * (1 + pre.growth),
                     pattern.points = pattern.points,
                     pattern.n = pattern.n,
                     pattern.height = pattern.height,
                     past.start = start.height * (1 + pre.growth) * (1 + pattern.growth),
                     past.coef = past.coef,
                     past.n = past.n,
                     past.end = start.height * (1 + pre.growth) * (1 + pattern.growth) * (1 + past.growth))
  
  rsd <- rsd
  #stn <- c(0.5, 1)
  ar1 <- ar1
  #ar1 <- c(0.1, 0.9)
  spike.prob <- spike.prob
  var.x <- var(seq(magnitude, magnitude * (1.4), length.out = total.n))
  params <- expand.grid(rsd = rsd, ar1 = ar1, spike.prob = spike.prob)
  series <- lapply(1:dim(params)[1], function(i) {
    noise <- noise_generator(n = length(x), sd = params[i, 'rsd'] / qnorm(0.95, 0, 1), spike.prob = params[i, 'spike.prob'],
                             season = NULL, ar1 = params[i, 'ar1'])
    list(
      y = x + noise,
      rsd = params[i, 'rsd'],
      spike.prob = params[i, 'spike.prob'],
      ar1 = params[i, 'ar1']
    )
  })
  if (plot) {
    plot.data <- do.call('rbind', lapply(c1, series, function(serie) {
      if (serie$spike.prob == 0.01) {
        data.frame(x = 1:length(serie$y), y = serie$y, rsd = serie$rsd, ar1 = serie$ar1)
      } else {
        NULL
      }
    }))
    print(ggplot(data = plot.data, aes(x = x, y = y)) +
            geom_rect(aes(xmin = pre.n + 1, xmax = pre.n + pattern.n + 1, ymin = -Inf, ymax = Inf), color = 'blue', alpha = 0.5) +
            geom_line() + facet_grid(rsd ~ ar1) + ggtitle(type))
  }
  return(list(series = series,
              type = type,
              pattern.window = c(pre.n + 1, pre.n + pattern.n + 1)))
}