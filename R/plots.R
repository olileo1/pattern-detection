sliding.features.plot.best.windows <- function(out, top.n = 5) {
  p <- lapply(1:min(nrow(out), top.n), function(i) {
    ggplot(data = data.frame(x = 1:length(out$y), y = out$y), aes(x = x, y = y)) +
      geom_rect(aes(xmin = out$windowresults[['win.start']][i],
                    xmax = out$windowresults[['win.end']][i],
                    ymin = -Inf, ymax = Inf)) +
      geom_line()
  })
  do.call('grid.arrange', c(p, ncol = 1))
}

sliding.features.scatter.plot <- function(out) {
  out$windowresults %>% ggplot(aes(x = pattern.coef.scaled, y = pattern.error.scaled)) + geom_point()
}
