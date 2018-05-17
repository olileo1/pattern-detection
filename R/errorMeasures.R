lrvar.smooth <- function(x) {lrvar(movingMedian(x, n.back = 5, n.ahead = 5))}
lrvar.smooth.trim <- function(x){lrvar(mid(movingMedian(x, n.back = 5, n.ahead = 5), alpha.begin = 0.1, alpha.end = 0.1))}
log.lrvar.smooth <- function(x) {log(lrvar(movingMedian(x, n.back = 5, n.ahead = 5)))}
log.lrvar.smooth.trim <- function(x){log(lrvar(mid(movingMedian(x, n.back = 5, n.ahead = 5), alpha.begin = 0.1, alpha.end = 0.1)))}