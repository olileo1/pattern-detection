source('dependencies.R')

file.dir <- 'R'
file.sources <- list.files(file.dir, pattern = '*.R')
sapply(paste(file.dir, file.sources, sep = '/'), source)
