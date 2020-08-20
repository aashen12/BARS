p <- function(lambda = 1, .k) {
  num <- .k * log(lambda)
  den <- log(exp(lambda) - 1) + lfactorial(.k)
  return(num - den)
} #RJMCMC poisson prior with default lambda = 1
