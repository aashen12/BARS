proposal_b <- function(k) {
  d <- 1/3 #probability of death
  b <- 1/3 #probability of birth
  # NOTE: P(change) = 1/3
  num <- log(d) - log(k)
  den <- log(b) - log(2)
  return(num - den)
} #proposal in RJMCMC ratio for a birth step