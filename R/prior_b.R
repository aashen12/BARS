prior_b <- function(k) {
  return(0.5 * (p(.k = k+1) - p(.k = k)))
} #prior in RJMCMC ratio for a birth step
