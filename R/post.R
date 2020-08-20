post <- function(tau_sq = 0.0001, g1 = 10001, g2 = 10000, Xcurr = X_curr, Xcand = X_cand, .Y = y, n = length(y)) {
  ts <- 1 / tau_sq
  V_curr <- solve(t(Xcurr) %*% Xcurr + tau_sq * diag(ncol(Xcurr)))
  V_cand <- solve(t(Xcand) %*% Xcand + tau_sq * diag(ncol(Xcand)))

  num <- 0.5 * determinant(ts * diag(ncol(V_curr)), logarithm = T)$mod + 0.5 * determinant(V_cand, logarithm = T)$mod
  den <- 0.5 * determinant(ts * diag(ncol(V_cand)), logarithm = T)$mod + 0.5 * determinant(V_curr, logarithm = T)$mod

  ahatcurr <- bhat(sig_sq = 1, X = Xcurr, y = y)
  ahatcand <- bhat(sig_sq = 1, X = Xcand, y = y)
  d_curr <- g2 + (t(.Y) %*% .Y) - (t(ahatcurr) %*% solve(V_curr) %*% ahatcurr)
  d_cand <- g2 + (t(.Y) %*% .Y) - (t(ahatcand) %*% solve(V_cand) %*% ahatcand)
  d_term <- (g1 + (n/2)) * (log(d_curr/d_cand))
  result <- (num - den) + d_term
  if(is.na(result)) browser()
  return(result) #simplifying everything to logs
} #RJMCMC posterior
