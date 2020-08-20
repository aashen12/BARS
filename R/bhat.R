bhat <- function(sig_sq, tau_sq = 10000, X, y) {
  p = ncol(X) - 1
  sig <- solve( (1/sig_sq) * (t(X) %*% X) + (1/tau_sq) * diag(p+1) )
  mu <- (1/sig_sq) * sig %*% t(X) %*% y
  mu
} #marginalized betahat for regression coefficients