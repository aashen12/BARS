p_sig <- function(a = 1, b = 1, X, beta,y) {
  n <- nrow(X)
  a_term <- a + (n/2)
  b_term <- 0.5 * (2*b + (t(y - (X %*% beta)) %*% (y - (X %*% beta))))
  1 / rgamma(1, shape = a_term, rate = b_term)
} # full conditional for sigma^2 in the regression equation

