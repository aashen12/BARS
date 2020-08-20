p_t <- function(sig_sq, betahat, X, y) {
  (-1/(2*sig_sq)) * t(y - (X %*% betahat)) %*% (y - (X %*% betahat))
} # full conditional for the knot locations, t
