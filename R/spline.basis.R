spline.basis <- function(nknot, knots, signs) {
  
  Xm <- matrix(NA, nknot, length(x))
  for(i in 1:nknot) {
    for(j in 1:length(x)) {
      Xm[i,j] <- pos(signs[i] * (x[j] - knots[i]))
    } 
  } 
  Xm <- t(Xm)
  Xm <- cbind(1, Xm)
  Xm
}