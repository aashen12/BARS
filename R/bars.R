# Documentation goes here.
# GAUSSIAN LIKELIHOOD
library(mvtnorm)
bars <- function(its, max_knot = 200, verbose = FALSE) {
  # Function that returns the results of an RJMCMC algorithm to fit
  # a linear basis spline to a univariate dataset.
  # Input of the function is simply the number of MCMC iterations, with
  # additional functions within the function.
  # The RJMCMC acceptance probabilities are evaluated on the log-scale for convenience


  ### RJ MCMC Conditionals ###

  ## Initializing Values for RJMCMC ##

  nknot <- 0 # First iteration must be a birth
  X_curr <- as.matrix(rep(1, length(x)))
  # first current X-matrix is just an intercept (matrix of 1's)

  mat_t <- matrix(NA, its, max_knot)
  mat_s <- matrix(NA, its, max_knot)

  mat_beta <- matrix(NA, its, max_knot)
  mat_beta[1,1] <- mean(y) #arbitrary starting value for beta (in this case beta_0)
  mat_sig <- rep(NA, its)
  mat_sig[1] <- 1

  for(it in 2:its) {

    # need to generate number of knots and produce X_cand
    # you must choose birth in first iteration
    # sample knot, sign
    # if accept, you have 1 BF
    # next step: B, D, or C
    choice <- samp(knots = nknot, max = max_knot)

    if(choice == 1) { # BIRTH
      candidate_t <- runif(1)
      candidate_s <- sample(c(-1,1), 1, replace = TRUE)

      basis_vec <- pos(candidate_s * (x - candidate_t))

      X_cand <- cbind(X_curr, basis_vec)
      #candidate X matrix is current binded with candidate basis vector

      ratio_rj <- post(Xcurr = X_curr, Xcand = X_cand) + prior_b(k = nknot) + proposal_b(k = nknot)

      accept_prob <- min(0, ratio_rj)
      if(is.na(ratio_rj)) browser()
      if(log(runif(1)) < accept_prob) {
        nknot <- nknot + 1
        mat_t[it,] <- mat_t[it-1,]
        mat_t[it,nknot] <- candidate_t
        mat_s[it,] <- mat_s[it-1,]
        mat_s[it,nknot] <- candidate_s
        X_curr <- X_cand
      } else {
        mat_t[it,] <- mat_t[it-1,]
        mat_s[it,] <- mat_s[it-1,]
      }
    }

    else if(choice == 2) { #DEATH
      pick <- sample(2:(nknot+1), 1) #pick = (column in X to delete) = knot number + 1 (due to intercept)
      X_cand <- X_curr[,-pick]

      ratio_rj <- post(Xcurr = X_curr, Xcand = X_cand) + prior_d(k = nknot) + proposal_d(k = nknot)

      accept_prob <- min(0, ratio_rj)
      if(is.na(ratio_rj)) browser()

      if(log(runif(1)) < accept_prob) {
        nknot <- nknot - 1
        X_curr <- X_cand
        mat_t[it,(1:nknot)] <- mat_t[it-1,(1:(nknot+1))[-(pick-1)]]
        # (pick - 1) = (column in mat_t to delete) = (knot number to delete)

        mat_s[it,(1:nknot)] <- mat_s[it-1,(1:(nknot+1))[-(pick-1)]]
      } else{
        mat_t[it,] <- mat_t[it-1,]
        mat_s[it,] <- mat_s[it-1,]
      }
    }

    else { #CHANGE
      X_cand <- X_curr
      pick <- sample(2:(nknot+1), 1)
      #pick = (column in X to change) = knot number + 1 (due to intercept)
      candidate_t <- runif(1)
      candidate_s <- sample(c(-1,1), 1, replace = TRUE)

      basis_vec <- pos(candidate_s * (x - candidate_t))
      X_cand[,pick] <- basis_vec

      ratio_change <- function(tau_sq = 0.01, g1 = 0.01, g2 = 0.01, n = length(y)) { # MAY NEED TO BE MODIFIED FOR BMARS
        Vprime <- solve(t(X_cand) %*% X_cand + tau_sq * diag(ncol(X_cand)))
        V <- solve(t(X_curr) %*% X_curr + tau_sq * diag(ncol(X_curr)))

        aprime <- Vprime %*% t(X_cand) %*% y
        a <- V %*% t(X_curr) %*% y

        dprime <- g2 + (t(y) %*% y) - (t(aprime) %*% solve(Vprime) %*% aprime)
        d <- g2 + (t(y) %*% y) - (t(a) %*% solve(V) %*% a)

        V_part <- 0.5 * determinant(Vprime)$mod -  0.5 * determinant(V)$mod
        d_part <- (g1 + n/2) * (log(d) - log(dprime))

        V_part + d_part
      }

      acc <- min(0, ratio_change())

      if(log(runif(1)) < acc) {
        X_curr <- X_cand
        mat_t[it,] <- mat_t[it-1,]
        mat_t[it, (pick-1)] <- candidate_t
        mat_s[it,] <- mat_s[it-1,]
        mat_s[it, (pick-1)] <- candidate_s
        # (pick - 1) = (column in mat_t to change) = (knot number to change)
      } else{
        mat_t[it,] <- mat_t[it-1,]
        mat_s[it,] <- mat_s[it-1,]
      }
    }

    mat_beta[it,1:(nknot+1)] <- p_beta(sig_sq = mat_sig[it-1], X = X_curr,y=y)
    mat_sig[it] <- p_sig(X = X_curr, beta = mat_beta[it,(1:(nknot+1))],y=y)

    if(verbose == TRUE) {
      if(it %% 1000 == 0) {
        cat("Iteration", it, "\t")
        cat("sigma^2 =",mat_sig[it],"\n")
      }
    }
  }
  list("beta"=mat_beta, "sig"=mat_sig, "knots"=mat_t, "signs"=mat_s, "X"=X_curr,
       "knots_total" = nknot)
}
