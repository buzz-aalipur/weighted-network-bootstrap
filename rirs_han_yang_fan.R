## =======================================================
## rirs_han_yang_fan.R
##
## Implements the RIRS test from Han, Yang, & Fan (2024)
## "Universal Rank Inference via Residual Subsampling..."
## Paper: https://arxiv.org/abs/1912.11583
## =======================================================

source("helpers.R") # For .as_sym_no_diag and .top_k0_approx

#' Implements the Han, Yang, & Fan (2024) RIRS test
#'
#' @param A The n x n adjacency matrix.
#' @param K0 The rank of the null hypothesis H0: K = K0.
#' @param alpha The significance level.
#' @param m The subsampling parameter. Per Section 4.1.1,
#'   P(Y_ij=1) = 1/m. Defaults to m = sqrt(n).
#' @return A list with the test statistic, p-value, and reject decision.
rirs_hyf_test_rank <- function(A, K0, alpha = 0.05, m = NULL) {
  
  # 1. Symmetrize, remove diag
  A <- .as_sym_no_diag(A)
  n <- nrow(A)
  
  # 2. Set default subsampling parameter m (from paper Sec 4.1.1)
  if (is.null(m)) {
    m <- sqrt(n)
  }
  
  # 3. Calculate residual matrix R (called W_hat in paper)
  A_k0 <- if (K0 > 0) .top_k0_approx(A, K0) else matrix(0, n, n)
  R <- A - A_k0
  
  # 4. Create Bernoulli subsampling matrix Y
  # Y_ij ~ Bern(1/m) for i < j
  Y <- matrix(0, n, n)
  inds <- which(upper.tri(Y), arr.ind = TRUE)
  draws <- rbinom(nrow(inds), size = 1, prob = 1/m)
  Y[inds] <- draws
  Y <- Y + t(Y) # Symmetrize
  
  # 5. Calculate Test Statistic T_n (Eq. 6)
  # Numerator: sqrt(m) * sum(R_ij * Y_ij)
  # Note: sum(R * Y) is sum(R_ij * Y_ij) over all i,j.
  # Since R_ii = Y_ii = 0, this is sum_{i!=j}
  num <- sqrt(m) * sum(R * Y)
  
  # Denominator: sqrt(2 * sum_{i!=j} R_ij^2)
  # sum(R*R) is sum_{i!=j} R_ij^2
  den <- sqrt(2 * sum(R * R))
  
  T_obs <- 0
  if (den > 1e-10) {
    T_obs <- num / den
  }
  
  # 6. Get p-value from N(0,1) (Theorem 3.1)
  pval <- 2 * (1 - pnorm(abs(T_obs)))
  crit <- qnorm(1 - alpha / 2)
  
  list(
    method = "RIRS-HYF",
    T_obs = T_obs,
    p_value = pval,
    reject = (abs(T_obs) > crit),
    n = n, K0 = K0, m = m
  )
}

#' Convenience wrapper to test H0: K = r0 - 1
rirs_hyf_lambda_r0_nonzero <- function(A, r0, alpha = 0.05, ...) {
  # Pass extra args (...) to rirs_hyf_test_rank
  out <- rirs_hyf_test_rank(A, K0 = r0 - 1, alpha = alpha, ...)
  list(p_value = out$p_value, reject = isTRUE(out$reject))
}