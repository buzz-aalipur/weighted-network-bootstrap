## =========================
## rirs.R
## Paper-faithful RIRS test
## =========================
source("helpers.R")

.as_sym_no_diag <- function(A) {
  A <- as.matrix(A)
  if (!is.numeric(A) || nrow(A) != ncol(A)) stop("A must be numeric square.")
  A[!is.finite(A)] <- 0
  A <- 0.5 * (A + t(A)); diag(A) <- 0; A
}

.top_k0_approx <- function(A, K0) {
  if (K0 <= 0) return(matrix(0, nrow(A), ncol(A)))
  ev <- eigen(A, symmetric = TRUE, only.values = FALSE)
  ord <- order(abs(ev$values), decreasing = TRUE)
  vals <- ev$values[ord]; vecs <- ev$vectors[, ord, drop = FALSE]
  K0 <- min(K0, length(vals))
  Uk <- vecs[, seq_len(K0), drop = FALSE]
  Dk <- diag(vals[seq_len(K0)], nrow = K0, ncol = K0)
  Uk %*% Dk %*% t(Uk)
}

.sample_offdiag_pairs <- function(n, m) {
  M <- n * (n - 1) / 2
  if (m > M) m <- M
  idx_i <- integer(M); idx_j <- integer(M); cnt <- 0L
  for (i in 1:(n-1)) {
    len <- n - i; rng <- (cnt + 1L):(cnt + len)
    idx_i[rng] <- i; idx_j[rng] <- (i+1):n; cnt <- cnt + len
  }
  sel <- sample.int(M, m, replace = FALSE)
  cbind(idx_i[sel], idx_j[sel])
}

.rirs_T_once <- function(R, pairs_ij) {
  x <- R[pairs_ij]; s <- sum(x); v <- sum(x * x)
  if (v <= 0) return(0)
  sqrt(nrow(pairs_ij)) * s / sqrt(2 * v)
}

rirs_test_rank <- function(A, K0, alpha = 0.05,
                           m = NULL,
                           calibration = c("normal", "wild"),
                           B = 999,
                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  calibration <- match.arg(calibration)
  A <- .as_sym_no_diag(A)
  n <- nrow(A)
  
  M_tot <- n * (n - 1) / 2
  if (is.null(m)) m <- min(M_tot, max(200L, as.integer(n * max(1, log(n))))) else m <- max(50L, min(as.integer(m), M_tot))
  
  A_k0 <- if (K0 > 0) .top_k0_approx(A, K0) else matrix(0, n, n)
  R <- A - A_k0; diag(R) <- 0
  
  pairs_ij <- .sample_offdiag_pairs(n, m)
  T_obs <- .rirs_T_once(R, pairs_ij)
  
  if (calibration == "normal") {
    pval <- 2 * (1 - pnorm(abs(T_obs)))
    crit <- qnorm(1 - alpha / 2)
    return(list(method = "RIRS (normal)", T_obs = T_obs, p_value = pval, reject = abs(T_obs) > crit,
                n = n, K0 = K0, m = m))
  }
  
  x <- R[pairs_ij]; v <- sum(x * x)
  if (v <= 0) return(list(method = "RIRS (wild)", T_obs = T_obs, p_value = 1, reject = FALSE, n = n, K0 = K0, m = m, B = B))
  Tb <- numeric(B)
  for (b in seq_len(B)) {
    eps <- sample(c(-1,1), length(x), replace = TRUE)
    s_b <- sum(eps * x)
    Tb[b] <- sqrt(nrow(pairs_ij)) * s_b / sqrt(2 * v)
  }
  pval <- (sum(abs(Tb) >= abs(T_obs)) + 1) / (B + 1)
  crit <- unname(quantile(abs(Tb), 1 - alpha/2, type = 8))
  list(method = "RIRS (wild)", T_obs = T_obs, p_value = pval, reject = abs(T_obs) > crit,
       n = n, K0 = K0, m = m, B = B)
}

## Convenience: r0-th spike nonzero? (H0: K = r0-1)
rirs_lambda_r0_nonzero <- function(A, r0, alpha = 0.05, calibration = "normal", B = 999) {
  out <- rirs_test_rank(A, K0 = r0 - 1, alpha = alpha, calibration = calibration, B = B)
  list(p_value = out$p_value, reject = isTRUE(out$reject))
}

