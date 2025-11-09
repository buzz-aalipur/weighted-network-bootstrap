## =======================================================
## my_methods.R (Refactored to match paper algorithms)
## =======================================================
source("helpers.R")

# --- Engine 1: Parametric Bootstrap (for Algs 2, 3, S1, S2, S3) ---
parametric_bootstrap_engine <- function(A, r0, alpha = 0.05, B = 500, S = NULL,
                                        k_neg_mode = "auto", k_pos_mode = r0 - 1L,
                                        method_type = c("test_positive", "test_negative", "estimate_se", "test_rdpg"),
                                        S_scale = "default", normalized = FALSE) {
  method_type <- match.arg(method_type)
  A <- (A + t(A))/2; diag(A) <- 0
  n <- nrow(A)
  
  # Handle S with scaling
  if (is.null(S)) {
    if (S_scale == "small") {
      S <- max(1L, floor(sqrt(n)))
    } else {
      # Use the formula from your paper's recommendation [cite: 3230-3231, 3595-3596]
      S <- max(1L, floor(sqrt(n * max(1, log(n)))))
    }
  }
  
  P_hat <- nbd_smooth(A, S)
  
  # Projector opts based on modes
  projector_opts <- list(k_neg = k_neg_mode, k_pos = k_pos_mode)
  
  if (method_type == "test_negative") {
    projector_opts$k_pos <- "all"  # S1: all positives [cite: 346, 358]
    projector_opts$k_neg <- r0 - 1L
  } else if (method_type == "test_rdpg") {
    projector_opts$k_neg <- 0L  # S2: no negatives [cite: 347-348, 371]
  } else if (method_type == "estimate_se") {
    # For SE, use a standard projection (default opts)
  }
  
  # Call the renamed projector
  P_tilde <- project_Phat(P_hat, r0 = r0, method = "default_projector", opts = projector_opts)
  h_tilde <- matrix_to_graphon(P_tilde)
  
  # Op based on type
  op <- if (method_type == "test_negative") "lambda_neg_r0" else "lambda_pos_r0"
  
  T_obs <- compute_test_stat(A, op = op, r0 = r0, normalized = normalized)
  
  T_boot <- numeric(B)
  for (b in seq_len(B)) {
    Ab <- sample_A_from_P(graphon_P(h_tilde, n, clamp = TRUE)$P, m = 1L, directed = FALSE)
    T_boot[b] <- compute_test_stat(Ab, op = op, r0 = r0, normalized = normalized)
  }
  
  # Outputs based on type
  if (method_type == "estimate_se") {
    se <- sd(T_boot)
    ci <- if (!is.null(alpha)) quantile(T_boot, probs = c(alpha/2, 1 - alpha/2)) else NULL
    return(list(se = se, ci = ci, T_obs = T_obs, T_boot = T_boot))
  } else {
    pval <- mean(T_obs < T_boot) 
    return(list(p_value = pval, reject = (pval < alpha), T_obs = T_obs, T_boot = T_boot))
  }
}

# --- Engine 2: Non-Parametric Bootstrap (for Alg 1) ---
Algorithm_1_NonParametric <- function(A, r0, alpha = 0.05, B = 500, S = NULL, normalized = FALSE) {
  
  A <- (A + t(A))/2; diag(A) <- 0
  n <- nrow(A)
  if (is.null(S)) {
    S <- max(1L, floor(sqrt(n * max(1, log(n)))))
  }
  
  # 1. Get Neighborhoods (list where N[[i]] = vector of neighbor indices for node i)
  # This helper function must be in helpers.R
  N <- get_neighborhoods(A, S)
  
  # 2. Get Observed Statistic (r0 is used as 'i' for SE)
  T_obs <- compute_test_stat(A, op = "lambda_pos_r0", r0 = r0, normalized = normalized)
  
  # 3. Bootstrap Loop
  T_boot <- numeric(B)
  for (b in 1:B) {
    A_boot <- matrix(0L, n, n)
    
    # Sample vertices with replacement [cite: 3203, 3552]
    v_b <- sample(1:n, n, replace = TRUE)
    
    # Reconstruct A_boot based on paper's rule [cite: 3194, 3204, 3537-3557]
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        # Get neighbor sets for the *sampled* vertices
        Ni_indices <- N[[v_b[i]]]
        Nj_indices <- N[[v_b[j]]]
        
        # Create the sampling pool {A_kl} where k in N(v_i), l in N(v_j), k != l
        pool <- c()
        for (k in Ni_indices) {
          pool <- c(pool, A[k, Nj_indices[Nj_indices != k]])
        }
        
        if (length(pool) == 0) {
          A_boot[i, j] <- 0 # Failsafe
        } else {
          # Sample from the pool [cite: 3195, 3204, 3542, 3556]
          A_boot[i, j] <- sample(pool, 1, replace = TRUE)
        }
        A_boot[j, i] <- A_boot[i, j] # Symmetrize
      }
    }
    T_boot[b] <- compute_test_stat(A_boot, op = "lambda_pos_r0", r0 = r0, normalized = normalized)
  }
  
  # 4. Return SE and CI [cite: 3205, 3558]
  se <- sd(T_boot)
  ci <- if (!is.null(alpha)) quantile(T_boot, probs = c(alpha/2, 1 - alpha/2)) else NULL
  return(list(se = se, ci = ci, T_obs = T_obs, T_boot = T_boot))
}


# --- Engine 3: Algorithm 4 (Magnitude Test) ---
Algorithm_4_MagnitudeTest <- function(A, r0, alpha = 0.05, B = 500, S = NULL, normalized = FALSE) {
  # Note: r0 is the *null rank* (R0) [cite: 3292, 3690]
  R0 <- r0 
  r0_stat <- R0 + 1 # We test the (R0+1)-th eigenvalue [cite: 3298]
  
  A <- (A + t(A))/2; diag(A) <- 0
  n <- nrow(A)
  if (is.null(S)) S <- max(1L, floor(sqrt(n * max(1, log(n)))))
  
  P_hat <- nbd_smooth(A, S)
  
  # --- Algorithm 4 Projection Logic ---
  ez <- eigen((P_hat + t(P_hat))/2, symmetric = TRUE)
  vals <- ez$values; vecs <- ez$vectors
  
  # 1. Get R0-th largest eigenvalue by MAGNITUDE [cite: 3294, 3709]
  abs_vals <- abs(vals)
  if (length(abs_vals) < R0 || R0 == 0) {
    threshold <- Inf # Keep nothing
  } else {
    threshold <- sort(abs_vals, decreasing = TRUE)[R0]
  }
  
  # 2. Create mask: Keep all eigenvalues where |val| >= threshold [cite: 3295, 3709]
  mask <- as.numeric(abs_vals >= (threshold - 1e-10))
  Lambda_tilde <- diag(mask * vals)
  P_tilde <- ensure_undirected_noloop(vecs %*% Lambda_tilde %*% t(vecs))
  h_tilde <- matrix_to_graphon(P_tilde)
  # --- End of Projection Logic ---
  
  # --- Test Statistic for Alg 4 (k-th largest magnitude) ---
  compute_stat_alg4 <- function(M, k) {
    vals <- eigen((M + t(M))/2, symmetric = TRUE, only.values = TRUE)$values
    abs_vals <- abs(vals)
    if (length(abs_vals) < k) return(0)
    return(sort(abs_vals, decreasing = TRUE)[k])
  }
  
  T_obs <- compute_stat_alg4(A, r0_stat)
  
  T_boot <- numeric(B)
  for (b in seq_len(B)) {
    Ab <- sample_A_from_P(graphon_P(h_tilde, n, clamp = TRUE)$P, m = 1L)
    T_boot[b] <- compute_stat_alg4(Ab, r0_stat)
  }
  
  pval <- mean(T_obs < T_boot) # Right-tailed test [cite: 3298]
  list(p_value = pval, reject = (pval < alpha), T_obs = T_obs, T_boot = T_boot)
}

## =======================================================
## Dispatcher: Map names to consolidated calls
## =======================================================
run_my_method <- function(which = c("1", "2", "3", "4",
                                    "S1", "S2", "S3",
                                    "SE_parametric", "SE_non_parametric"),
                          A, r0, alpha = 0.05, B = 500, S = NULL, normalized = FALSE) {
  
  which <- as.character(which)
  which <- match.arg(which)
  
  # Pass common parameters
  params <- list(A = A, r0 = r0, alpha = alpha, B = B, S = S, normalized = normalized)
  
  if (which == "1" || which == "SE_non_parametric") {
    # Algorithm 1: Non-parametric SE [cite: 3200-3205]
    return(do.call(Algorithm_1_NonParametric, params))
    
  } else if (which == "4") {
    # Algorithm 4: Magnitude Test [cite: 3294-3298]
    return(do.call(Algorithm_4_MagnitudeTest, params))
    
  } else if (which == "2" || which == "S3" || which == "SE_parametric") {
    # Algorithm 2/S3: Parametric SE [cite: 3209-3214, 382-402]
    params$method_type <- "estimate_se"
    
  } else if (which == "3") {
    # Algorithm 3: Positive Test (keep negs) [cite: 3281-3288]
    params$method_type <- "test_positive"
    params$k_neg_mode <- "auto" # Keep all negatives
    params$k_pos_mode <- r0 - 1 # Keep r0-1 positives
    
  } else if (which == "S1") {
    # Algorithm S1: Negative Test [cite: 352-366]
    params$method_type <- "test_negative"
    
  } else if (which == "S2") {
    # Algorithm S2: RDPG-style Positive Test [cite: 367-380]
    params$method_type <- "test_rdpg"
  }
  
  # Call the main parametric engine with the correct parameters
  do.call(parametric_bootstrap_engine, params)
}