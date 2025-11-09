## =========================
## helpers.R (Final Version)
## Shared utilities + graphons
## =========================

## ---- Nice-to-haves ----
`%||%` <- function(a, b) if (is.null(a)) b else a
safelog <- function(...) cat(sprintf(...), "\n")

.strict_scalar <- function(val) {
  if (length(val) == 0) return(0)
  v <- suppressWarnings(as.numeric(val))
  v <- v[is.finite(v)]
  if (length(v) == 0) return(0)
  v[1]
}

ensure_undirected_noloop <- function(M, clip = TRUE) {
  M <- (M + t(M)) / 2
  diag(M) <- 0
  if (clip) {
    M[M < 0] <- 0
    M[M > 1] <- 1
  }
  M
}

## ---- Graphon helpers ----
graphon_P <- function(graphon, n, xi = NULL, clamp = TRUE, symmetric = TRUE) {
  if (is.matrix(graphon)) {
    P <- graphon
    stopifnot(nrow(P) == ncol(P))
    return(list(P = P, xi = NULL))
  }
  stopifnot(is.function(graphon))
  if (is.null(xi)) xi <- runif(n)
  
  P <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in i:n) {
      val <- .strict_scalar(graphon(xi[i], xi[j]))
      P[i, j] <- val
      if (i != j) P[j, i] <- val
    }
  }
  
  if (clamp) P <- pmin(pmax(P, 0), 1)
  list(P = P, xi = xi)
}
sample_A_from_P <- function(P, m = 1L, directed = FALSE, self_loops = FALSE) {
  stopifnot(is.matrix(P), nrow(P) == ncol(P), m >= 0, is.finite(m))
  n <- nrow(P)
  if (directed) {
    A <- matrix(rbinom(n * n, size = m, prob = as.vector(P)), n, n)
    if (!self_loops) diag(A) <- 0
    storage.mode(A) <- "integer"
    return(A)
  } else {
    A <- matrix(0L, n, n)
    inds <- which(upper.tri(A, diag = self_loops), arr.ind = TRUE)
    draws <- rbinom(nrow(inds), size = m, prob = P[inds])
    A[inds] <- draws
    A <- A + t(A)
    if (!self_loops) diag(A) <- 0L
    storage.mode(A) <- "integer"
    return(A)
  }
}

simulate_graphon <- function(graphon, n, m = 1L, return_A = TRUE,
                             xi = NULL, clip = TRUE, rng_seed = NULL,
                             directed = FALSE, self_loops = FALSE) {
  stopifnot(is.function(graphon), n >= 1, m >= 1)
  if (!is.null(rng_seed)) set.seed(rng_seed)
  gp <- graphon_P(graphon = graphon, n = n, xi = xi, clamp = clip, symmetric = !directed)
  P  <- gp$P; xi <- gp$xi
  if (!self_loops) diag(P) <- 0
  if (!return_A) return(list(P = P, A = NULL, xi = xi))
  A <- sample_A_from_P(P, m = m, directed = directed, self_loops = self_loops)
  list(P = P, A = A, xi = xi)
}

matrix_to_graphon <- function(M) {
  n <- nrow(M)
  function(x, y) {
    i <- pmin(pmax(1, findInterval(x, seq(0, 1, length.out = n + 1))), n)
    j <- pmin(pmax(1, findInterval(y, seq(0, 1, length.out = n + 1))), n)
    p <- M[i, j]
    if (!is.finite(p)) p <- 0
    if (p < 0) p <- 0
    if (p > 1) p <- 1
    p
  }
}

## ---- Eigen/stat helpers ----
eig_sym <- function(M, k = NULL) {
  stopifnot(is.matrix(M), nrow(M) == ncol(M))
  M <- (M + t(M)) / 2
  n <- nrow(M)
  if (is.null(k) || k >= n - 1L) {
    ez <- eigen(M, symmetric = TRUE)
    return(list(values = ez$values, vectors = ez$vectors))
  } else {
    ez <- RSpectra::eigs_sym(M, k, which = "LA")
    return(list(values = as.numeric(ez$values), vectors = as.matrix(ez$vectors)))
  }
}

kth_positive_eigen <- function(vals, k) {
  pos <- sort(vals[vals > 0], decreasing = TRUE)
  if (length(pos) < k) return(0)
  pos[k]
}

kth_negative_eigen <- function(vals, k) {
  neg <- vals[vals < -1e-10] 
  if (length(neg) < k) return(0)
  sort(neg, decreasing = FALSE)[k] 
}

# CORRECTED compute_test_stat function
compute_test_stat <- function(A, op = c("lambda_pos_r0", "lambda_neg_r0"), r0 = 1, fast = TRUE, normalized = FALSE) {
  op <- match.arg(op)
  M  <- (A + t(A)) / 2
  n  <- nrow(M)
  
  if (fast) {
    k_target <- max(2L, r0 + 2L)
    k <- min(k_target, max(1L, n - 2L))
    
    if (op == "lambda_pos_r0") {
      rspectra_which <- "LA" # Largest Algebraic (most positive)
    } else {
      rspectra_which <- "SA" # Smallest Algebraic (most negative)
    }
    
    vals <- tryCatch(
      RSpectra::eigs_sym(M, k = k, which = rspectra_which)$values,
      error = function(e) {
        eigen(M, symmetric = TRUE, only.values = TRUE)$values
      }
    )
    
  } else {
    vals <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  }
  
  if (op == "lambda_pos_r0") {
    eig_val <- kth_positive_eigen(vals, r0)
  } else if (op == "lambda_neg_r0") {
    eig_val <- kth_negative_eigen(vals, r0)
  } else {
    stop("Unknown test op: ", op)
  }
  
  if (normalized) eig_val <- eig_val / n
  eig_val
}

## ---- Neighborhood smoothing ----
Distance <- function(A) {
  n <- nrow(A); D <- matrix(0, n, n)
  A2 <- (A %*% A) / n
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      v <- abs(A2[i, ] - A2[j, ])
      v[i] <- 0; v[j] <- 0
      D[i, j] <- D[j, i] <- max(v)
    }
  }
  D
}

ssmallest_mask <- function(v, S) {
  idx <- order(v, decreasing = FALSE)[1:S]
  m <- numeric(length(v)); m[idx] <- 1; m
}

# New helper to get the actual neighbor indices as a list
# N[[i]] will be a vector of neighbor indices for node i
get_neighborhoods <- function(A, S) {
  n <- nrow(A)
  D <- Distance(A)
  diag(D) <- 0 # Set diag to 0 to include self in N(i) [cite: 3197, 3544]
  
  # apply over rows (margin=1)
  neighbor_indices <- apply(D, 1, function(row_dists) {
    # Get indices of the S smallest distances
    order(row_dists, decreasing = FALSE)[1:S]
  })
  # Return as a list of vectors
  as.list(as.data.frame(neighbor_indices))
}


nbd_smooth <- function(A, S) {
  n <- nrow(A)
  D <- Distance(A)
  diag(D) <- max(D) + 1 # Exclude self from smoothing average
  K <- t(apply(D, 1, ssmallest_mask, S = S))
  diag(K) <- 0
  
  row_sums <- rowSums(K)
  zero_sum_rows <- which(row_sums == 0)
  if (length(zero_sum_rows) > 0) {
    # Failsafe: if a node has no other neighbors (e.g. S=1), 
    # its smoothed row is just its own row.
    K[zero_sum_rows, zero_sum_rows] <- 1 
    row_sums[zero_sum_rows] <- 1
  }
  
  K <- K / row_sums
  P <- K %*% A
  (P + t(P)) / 2
}

## ---- Projector registry (Main projector) ----

# Renamed from .project_alg3
.project_algs_123S1S2 <- function(P_hat, r0, opts = list()) {
  k_neg_mode <- opts$k_neg %||% "all"
  tol_pos    <- opts$tol_pos %||% 1e-10
  k_pos_over <- opts$k_pos_override
  k_pos_mode <- opts$k_pos %||% (r0 - 1L)  # Default to r0-1 for pos test
  
  P_e <- (P_hat + t(P_hat)) / 2
  n   <- nrow(P_e)
  
  ev <- eigen(P_e, symmetric = TRUE, only.values = TRUE)$values
  
  # Get the total number of positive eigenvalues
  total_positive_eigs <- sum(ev > tol_pos)
  
  if (identical(k_pos_mode, "all")) {
    k_pos <- total_positive_eigs
  } else if (!is.null(k_pos_over)) {
    k_pos <- as.integer(k_pos_over)
  } else {
    k_pos <- max(0L, as.integer(k_pos_mode)) # This is k_pos = r0 - 1
    
    # Your rule: if not enough positive eigs, set all to 0
    if (total_positive_eigs < k_pos) {
      k_pos <- 0L 
    }
  }
  k_pos <- max(0L, min(k_pos, n - 2L))
  
  
  if (identical(k_neg_mode, "all")) {
    k_neg <- min(n - 2L, opts$k_neg_cap %||% (n - 2L))
  } else if (identical(k_neg_mode, "auto")) {
    k_neg <- sum(ev < -tol_pos) # Keep all negative eigs
    k_neg <- max(0L, min(k_neg, n - 2L))
  } else {
    k_neg <- max(0L, min(as.integer(k_neg_mode), n - 2L))
  }
  
  Initial <- matrix(0, n, n)
  
  if (k_pos > 0) {
    z_pos <- tryCatch(RSpectra::eigs_sym(P_e, k = k_pos, which = "LA"),
                      error = function(e) eigen(P_e, symmetric = TRUE))
    if (!is.null(z_pos$vectors)) {
      pos_val <- as.numeric(z_pos$values); pos_vec <- as.matrix(z_pos$vectors)
    } else {
      ord <- order(z_pos$values, decreasing = TRUE)[seq_len(k_pos)]
      pos_val <- z_pos$values[ord]; pos_vec <- z_pos$vectors[, ord, drop = FALSE]
    }
    Initial <- Initial + pos_vec %*% diag(pos_val, length(pos_val)) %*% t(pos_vec)
  }
  
  if (k_neg > 0) {
    # Use "SA" (Smallest Algebraic) to get largest-magnitude negative eigs
    z_neg <- tryCatch(RSpectra::eigs_sym(P_e, k = k_neg, which = "SA"),
                      error = function(e) eigen(P_e, symmetric = TRUE))
    if (!is.null(z_neg$vectors)) {
      neg_val <- as.numeric(z_neg$values); neg_vec <- as.matrix(z_neg$vectors)
    } else {
      ord <- order(z_neg$values, decreasing = FALSE)[seq_len(k_neg)]
      neg_val <- z_neg$values[ord]; neg_vec <- z_neg$vectors[, ord, drop = FALSE]
    }
    Initial <- Initial + neg_vec %*% diag(neg_val, length(neg_val)) %*% t(neg_vec)
  }
  
  (Initial + t(Initial)) / 2
}

PROJECTORS <- new.env(parent = emptyenv())
# Register the function with its new name
PROJECTORS$map <- list(default_projector = .project_algs_123S1S2)

register_projector <- function(name, fun) {
  stopifnot(is.character(name), length(name) == 1, is.function(fun))
  PROJECTORS$map[[name]] <- fun
  invisible(TRUE)
}

project_Phat <- function(P_hat, r0, method = "default_projector", opts = list()) {
  # This makes 'method' argument optional, it will always use the default
  f <- PROJECTORS$map[[method[1]]]
  if (!is.function(f)) stop("No projector registered for '", method[1], "'.")
  f(P_hat, r0, opts)
}


## ---- Graphons (your list) ----
G1  <- function(x,y){ K <- floor(log(1e6)); g <- seq(0,1,by=(1/K)); xp <- sum(g < x); yp <- sum(g < y); if (xp==yp) xp/(K+1) else 0.3/(K+1) }
G2  <- function(x,y){ (sin(5*pi*(x+y-1)+1)/2)+0.5 }
G3  <- function(x,y){ 1-(1 + exp(15*(0.8*abs(x-y))^(4/5)-0.1))^(-1) }
G4  <- function(x,y){ ((x^2 + y^2)/3)*cos(1/(x^2 + y^2)) + 0.15 }
G5  <- function(x,y){ x*y }
G6  <- function(x,y){ exp(-(x^(0.7)+ y^(0.7))) }
G7  <- function(x,y){ (1/4)*(x^2+y^2+x^(1/2)+y^(1/2)) }
G8  <- function(x,y){ (1/2)*(x+y) }
G9  <- function(x,y){ 1/(1+exp(-10*(x^2+y^2))) }
G10 <- function(x,y){ abs(x-y) }
G11 <- function(x,y){ 1/(1+(exp(-(max(x,y)^2+min(x,y)^4)))) }
G12 <- function(x,y){ exp(-max(x,y)^(3/4)) }
G13 <- function(x,y){ exp(-(1/2)*(min(x,y)+x^(1/2)+y^(1/2))) }
G14 <- function(x,y){ log(1 + 0.5*max(x, y)) }

## Optional: choose by name
GRAPHONS <- list(G1=G1,G2=G2,G3=G3,G4=G4,G5=G5,G6=G6,G7=G7,G8=G8,G9=G9,G10=G10,G11=G11,G12=G12,G13=G13,G14=G14)