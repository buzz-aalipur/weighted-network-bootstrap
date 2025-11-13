## =========================
## replicate_figure_3.R
## =========================
library(ggplot2)
library(tidyr)
library(pracma)
library(dplyr)
library(scales)
setwd("D:/Courses/U Cincinnati/Research/Bootstrap Based Inference on Eigenvalues of Network_Submitted to Statistics and Computing/Feedbacks/Simulations")

# Load all the functions
source("helpers.R")
source("my_methods.R") # Load the CORRECTED version
source("rirs.R")

# --- Simulation Parameters ---
N_sim <- 200       # n=500 [cite: 2163]
MC_reps <- 5   # MC=300 [cite: 2163]
B_boot <- 5    # B=100 [cite: 2163]
ALPHA <- 0.05    # alpha=0.05 [cite: 2163]
R0 <- 2          # Null rank is 2 [cite: 2163]

# --- Graphon Construction ---
# Define constants from paper [cite: 2158-2159]
random_matrix <- cbind(c(1000,900,800,1000),matrix(rnorm(12), nrow = 4, ncol = 3))
gram = gramSchmidt(random_matrix)
Qnormal = gram$Q # matrix of orthon
u1 <- Qnormal[,1]
u2 <- Qnormal[,2]
u3 <- Qnormal[,3]

mu1 <- 500/230
mu2 <- 130/230

# Helper function to create the 4x4 SBM and the graphon
create_graphon_from_mu <- function(mu) {
  B_mat <- mu1 * (u1 %o% u1) + 
    mu2 * (u2 %o% u2) + 
    mu * (u3 %o% u3)
  
  # The eigenvalues of the integral operator are lambda(B_mat) * (block_size)
  # Here, block_size is 1/4 (for 4 blocks)
  true_eigs <- eigen(B_mat, symmetric = TRUE, only.values = TRUE)$values / 4.0
  
  # Find the 3rd eigenvalue (which is mu/4)
  # It's not guaranteed to be the 3rd, so let's find it
  eigs_sorted_mag <- sort(abs(true_eigs), decreasing = TRUE)
  true_eig3_val <- 0
  if(length(eigs_sorted_mag) >= 3) {
    # Find the eigenvalue corresponding to the 3rd largest magnitude
    cutoff <- eigs_sorted_mag[3]
    idx <- which(abs(true_eigs) >= (cutoff - 1e-10))[3]
    true_eig3_val <- true_eigs[idx]
  }
  
  # Create the step-function (graphon)
  h_func <- function(x, y) {
    # Assign block index (1, 2, 3, or 4)
    i <- floor(x * 4) + 1
    j <- floor(y * 4) + 1
    # Handle x=1 or y=1
    i <- pmin(i, 4)
    j <- pmin(j, 4)
    return(B_mat[i, j])
  }
  
  return(list(graphon = h_func, true_eig_val = true_eig3_val))
}

# --- Set up the simulation sweep ---
# A mu sweep from -0.2 to 0.2 gives an eigenvalue sweep from -0.05 to 0.05
mu_sweep <- seq(-0.2, 0.2, by = 0.04) 
results_list <- list()

cat("Starting simulation to replicate Figure 3...\n")
cat(sprintf("n=%d | MC Reps=%d | Bootstrap B=%d | Points=%d\n", N_sim, MC_reps, B_boot, length(mu_sweep)))

# --- 5. Run the Main Simulation Loop ---
for (mu in mu_sweep) {
  
  graphon_setup <- create_graphon_from_mu(mu)
  graphon_alt <- graphon_setup$graphon
  true_eig_val <- graphon_setup$true_eig_val
  
  cat(sprintf("Running for mu = %.2f (True 3rd eig ≈ %.4f)\n", mu, true_eig_val))
  
  # Initialize rejection counters for each method
  alg_pos_rejects <- 0
  alg_neg_rejects <- 0
  alg4_rejects <- 0
  
  pb <- txtProgressBar(min = 0, max = MC_reps, style = 3, width = 50)
  
  for (r in 1:MC_reps) {
    # 1. Simulate the network
    A <- simulate_graphon(graphon_alt, n = N_sim, m = 1L)$A
    
    # --- Run all 3 tests using the new names ---
    
    # Test 1: Alg 3 (Positive Test, H0: lambda_3+ = 0)
    # This is correctly implemented by "1"
    # res_alg_pos <- run_my_method(which = "3", 
    #                              A = A, 
    #                              r0 = 3, # Test the 3rd positive eigenvalue
    #                              alpha = ALPHA, 
    #                              B = B_boot)
    # if (res_alg_pos$reject) {
    #   alg_pos_rejects <- alg_pos_rejects + 1
    # }
    # 
    #Test 2: Alg S1 (Negative Test, H0: lambda_1- = 0)
    res_alg_neg <- run_my_method(which = "S1",
                                 A = A,
                                 r0 = 1, # Test the 1st negative eigenvalue
                                 alpha = ALPHA,
                                 B = B_boot)
    if (res_alg_neg$reject) {
      alg_neg_rejects <- alg_neg_rejects + 1
    }
    # 
    #Test 3: Alg 4 (Magnitude Test, H0: R+ + R- <= 2)
    
    # R0_NULL = 2
    # res_alg4 <- run_my_method(which = "4",
    #                           A = A,
    #                           r0 = R0_NULL, # Pass the null rank R0=2
    #                           alpha = ALPHA,
    #                           B = B_boot)
    # if (res_alg4$reject) {
    #   alg4_rejects <- alg4_rejects + 1
    # }
    
    setTxtProgressBar(pb, r)
  }
  close(pb)
  
  # Store results
  #power_alg_pos <- alg_pos_rejects / MC_reps
  power_alg_neg <- alg_neg_rejects / MC_reps
  #power_alg4 <- alg4_rejects / MC_reps
  
  # results_list[[paste0(mu, "_pos")]] <- data.frame(
  #   true_eig_val = true_eig_val,
  #   mu = mu,
  #   Method = "Alg 3 (Pos Test, H₀: λ₃⁺=0)",
  #   Power = power_alg_pos
  # )
  results_list[[paste0(mu, "_neg")]] <- data.frame(
    true_eig_val = true_eig_val,
    mu = mu,
    Method = "Alg S1 (Neg Test, H₀: λ₁⁻=0)",
    Power = power_alg_neg
  )
  # results_list[[paste0(mu, "_mag")]] <- data.frame(
  #   true_eig_val = true_eig_val,
  #   mu = mu,
  #   Method = "Alg 4 (Mag Test, H₀: Rank≤2)",
  #   Power = power_alg4
  # )
}

cat("Simulation finished.\n")

# --- 6. Process and Plot Results ---
results_df <- do.call(rbind, results_list)
rownames(results_df) <- NULL

# [cite_start]Create the plot, styled to match Figure 3 [cite: 204, 434-437]
fig3_replication_plot <- ggplot(results_df, aes(x = true_eig_val, y = Power, 
                                                color = Method, linetype = Method, shape = Method)) +
  geom_line(size = 0.8) +
  geom_point(size = 2.5, fill = "white") +
  geom_hline(yintercept = ALPHA, linetype = "dashed", color = "red", size = 0.5) +
  annotate("text", x = min(results_df$true_eig_val), y = ALPHA + 0.04, 
           label = "5% Significance Level", hjust = 0, size = 3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(-0.05, 1.05)) +
  scale_x_continuous(breaks = seq(-0.06, 0.06, by = 0.02)) +
  labs(
    title = "Replication of Figure 3: Power vs. True Eigenvalue",
    subtitle = sprintf("n=%d, MC=%d, B=%d", N_sim, MC_reps, B_boot),
    x = "True Value of 3rd Non-Zero Eigenvalue (λ₃⁺ or λ₁⁻)",
    y = "Rejection Rate (Power)"
  ) +
  scale_color_manual(values = c("Alg 3 (Pos Test, H₀: λ₃⁺=0)" = "green3", 
                                "Alg S1 (Neg Test, H₀: λ₁⁻=0)" = "gray50", 
                                "Alg 4 (Mag Test, H₀: Rank≤2)" = "black")) +
  scale_linetype_manual(values = c("Alg 3 (Pos Test, H₀: λ₃⁺=0)" = "solid", 
                                   "Alg S1 (Neg Test, H₀: λ₁⁻=0)" = "solid", 
                                   "Alg 4 (Mag Test, H₀: Rank≤2)" = "solid")) +
  scale_shape_manual(values = c("Alg 3 (Pos Test, H₀: λ₃⁺=0)" = 22, # Square
                                "Alg S1 (Neg Test, H₀: λ₁⁻=0)" = 21, # Circle
                                "Alg 4 (Mag Test, H₀: Rank≤2)" = 19)) + # Solid Circle
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.title = element_blank())

# Save the plot
ggsave("Figure3_Replication.png", fig3_replication_plot, width = 8, height = 7, dpi = 300)

cat("Plot saved as Figure3_Replication.png\n")
print(fig3_replication_plot)
