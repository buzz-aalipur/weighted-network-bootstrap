## =======================================================
## Replicate Figure 3 from "Weighted_Network_Bootstrap" (V3)
##
## This script tests the power of multiple hypotheses
## at each point along the eigenvalue sweep.
## =======================================================

# --- 1. Load Libraries and Source Files ---
library(pracma)
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

# Load all R files from your working directory
source("helpers.R")
source("my_methods.R") 
source("rirs_han_yang_fan.R") # Load the new RIRS
source("rirs_Levin_Levina.R")
# --- 2. Simulation Parameters ---
N_sim <- 400       # n=500 [cite: 204, 431-432]
MC_reps <- 30     # MC=300 [cite: 204, 431-432]
B_boot <- 30      # B=100 [cite: 204, 431-432]
ALPHA <- 0.05      # alpha=0.05 [cite: 204, 431-432]
R0_NULL <- 2       # The null rank for Algorithm 4 is H0: R+ + R- <= 2 [cite: 204, 431-432]

# --- 3. Graphon Construction ---
# Define constants from paper [cite: 204, 427-428]
set.seed(200)
random_matrix <- cbind(c(1000,900,800,1000),matrix(rnorm(12), nrow = 4, ncol = 3))
gram = gramSchmidt(random_matrix)
Qnormal = gram$Q # matrix of orthon
u1 <- Qnormal[,1]
u2 <- Qnormal[,2]
u3 <- Qnormal[,3]

mu1 <- 500/230
mu2 <- 130/230

create_graphon_from_mu <- function(mu) {
  B_mat <- mu1 * (u1 %o% u1) + 
    mu2 * (u2 %o% u2) + 
    mu * (u3 %o% u3)
  true_eig_val <- mu / 4.0
  
  h_func <- matrix_to_graphon(B_mat)
  return(list(graphon = h_func, true_eig_val = true_eig_val))
}

# --- 4. Set up the simulation sweep ---
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
  
  # Initialize rejection counters
  alg_pos_rejects <- 0
  #alg_neg_rejects <- 0
  #alg4_rejects <- 0
  rirs_hyf_rejects <- 0 # <-- ADDED
  
  pb <- txtProgressBar(min = 0, max = MC_reps, style = 3, width = 50)
  
  for (r in 1:MC_reps) {
    A <- simulate_graphon(graphon_alt, n = N_sim, m = 1L)$A
    
    # Test 1: Alg 3 (Positive Test, H0: lambda_3+ = 0) -> calls "1"
    res_alg_pos <- run_my_method(which = "3", 
                                 A = A, r0 = 3, alpha = ALPHA, B = B_boot)
    if (res_alg_pos$reject) alg_pos_rejects <- alg_pos_rejects + 1
    
    # Test 2: Alg S1 (Negative Test, H0: lambda_1- = 0) -> calls "S1"
    # res_alg_neg <- run_my_method(which = "S1",
    #                              A = A, r0 = 1, alpha = ALPHA, B = B_boot)
    # if (res_alg_neg$reject) alg_neg_rejects <- alg_neg_rejects + 1
    # 
    # Test 3: Alg 4 (Magnitude Test, H0: R+ + R- <= 2) -> calls "4"
    # res_alg4 <- run_my_method(which = "4",
    #                           A = A, r0 = R0_NULL, alpha = ALPHA, B = B_boot)
    # if (res_alg4$reject) alg4_rejects <- alg4_rejects + 1
    # 
    # Test 4: RIRS-HYF (Rank Test, H0: K=2) -> calls rirs_hyf
    res_rirs_hyf <- rirs_hyf_lambda_r0_nonzero(A, r0 = 3, alpha = ALPHA) # r0=3 means K0=2
    if (res_rirs_hyf$reject) rirs_hyf_rejects <- rirs_hyf_rejects + 1
    
    setTxtProgressBar(pb, r)
  }
  close(pb)
  
  # Store results
  results_list[[paste0(mu, "_pos")]] <- data.frame(
    true_eig_val = true_eig_val, Method = "Alg 3 (Pos Test, H₀: λ₃⁺=0)",
    Power = alg_pos_rejects / MC_reps
  )
  # results_list[[paste0(mu, "_neg")]] <- data.frame(
  #   true_eig_val = true_eig_val, Method = "Alg S1 (Neg Test, H₀: λ₁⁻=0)",
  #   Power = alg_neg_rejects / MC_reps
  # )
  # results_list[[paste0(mu, "_mag")]] <- data.frame(
  #   true_eig_val = true_eig_val, Method = "Alg 4 (Mag Test, H₀: Rank≤2)",
  #   Power = alg4_rejects / MC_reps
  # )
  results_list[[paste0(mu, "_hyf")]] <- data.frame(
    true_eig_val = true_eig_val, Method = "RIRS-HYF (Rank Test, H₀: Rank≤2)",
    Power = rirs_hyf_rejects / MC_reps
  )
}

cat("Simulation finished.\n")

# --- 6. Process and Plot Results ---
# --- 6. Process and Plot Results ---
results_df <- do.call(rbind, results_list)
rownames(results_df) <- NULL

# --- THIS IS THE UPDATED PLOTTING CODE ---

# Create the plot, styled to match Figure 3
fig_replication_plot <- ggplot(results_df, aes(x = true_eig_val, y = Power, 
                                               color = Method, linetype = Method, shape = Method)) +
  geom_line(size = 0.8) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = ALPHA, linetype = "dashed", color = "red", size = 0.5) +
  annotate("text", x = min(results_df$true_eig_val), y = ALPHA + 0.04, 
           label = "5% Significance Level", hjust = 0, size = 3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(-0.05, 1.05)) +
  scale_x_continuous(breaks = seq(-0.06, 0.06, by = 0.02)) +
  labs(
    title = "Algorithm 3 versus Rank Test",
    subtitle = sprintf("n=%d, MC=%d, B=%d", N_sim, MC_reps, B_boot),
    x = "True Value of 3rd Non-Zero Eigenvalue (λ₃⁺ or λ₁⁻)",
    y = "Rejection Rate (Power)",
    color = "Hypothesis Test", # Legend title
    shape = "Hypothesis Test", # Legend title
    linetype = "Hypothesis Test" # Legend title
  ) +
  
  # --- ADDED RIRS-HYF to all scales ---
  scale_color_manual(
    values = c("Alg 3 (Pos Test, H₀: λ₃⁺=0)" = "green3", 
               "Alg S1 (Neg Test, H₀: λ₁⁻=0)" = "gray50", 
               "Alg 4 (Mag Test, H₀: Rank≤2)" = "black",
               "RIRS-HYF (Rank Test, H₀: Rank≤2)" = "blue"), # <-- ADDED
    labels = c("Alg 3 (Pos Test, H₀: λ₃⁺=0)" = "Algorithm 3 (Positive Test)", 
               "Alg S1 (Neg Test, H₀: λ₁⁻=0)" = "Algorithm S1 (Negative Test)", 
               "Alg 4 (Mag Test, H₀: Rank≤2)" = "Algorithm 4 (Magnitude Test)",
               "RIRS-HYF (Rank Test, H₀: Rank≤2)" = "RIRS (Han, Yang, Fan)") # <-- ADDED
  ) +
  scale_linetype_manual(
    values = c("Alg 3 (Pos Test, H₀: λ₃⁺=0)" = "solid", 
               "Alg S1 (Neg Test, H₀: λ₁⁻=0)" = "solid", 
               "Alg 4 (Mag Test, H₀: Rank≤2)" = "solid",
               "RIRS-HYF (Rank Test, H₀: Rank≤2)" = "dashed"), # <-- ADDED
    labels = c("Alg 3 (Pos Test, H₀: λ₃⁺=0)" = "Algorithm 3 (Positive Test)", 
               "Alg S1 (Neg Test, H₀: λ₁⁻=0)" = "Algorithm S1 (Negative Test)", 
               "Alg 4 (Mag Test, H₀: Rank≤2)" = "Algorithm 4 (Magnitude Test)",
               "RIRS-HYF (Rank Test, H₀: Rank≤2)" = "RIRS (Han, Yang, Fan)") # <-- ADDED
  ) +
  scale_shape_manual(
    values = c("Alg 3 (Pos Test, H₀: λ₃⁺=0)" = 22, # Square
               "Alg S1 (Neg Test, H₀: λ₁⁻=0)" = 21, # Circle
               "Alg 4 (Mag Test, H₀: Rank≤2)" = 19, # Solid Circle
               "RIRS-HYF (Rank Test, H₀: Rank≤2)" = 17), # Triangle
    labels = c("Alg 3 (Pos Test, H₀: λ₃⁺=0)" = "Algorithm 3 (Positive Test)", 
               "Alg S1 (Neg Test, H₀: λ₁⁻=0)" = "Algorithm S1 (Negative Test)", 
               "Alg 4 (Mag Test, H₀: Rank≤2)" = "Algorithm 4 (Magnitude Test)",
               "RIRS-HYF (Rank Test, H₀: Rank≤2)" = "RIRS (Han, Yang, Fan)") # <-- ADDED
  ) +
  
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold")
  )

# Save the plot
ggsave("Figure3_Replication_with_HYF.png", fig_replication_plot, width = 8, height = 7, dpi = 300)

cat("Plot saved as Figure3_Replication_with_HYF.png\n")
print(fig_replication_plot)