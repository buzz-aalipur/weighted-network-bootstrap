## =========================
## hypothesis_test.R
## Main harness to compare methods
## =========================
setwd("D:/Courses/U Cincinnati/Research/Bootstrap Based Inference on Eigenvalues of Network_Submitted to Statistics and Computing/Feedbacks/Simulations")
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
source("helpers.R")
source("my_methods.R")
source("rirs.R")
source("rirs_han_yang_fan.R") # <-- ADDED

## -- choose a graphon by name or pass a custom function
get_graphon <- function(name_or_fun) {
  if (is.function(name_or_fun)) return(name_or_fun)
  if (is.character(name_or_fun) && name_or_fun %in% names(GRAPHONS)) return(GRAPHONS[[name_or_fun]])
  stop("Unknown graphon: ", name_or_fun)
}

## -- run all tests once on the same A
compare_once <- function(A, r0, alpha, my_choice = "1", method_type = c("test_positive", "test_negative", "estimate_se"),
                         B_boot = 500, S = NULL, normalized = FALSE,
                         rirs_calib = "normal", B_rirs = 999) {
  
  method_type <- match.arg(method_type)
  
  if (method_type == "estimate_se") {
    # Run the SE method
    out_my <- run_my_method(which = my_choice, A = A, r0 = r0, alpha = alpha, B = B_boot, S = S, normalized = normalized)
    return(data.frame(method = my_choice, se = out_my$se, ci_low = out_my$ci[1], ci_high = out_my$ci[2]))
    
  } else {
    # Run your method
    out_my   <- run_my_method(which = my_choice, A = A, r0 = r0, alpha = alpha, B = B_boot, S = S, normalized = normalized)
    
    # RIRS (Levin & Levina) only runs for positive rank tests
    if (method_type == "test_positive") {
      out_rirs <- rirs_lambda_r0_nonzero(A, r0 = r0, alpha = alpha, calibration = rirs_calib, B = B_rirs)
      out_rirs_hyf <- rirs_hyf_lambda_r0_nonzero(A, r0 = r0, alpha = alpha)
      
      return(data.frame(
        method = c(my_choice, "RIRS (L&L)", "RIRS (HYF)"),
        p_value = c(out_my$p_value, out_rirs$p_value, out_rirs_hyf$p_value),
        reject  = c(out_my$reject,  out_rirs$reject, out_rirs_hyf$reject)
      ))
      
    } else {
      # For negative tests, only run your method
      return(data.frame(method = my_choice, p_value = out_my$p_value, reject = out_my$reject))
    }
  }
}

compare_power_extended <- function(
    graphon_null, graphon_alt,
    n = 120, r0 = 2, reps = 200, alpha = 0.05,
    my_choice = "1", method_type = "test_positive",
    B_boot = 500, S = NULL, normalized = FALSE,
    rirs_calib = "normal", B_rirs = 999, base_seed = 1
) {
  method_type <- match.arg(method_type)
  set.seed(base_seed)
  
  # Initialize results data frame
  res <- data.frame()
  
  # Helper to run for a graphon
  run_for_graphon <- function(graphon, setting_name) {
    cat(sprintf("Running setting: %s\n", setting_name))
    pb <- txtProgressBar(min = 0, max = reps, style = 3, width = 50)
    
    for (t in seq_len(reps)) {
      A <- simulate_graphon(get_graphon(graphon), n = n, m = 1L, return_A = TRUE, rng_seed = sample.int(1e9,1))$A
      
      # Run your method
      out_my <- run_my_method(which = my_choice, A = A, r0 = r0, alpha = alpha, B = B_boot, S = S, normalized = normalized)
      
      if (method_type == "estimate_se") {
        res_row <- data.frame(setting = setting_name, rep = t,
                              method = my_choice, p_value = NA, reject = NA,
                              se = out_my$se, ci_low = out_my$ci[1], ci_high = out_my$ci[2])
        res <<- rbind(res, res_row)
        
      } else if (method_type == "test_positive") {
        # Run both RIRS competitors
        out_rirs <- rirs_lambda_r0_nonzero(A, r0 = r0, alpha = alpha, calibration = rirs_calib, B = B_rirs)
        out_rirs_hyf <- rirs_hyf_lambda_r0_nonzero(A, r0 = r0, alpha = alpha)
        
        res_rows <- data.frame(
          setting = setting_name, rep = t,
          method = c(my_choice, "RIRS (L&L)", "RIRS (HYF)"),
          p_value = c(out_my$p_value, out_rirs$p_value, out_rirs_hyf$p_value),
          reject  = c(out_my$reject,  out_rirs$reject, out_rirs_hyf$reject),
          se = NA, ci_low = NA, ci_high = NA
        )
        res <<- rbind(res, res_rows)
        
      } else { # e.g., "test_negative"
        res_row <- data.frame(setting = setting_name, rep = t,
                              method = my_choice, p_value = out_my$p_value, reject = out_my$reject,
                              se = NA, ci_low = NA, ci_high = NA)
        res <<- rbind(res, res_row)
      }
      setTxtProgressBar(pb, t)
    }
    close(pb)
  }
  
  ## Null
  run_for_graphon(graphon_null, "Null")
  
  ## Alt
  run_for_graphon(graphon_alt, "Alt")
  
  # Summarize based on type
  summarize <- function(df, is_se = FALSE) {
    if (is_se) {
      df %>% 
        group_by(method) %>%
        summarize(
          n = n, r0 = r0, alpha = alpha, reps = n(),
          mean_se = mean(se, na.rm = TRUE),
          med_ci_width = median(ci_high - ci_low, na.rm = TRUE),
          .groups = "drop"
        )
    } else {
      df %>%
        group_by(method) %>%
        summarize(
          n = n, r0 = r0, alpha = alpha, reps = n(),
          rejection_rate = mean(reject, na.rm = TRUE),
          median_p_value = median(p_value, na.rm = TRUE),
          .groups = "drop"
        )
    }
  }
  
  list(
    params   = list(n=n, r0=r0, alpha=alpha, my_choice=my_choice, method_type=method_type,
                    rirs_calib = rirs_calib, B_boot=B_boot, B_rirs=B_rirs, normalized=normalized),
    size_null = summarize(subset(res, setting=="Null"), is_se = (method_type == "estimate_se")),
    power_alt = summarize(subset(res, setting=="Alt"), is_se = (method_type == "estimate_se")),
    results   = res
  )
}



