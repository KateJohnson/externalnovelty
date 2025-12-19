# ==============================================================================
# LOAD EPIC (NOT A TEST)
# ==============================================================================

library(epicUS)
library(tidyverse)
library(dplyr)
library(copula)

settings <- get_default_settings()
settings$record_mode <- 2
settings$n_base_agents <- 1.5e4
init_session(settings = settings)
input <- get_input()
time_horizon <- 25
input$values$global_parameters$time_horizon <- time_horizon
run(input = input$values)
op <- Cget_output()
all_events2 <- Cget_output_ex()
all_events2 <- as.data.frame(Cget_all_events_matrix())

terminate_session()

# ==============================================================================
#  EPIC POPULATION TEST (COHORT CREATION)
# ==============================================================================
# Test objective: This section confirms only patients with a confirmed diagnosis 
# of COPD are captured and simulated patients have at least 1 year of follow up 
# post-diagnosis to be eligible for matching

# Create the calendar year variable to capture temporal trends
all_events2 <- all_events2 %>%
  mutate(calendar_year = floor(2015 + time_at_creation + local_time)) # CORRECTION: add floor to prevent rounding up partial years
         
# Identify patients that fulfill eligibility criteria
cand <- all_events2 %>%
  # Filter to those that are diagnosed
  filter(diagnosis > 0, gold > 0, FEV1_baseline > 0) %>%
  
  # Select ONLY the Index Date row (The first confirmed diagnosis event)
  group_by(id) %>%
  slice_min(local_time, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  
  # Define 'index_time' 
  mutate(index_time = calendar_year) %>%
  
  # Apply Inclusion Criteria
  filter(
    # Index Date must be in 2016 or 2017
    index_time >= 2016 & index_time < 2018,
    
    # Look back 1 year for history
    # If local_time >= 1, it guarantees >1 year of history exists prior to this row.
    local_time >= 1 
  )

# The following is a test to confirm the above is correctly creating the EPIC cohort
test_cand <- function(cand_df) {
  
  message("--- Sanity Check ---")
  
  # 1. CHECK TOTALS
  total_patients <- nrow(cand_df)
  message(paste("Total Candidates Selected:", total_patients))
  
  # 2. CHECK DATES (Should be 2016 and 2017)
  bad_dates <- cand_df %>% filter(index_time < 2016 | index_time >= 2018)
  
  if(nrow(bad_dates) == 0) {
    message("✅ DATE CHECK: PASS (All index dates are within 2016-2017)")
  } else {
    message(paste("❌ DATE CHECK: FAIL (Found", nrow(bad_dates), "rows outside the window)"))
    print(head(bad_dates %>% select(id, index_time)))
  }
  
  # 3. CHECK HISTORY (Should be local_time >= 1)
  bad_history <- cand_df %>% filter(local_time < 1)
  
  if(nrow(bad_history) == 0) {
    message("✅ HISTORY CHECK: PASS (All patients have >1 year history)")
  } else {
    message(paste("❌ HISTORY CHECK: FAIL (Found", nrow(bad_history), "rows with insufficient history)"))
    print(head(bad_history %>% select(id, local_time)))
  }
  
  # 4. VISUAL SAMPLE (Look at 5 random rows yourself)
  message("\n--- Visual Sample (Do these look right?) ---")
  sample_n(cand_df, min(5, nrow(cand_df))) %>% 
    select(id, index_time, local_time, gold, diagnosis) %>%
    print()
}

# Run the check
test_cand(cand)


# ==============================================================================
#  EPIC COHORT SETUP FOR COPULA (NOT A TEST)
# ==============================================================================
# This section details data preparation for copula (this is a data prep step and NOT a test)

# If patient has under 0.1 pack years and is currently a non-smoker they are considered non-smoker
PY_EVER_THRESH <- 0.1    

cand <- cand %>%
  mutate(
    smoking    = as.integer(as.character(smoking_status)),
    pack_years = as.numeric(pack_years),
    
    gold       = as.integer(gold),
    exac_mod   = as.integer(exac_history_n_moderate),
    exac_sev   = as.integer(exac_history_n_severe_plus)
  ) %>%
  mutate(
    smoking_status = dplyr::case_when(
      smoking == 1L ~ 1L,  # Current
      
      # "Never Smoker" logic: Must be non-current AND pack_years below threshold.
      pack_years <= PY_EVER_THRESH ~ 0L,
      
      # Everyone else is Former
      TRUE ~ 2L
    )
  )

# Load NOVELTY RDS data (can be found on the github repo)
target <- readRDS("novelty_target_COPD.rds") #CORRECTION: add readRDS and File name wrong

# Load Parameters from target 
lev_smoking    <- target$smoking_levels
prob_smoking   <- target$smoking_probs 
lev_gold       <- target$GOLD_levels
prob_gold      <- target$GOLD_probs     
p_sex          <- target$p_sex_baseline 
p_mmrc         <- target$p_mMRC_baseline 
ml_age         <- target$meanlog_ageCOPD_baseline
sl_age         <- target$sdlog_ageCOPD_baseline
lambda_exac_mod <- target$lambda_exacerbations_moderate_baseline
lambda_exac_sev <- target$lambda_exacerbations_severeplus_baseline

# ==============================================================================
#  HELPER FUNCTION SET UP FOR COPULA (NOT A TEST)
# ==============================================================================

clip01 <- function(U, eps = 1e-9) {pmin(pmax(U, eps), 1 - eps)}

pit_categorical_from_probs <- function(values, level_order, level_probs) { # this function assigns numeric values
  #to the smoking status with values uniformly distributed between the prop of smoking categories observed in NOVELTY
  idx <- match(as.integer(values), as.integer(level_order))
  if (anyNA(idx)) stop("Found NAs or values not in level_order.")
  p <- as.numeric(level_probs); p <- p / sum(p)
  cum_all <- cumsum(p); lower_all <- c(0, head(cum_all, -1))
  lower <- lower_all[idx]; upper <- cum_all[idx]
  clip01(stats::runif(length(values), min = lower, max = upper))
}

pit_bernoulli <- function(y01, p) { # this function assigns numeric values to sex the values uniformly distributed between 0 and 1-p_sex when sex in EPIC
  # = female, and 1-p_sex and 1 when sex = male.
  out <- rep(NA_real_, length(y01))
  i0 <- y01 == 0; i1 <- y01 == 1
  if(sum(i0, na.rm=T)>0) out[i0] <- runif(sum(i0, na.rm=T), 0, 1-p)
  if(sum(i1, na.rm=T)>0) out[i1] <- runif(sum(i1, na.rm=T), 1-p, 1)
  clip01(out)
}


pit_poisson <- function(k, lambda) {
  if (any(lambda <= 0 | !is.finite(lambda))) stop("Lambda must be positive.")
  upper <- ppois(k, lambda = lambda)
  lower <- ppois(k - 1, lambda = lambda) # Returns 0 if k=0
  u <- runif(length(k), min = lower, max = upper)
  clip01(u)
}

# ==============================================================================
#  COPULA HELPER FUNCTION TEST
# ==============================================================================
# This section transforms the EPIC characteristics into uniform probabilities (U) 
# based on NOVELTY distributions. These U values are the required inputs for the 
# Copula model, which generates the likelihood scores used for weighting.
#
# Test objective: This section confirms that the helper functions for the copula 
# are working as intended; sanity check of calculated versus expected values 
# for the matching characteristics.
# Matching characteristics include: Gender, Age, mMRC, prior moderate exacerbations, 
# prior severe exacerbations, GOLD stage, and smoking status.


message("\n---------------------------------------------------")
message("TEST: SMOKING STATUS")
message("---------------------------------------------------")
smk_cuts <- cumsum(prob_smoking) 
u_smoking <- pit_categorical_from_probs(cand$smoking_status, lev_smoking, prob_smoking)

message("Target Probs: ", paste(round(prob_smoking, 3), collapse=", "))

# Status 0 (Never) -> Expect 0 to cut[1]
u_sub <- u_smoking[cand$smoking_status == 0]
pass  <- all(u_sub >= 0 & u_sub <= smk_cuts[1] + 1e-9)
tag   <- if(pass) "[PASS]" else "[FAIL]"
message(sprintf("%s Status 0 (Never):   %s (Expect 0 - %.4f)", 
                tag, paste(round(range(u_sub), 4), collapse=" - "), smk_cuts[1]))

# Status 1 (Current) -> Expect cut[1] to cut[2]
u_sub <- u_smoking[cand$smoking_status == 1]
pass  <- all(u_sub >= smk_cuts[1] - 1e-9 & u_sub <= smk_cuts[2] + 1e-9)
tag   <- if(pass) "[PASS]" else "[FAIL]"
message(sprintf("%s Status 1 (Current): %s (Expect %.4f - %.4f)", 
                tag, paste(round(range(u_sub), 4), collapse=" - "), smk_cuts[1], smk_cuts[2]))

# Status 2 (Former) -> Expect cut[2] to 1
u_sub <- u_smoking[cand$smoking_status == 2]
pass  <- all(u_sub >= smk_cuts[2] - 1e-9 & u_sub <= 1 + 1e-9)
tag   <- if(pass) "[PASS]" else "[FAIL]"
message(sprintf("%s Status 2 (Former):  %s (Expect %.4f - 1)", 
                tag, paste(round(range(u_sub), 4), collapse=" - "), smk_cuts[2]))


message("\n---------------------------------------------------")
message("TEST: GOLD STAGE")
message("---------------------------------------------------")
gold_cuts <- cumsum(prob_gold) 
u_gold    <- pit_categorical_from_probs(cand$gold, lev_gold, prob_gold)

message("Target Probs: ", paste(round(prob_gold, 3), collapse=", "))

# GOLD 1
u_sub <- u_gold[cand$gold == 1]
pass  <- all(u_sub >= 0 & u_sub <= gold_cuts[1] + 1e-9)
tag   <- if(pass) "[PASS]" else "[FAIL]"
message(sprintf("%s GOLD 1: %s (Expect 0 - %.4f)", 
                tag, paste(round(range(u_sub), 4), collapse=" - "), gold_cuts[1]))

# GOLD 2
u_sub <- u_gold[cand$gold == 2]
pass  <- all(u_sub >= gold_cuts[1] - 1e-9 & u_sub <= gold_cuts[2] + 1e-9)
tag   <- if(pass) "[PASS]" else "[FAIL]"
message(sprintf("%s GOLD 2: %s (Expect %.4f - %.4f)", 
                tag, paste(round(range(u_sub), 4), collapse=" - "), gold_cuts[1], gold_cuts[2]))

# GOLD 3
u_sub <- u_gold[cand$gold == 3]
pass  <- all(u_sub >= gold_cuts[2] - 1e-9 & u_sub <= gold_cuts[3] + 1e-9)
tag   <- if(pass) "[PASS]" else "[FAIL]"
message(sprintf("%s GOLD 3: %s (Expect %.4f - %.4f)", 
                tag, paste(round(range(u_sub), 4), collapse=" - "), gold_cuts[2], gold_cuts[3]))

# GOLD 4
u_sub <- u_gold[cand$gold == 4]
pass  <- all(u_sub >= gold_cuts[3] - 1e-9 & u_sub <= 1 + 1e-9)
tag   <- if(pass) "[PASS]" else "[FAIL]"
message(sprintf("%s GOLD 4: %s (Expect %.4f - 1)", 
                tag, paste(round(range(u_sub), 4), collapse=" - "), gold_cuts[3]))


message("\n---------------------------------------------------")
message("TEST: SEX")
message("---------------------------------------------------")
u_sex <- pit_bernoulli(cand$female, p = p_sex)
sex_cut <- 1 - p_sex

message("Target Prob Female: ", p_sex, " | Cutoff: ", round(sex_cut, 4))

# Males (0) -> 0 to Cutoff
u_sub <- u_sex[cand$female == 0]
pass  <- all(u_sub >= 0 & u_sub <= sex_cut + 1e-9)
tag   <- if(pass) "[PASS]" else "[FAIL]"
message(sprintf("%s Males (0):   %s (Expect 0 - %.4f)", 
                tag, paste(round(range(u_sub), 4), collapse=" - "), sex_cut))

# Females (1) -> Cutoff to 1
u_sub <- u_sex[cand$female == 1]
pass  <- all(u_sub >= sex_cut - 1e-9 & u_sub <= 1 + 1e-9)
tag   <- if(pass) "[PASS]" else "[FAIL]"
message(sprintf("%s Females (1): %s (Expect %.4f - 1)", 
                tag, paste(round(range(u_sub), 4), collapse=" - "), sex_cut))


message("\n---------------------------------------------------")
message("TEST: mMRC")
message("---------------------------------------------------")

if (length(p_mmrc) > 1) {
  # Categorical Case
  mmrc_cuts <- cumsum(p_mmrc)
  u_mmrc <- pit_categorical_from_probs(cand$dyspnea, names(p_mmrc), p_mmrc)
  message("Target Probs (Cat): ", paste(round(p_mmrc, 3), collapse=", "))
  
  # Simple full range check for categorical (granular checks similar to GOLD possible)
  rng <- range(u_mmrc)
  pass <- rng[1] > 0 & rng[2] < 1
  tag  <- if(pass) "[PASS]" else "[FAIL]"
  message(sprintf("%s mMRC U-Range: %s", tag, paste(round(rng, 4), collapse=" - ")))
  
} else {
  # Binary Case
  u_mmrc <- pit_bernoulli(cand$dyspnea, p = p_mmrc)
  mmrc_cut <- 1 - p_mmrc
  message("Target Prob (Bin): ", p_mmrc, " | Cutoff: ", round(mmrc_cut, 4))
  
  # Score 0 -> 0 to Cutoff
  u_sub <- u_mmrc[cand$dyspnea == 0]
  pass  <- all(u_sub >= 0 & u_sub <= mmrc_cut + 1e-9)
  tag   <- if(pass) "[PASS]" else "[FAIL]"
  message(sprintf("%s Score 0: %s (Expect 0 - %.4f)", 
                  tag, paste(round(range(u_sub), 4), collapse=" - "), mmrc_cut))
  
  # Score 1 (if exists)
  if(any(cand$dyspnea == 1)) {
    u_sub <- u_mmrc[cand$dyspnea == 1]
    pass  <- all(u_sub >= mmrc_cut - 1e-9 & u_sub <= 1 + 1e-9)
    tag   <- if(pass) "[PASS]" else "[FAIL]"
    message(sprintf("%s Score 1: %s (Expect %.4f - 1)", 
                    tag, paste(round(range(u_sub), 4), collapse=" - "), mmrc_cut))
  }
}

message("\n---------------------------------------------------")
message("TEST: AGE (Check 0-1 range)")
message("---------------------------------------------------")

# Age
u_age <- clip01(plnorm(pmax(cand$age_at_COPD, 1e-6), meanlog = ml_age, sdlog = sl_age)) # CORRECTION use the age at index. Define a new variable in cand
# when age = age_at_creation + local_time
rng <- range(u_age)
pass <- rng[1] > 0 & rng[2] < 1
tag <- if(pass) "[PASS]" else "[FAIL]"
message(sprintf("%s Age U-Range: %s", tag, paste(round(rng, 4), collapse=" - ")))

message("\n---------------------------------------------------")
message("TEST: EXACERBATIONS (Check 0-1 range)")
message("---------------------------------------------------")

# Moderate Exacerbations (Jittered)
u_ex_mod <- pit_poisson(cand$exac_mod, lambda = lambda_exac_mod)
rng_mod  <- range(u_ex_mod)
# Check if all values are strictly between 0 and 1
pass_mod <- all(u_ex_mod > 0 & u_ex_mod < 1)
tag_mod  <- if(pass_mod) "[PASS]" else "[FAIL]"

message(sprintf("%s ExacMod U-Range: %s", tag_mod, paste(round(rng_mod, 4), collapse=" - ")))


# Severe Exacerbations (Jittered)
u_ex_sev <- pit_poisson(cand$exac_sev, lambda = lambda_exac_sev)
rng_sev  <- range(u_ex_sev)
# Check if all values are strictly between 0 and 1
pass_sev <- all(u_ex_sev > 0 & u_ex_sev < 1)
tag_sev  <- if(pass_sev) "[PASS]" else "[FAIL]"

message(sprintf("%s ExacSev U-Range: %s", tag_sev, paste(round(rng_sev, 4), collapse=" - ")))


message("\n---------------------------------------------------")
message("TEST: VISUALIZATION")
message("---------------------------------------------------")
layout(matrix(c(1,2), nrow=1))

hist(u_smoking, breaks=50, col="skyblue", main="Smoking U-Vals", xlab="U")
abline(v=smk_cuts[1:2], col="red", lwd=2, lty=2)

hist(u_gold, breaks=50, col="orange", main="GOLD U-Vals", xlab="U")
abline(v=gold_cuts[1:3], col="red", lwd=2, lty=2)

layout(1)


# ==============================================================================
#  COPULA MATCHING: JOINT LIKELIHOOD CALCULATION
# ==============================================================================
#
# PURPOSE:
#  -This script takes the "Percentile Transformed" agents (U-values calculated in the previous step)
#  and calculates a "Match Score" (Log Likelihood) for each agent.
#
# HOW IT WORKS:
#  Matrix Construction:
#   - Assembles the individual U-values (Age, Sex, Smoking, etc.) into a single matrix.
#   - Ensures the column order matches the Target Correlation Matrix (Sigma) exactly.
#
# Correlation Scoring (The Copula):
#   - Uses the 'dCopula' function to calculate how well an agent's combination of characteristics
#       matches the Target's correlation structure.
#   - Example: A 80-year-old smoker with GOLD 4 gets a HIGH score (Matched Correlation).
#   - Example: A 40-year-old non-smoker with GOLD 4 gets a LOW score (Correlation Mismatch).
#
# Total Score Calculation:
#   - Combines the Correlation Score (Copula) with the Individual Trait Scores (Marginals).
#   - Result: 'll_joint' (Total Log Likelihood).
#
# OUTPUT:
#  - 'll_joint': A single number for each agent. High numbers = Good Match. Low numbers = Bad Match.
#  - These scores are converted to weights (w = exp(ll_joint)) in the next step (Raking).
#
# ==============================================================================

# ==============================================================================
# SETUP NOVELTY PARAMATERS FROM THE OUTPUT FILE PROVIDED BY LAURA (NOT A TEST)
# ==============================================================================

# Helper for categorical log-prob
log_pmf_cat <- function(x, levels, probs) {
  idx <- match(as.integer(x), as.integer(levels))
  p   <- probs[idx]
  p <- pmax(p, 1e-9) # Safety: add epsilon to avoid log(0)
  log(p)
}

TAU_CAT <- 1.0

# Obtain sigma and var_order from the NOVELTY output (this contains the joint relationships between variables)
Sigma     <- target$Sigma
var_order <- target$var_order

# Apply names to Sigma so we can index it by name later
dimnames(Sigma) <- list(var_order, var_order)

# ==============================================================================
# COMBINE ALL "U" VALUES FROM EPIC TEST
# ==============================================================================

message("\n---------------------------------------------------")
message("COMBINE 'U' VALUES for EPIC")
message("---------------------------------------------------")
# TEST OBJECTIVE: DATA STRUCTURE VALIDITY
# What this tests:
# Do we have U-values for every single agent? (No NAs allowed)
# Do we have the correct number of columns before alignment?
# Are there any impossible values (e.g., exactly 0 or 1) that would break the math?
#
# Why it matters:
# The 'dCopula' function does not accept any NA or infinite value as it will cause the entire simulation run to crash or return invalid weights for everyone.

# Build the matrix where we combine all the 'U values from EPIC
U_epic <- cbind(
  sex_baseline                  = u_sex,
  ageCOPD_baseline              = u_age, # This should be updated after ageCOPD_baseline is change to age at index
  smoking_baseline              = u_smoking,
  GOLDgrade_baseline            = u_gold,
  mMRC_baseline                 = u_mmrc,
  exacerbations_moderate_baseline   = u_ex_mod,
  exacerbations_severeplus_baseline = u_ex_sev
)

message(sprintf("Matrix Size Created: %d rows x %d cols\n", nrow(U_epic), ncol(U_epic)))

# NA Check
na_count <- sum(is.na(U_epic))
if (na_count > 0) stop("FAIL: Matrix contains NAs. Check your 'u_' variables.")


# ==============================================================================
# COPULA ALIGNMENT TEST
# ==============================================================================

message("\n---------------------------------------------------")
message("PART B: ALIGNMENT & CORRELATION")
message("---------------------------------------------------")
# TEST OBJECTIVE: VARIABLE MATCHING
# What this tests:
# Does the simulation have every variable the Target Matrix expects?
# Are the columns in the exact same order as the Target Matrix?
#
# Why it matters:
# The Copula math relies purely on column order and does not recognize any variable names
# If Column 1 in the matrix is 'Age' but Column 1 in the Target is 'Sex', the model will try to correlate Age with Sex, producing incorrect weights

# Check Variables
sim_vars    <- colnames(U_epic)
target_vars <- var_order 

missing <- setdiff(target_vars, sim_vars)
extra   <- setdiff(sim_vars, target_vars)

message("Missing Variables: ", if(length(missing)==0) "None (PASS)" else paste(missing, collapse=", "))
message("Extra Variables:   ", if(length(extra)==0)   "None (PASS)" else paste(extra, collapse=", "))

if(length(missing) > 0) stop("FAIL: Simulation is missing variables required by Target Sigma!")

# Align U to Sigma Order
U_final <- U_epic[, target_vars, drop = FALSE]

# Prepare Correlation Matrix (R)
# Now this line will work because Sigma has names
R_target <- cov2cor(as.matrix(Sigma[target_vars, target_vars]))
R_final  <- as.matrix(Matrix::nearPD(R_target, corr = TRUE)$mat)

message("Alignment Status:  [PASS] U matrix matches Sigma structure.\n")


# ==============================================================================
# COPULA CALCULATION TEST
# ==============================================================================

message("\n---------------------------------------------------")
message("JOINT LIKELIHOOD CALCULATION")
message("---------------------------------------------------")
# TEST OBJECTIVE: CORRELATION FIT
# What this tests:
# The "Structure" Score: How well does this agent's combination of traits match the Target correlations?
#   - High Score (+): Traits fit the correlation (e.g. Old + Sick).
#   - Low Score (-): Traits violate correlation (e.g. Young + Sick).
# Mathematical Stability: Are the resulting densities finite numbers (not -Inf)?
#
# Why it matters:
# This is the core of the method as it identifies agents that represent the NOVELTY cohort.
# If this produces -Inf, it means the simulation produced an "impossible" agent.

# Define Copula Object
gcop <- copula::normalCopula(param = copula::P2p(R_final), dim = ncol(R_final), dispstr = "un")

# Calculate Density
message("Calculating Copula Density...")
ll_cop <- copula::dCopula(U_final, gcop, log = TRUE)

# Check Results
range_cop <- range(ll_cop)
message(sprintf("Copula LL Range:   %.2f to %.2f", range_cop[1], range_cop[2]))
message("                   (Expect: Finite numbers, typically -20 to +10)\n")

if (any(!is.finite(ll_cop))) {
  warning("FAIL: Infinite/NaN values found! (Weights will collapse)")
} else {
  message("Copula Values:     [PASS] All values are finite.\n")
}

# ==============================================================================
# TOTAL SCORE CALCULATION TEST
# ==============================================================================

message("\n---------------------------------------------------")
message("TOTAL SCORE (FINAL WEIGHT INPUT)")
message("---------------------------------------------------")
# TEST OBJECTIVE: FINAL WEIGHT PREPARATION
# What this tests:
# The combination of "Correlation Fit" (Part C) and "Individual Probability" (Marginals).
# Final Score = ll_cop + ll_marginals
#
# Why it matters:
# This 'Total Joint LL' is the raw number that gets converted into the final weight.
# weight = exp(Total Joint LL)
# If this number is valid, the weighting engine is ready to run.

# Category Log-Likelihoods
ll_cat <- 
  log_pmf_cat(cand$smoking_status, levels = lev_smoking, probs = prob_smoking) +
  log_pmf_cat(cand$gold,           levels = lev_gold,    probs = prob_gold)

# Marginal Log-Likelihoods
ll_marg <-
  dlnorm(pmax(cand$age_at_COPD, 1e-6), meanlog = ml_age, sdlog = sl_age, log = TRUE) +
  dbinom(cand$female, size = 1, prob = p_sex, log = TRUE) +
  dbinom(cand$dyspnea, size = 1, prob = p_mmrc, log = TRUE) +
  dpois(cand$exac_mod, lambda = lambda_exac_mod, log = TRUE) +
  dpois(cand$exac_sev, lambda = lambda_exac_sev, log = TRUE) +
  # Add Weighted Categories
  (TAU_CAT * ll_cat)

# Total Score
ll_joint <- ll_cop + ll_marg

# --- RESULTS ---
message(sprintf("Marginal LL Range: %.2f to %.2f", min(ll_marg), max(ll_marg)))
message(sprintf("Total Joint LL:    %.2f to %.2f", min(ll_joint), max(ll_joint)))

message("\n[SUCCESS] Weights are ready to be calculated from 'Total Joint LL'.")

# ==============================================================================
# CORRELATION RECOVERY TEST
# ==============================================================================

message("\n---------------------------------------------------")
message("CORRELATION RECOVERY CHECK")
message("---------------------------------------------------")
# DOCUMENTATION:
# This compares the empirical correlation of your simulated agents (U_final) against the Target Correlation Matrix (Sigma).
#
# Goal: We want the difference to be small (closer to 0).
# Interpretation:
# - Small difference (< 0.1): Your simulation naturally mimics the target structure.
# - Large difference (> 0.3): Your simulation has structural differences (e.g., 
#     maybe in EPIC, age and disease aren't as strongly linked as in NOVELTY).
#     The Copula weights will have to work very hard to fix this.

# Calculate Empirical Correlation of the Simulation
# (Using Spearman to capture rank relationships, or Pearson for raw U values)
emp_R <- cor(U_final, method = "pearson")

# Prepare Target Correlation (Subset to match current vars)
target_R <- cov2cor(as.matrix(Sigma[colnames(U_final), colnames(U_final)]))

# Calculate Absolute Difference Matrix
diff_R <- abs(emp_R - target_R)

# Report Metrics
max_diff <- max(diff_R)
avg_diff <- mean(diff_R)

# Find the worst pair to report specifically
worst_idx <- which(diff_R == max_diff, arr.ind = TRUE)[1, ]
var1 <- rownames(diff_R)[worst_idx[1]]
var2 <- colnames(diff_R)[worst_idx[2]]

message(sprintf("Max Correlation Diff:  %.4f (Between %s & %s)", max_diff, var1, var2))
message(sprintf("Avg Correlation Diff:  %.4f", avg_diff))

if (max_diff > 0.3) {
  warning("FAIL: Simulation correlation structure deviates significantly from Target (> 0.3).")
} else {
  message("STATUS: [PASS] Simulation structure roughly matches Target.")
}

message("\n=== FULL DIAGNOSTIC SUITE COMPLETE ===")