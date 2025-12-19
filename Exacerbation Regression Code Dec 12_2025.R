# ==============================================================================
# Exacerbation Regression Analysis for EPIC
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(glmmTMB)
  library(broom.mixed)
  library(tibble)
})

# ==============================================================================
# LOAD DATA 
# ==============================================================================

# 'combined_histories' contains the full lifetime of every patient (creation to death).
#  We cannot only use this dataset because it includes years of history and does not contain the index date 
if (!exists("combined_histories")) {
  if (file.exists("combined_histories.rds")) {
    message("Loading 'combined_histories' from RDS file...")
    combined_histories <- readRDS("combined_histories.rds")
  } else {
    stop("Error: 'combined_histories' is missing.")
  }
}

#  We will also need to use epic_selected which contains the list of patients chosen by the Copula.
#  This dataset also holds the index date for each patient
if (!exists("epic_selected")) {
  if (file.exists("epic_selected.rds")) {
    message("Loading 'epic_selected' from RDS file...")
    epic_selected <- readRDS("epic_selected.rds")
  } else {
    stop("Error: 'epic_selected.rds' missing.")
  }
}

# ==============================================================================
# LINKING DATASETS 
# ==============================================================================

# Look up the selected patients using the 'id' + '.run_id' columns from epic_selected
# The "id" represents the patient 'id' and the '.run_id' represents which EPIC run the patient was in
# Because we run EPIC many times to find the appropriate match of patients,
# looking for the 'id' + '.run_id' combination is to identify each unique patient included in the analysis
index_lookup <- epic_selected %>%
  select(id, .run_id, index_time = local_time) %>%
  mutate(
    id = as.character(id),
    .run_id = as.integer(.run_id)
  )

# Prepare the combined_histories dataset to join with epic_selected
  history_EPIC <- combined_histories %>%
  # Ensure variable types match
  mutate(
    id = as.character(id),
    .run_id = as.integer(.run_id),
    local_time = as.numeric(local_time),
    medication_status = suppressWarnings(as.numeric(medication_status)),
    # Ensure event counts are numeric for calculation later
    exac_history_n_moderate = as.numeric(exac_history_n_moderate),
    exac_history_n_severe_plus = as.numeric(exac_history_n_severe_plus)
  ) %>%
  
  # We link the history file to the lookup dataset using 'id' and the '.run_id'
  # This adds the column 'index_time' to every row of the patient's history.
  inner_join(index_lookup, by = c("id", ".run_id")) %>%
  
  # We delete all rows where 'local_time' is less than 'index_time'.
  # This removes the Pre-Index History so we don't analyze it.
  filter(local_time >= index_time) %>%
  
  # Create a unique ID for grouping and sort chronologically
  mutate(unique_id = paste0(id, "_run", .run_id)) %>%
  arrange(unique_id, local_time)

# ==============================================================================
# CALCULATE NUMBER OF EXACERBATION EVENTS
# ==============================================================================

exac_data_EPIC <- history_EPIC %>%
  group_by(unique_id) %>%
  mutate(
    # --- Calculate Number of New Events ---
    # The exac history data is cumulative (e.g., 5 events, then 6 events).
    # We use diff() to determine the change (ex. 6 - 5 = 1 new event).
    events_mod_new = c(0, diff(exac_history_n_moderate)),
    events_sev_new = c(0, diff(exac_history_n_severe_plus)),
    
    # --- Calculate Follow-up Time ---
    # We calculate the time elapsed since the previous record.
    followup_time  = c(0, diff(local_time))
  ) %>%
  
  # --- Remove first row ---
  # The first row for every patient is their status at the Index Date.
  # Because of the logic above, it has 0 follow-up and 0 new events.
  # We delete this row so we only analyze the follow-up period (Index -> Future).
  slice(-1) %>%
  ungroup()

# ==============================================================================
#  PREPARE DATASET FOR MODEL
# ==============================================================================

# Group 0 & 1 together to be the reference (No therapy or SABA only); 
# Group 4= Monotherapy (LAMA/LABA/ICS), Group 6= Dual therapy (LAMA + LABA, ICS + LABA), 14 (ICS + LAMA + LABA)
target_meds <- c("1", "4", "6", "14")

analysis_EPIC<- exac_data_EPIC  %>%
  mutate(
    med_group = case_when(
      medication_status %in% c(0, 1) ~ "1",
      medication_status == 4         ~ "4",
      medication_status == 6         ~ "6",
      medication_status == 14        ~ "14",
      TRUE                           ~ NA_character_
    )
  ) %>%
  filter(!is.na(med_group)) %>%
  mutate(med_group = factor(med_group, levels = target_meds)) %>%
  
  # Matches the NOVELTY analysis of splitting patient history into separate yearly
  # intervals (Year 1, 2, 3) rather than summarizing a single lifetime value.
  # Creates separate rows if a patient changes medication mid-year, ensuring follow-up 
  # is assessed using the correct drug.
  mutate(analysis_year = floor(local_time)) %>%
  group_by(unique_id, analysis_year, med_group) %>%
  summarise(
    events_mod = sum(events_mod_new, na.rm = TRUE), # Sum of new events
    events_sev = sum(events_sev_new, na.rm = TRUE), # Sum of new events
    exposure   = sum(followup_time, na.rm = TRUE),  # Total time exposed
  ) %>%
  # Removes rows with neglibible follow-up duration
  filter(exposure > 0.001)

# Check: Ensure we have reasonable numbers before running the model
cat("\n--- QC: Summary Data ---\n")
print(analysis_EPIC %>% group_by(med_group) %>% ## CORRECTION: Change from analysis cohort to analysis_EPIC
        summarise(Patients = n_distinct(unique_id), 
                  Total_Years = round(sum(exposure), 1), 
                  Total_Mod_Exac = sum(events_mod)))

# ==============================================================================
# MODEL EXABERATION RATE (PER PATIENT PER YEAR)
# ==============================================================================

# To match the NOVELTY regression model
message("\nRunning Moderate Exacerbation Model...")
mod_exac_results <- glmmTMB(
  events_mod ~ med_group + offset(log(exposure)) + (1 | unique_id),
  data = analysis_EPIC,
  family = nbinom1
)

message("Running Severe Exacerbation Model...")
sev_exac_results <- glmmTMB(
  events_sev ~ med_group + offset(log(exposure)) + (1 | unique_id),
  data = analysis_EPIC,
  family = nbinom1 
)

# ==============================================================================
# RESULTS OUTPUT
# ==============================================================================

# Extract and format Rate Ratios (RR)
format_results <- function(model, label) {
  broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    mutate(
      Rate_Ratio = exp(estimate),   # Convert Log-Odds to Rate Ratio
      RR_Low     = exp(conf.low),
      RR_High    = exp(conf.high)
    ) %>%
    select(Term = term, Rate_Ratio, RR_Low, RR_High, P_Value = p.value) %>%
    mutate(Outcome = label) %>%
    # Rounding for readability
    mutate(across(c(Rate_Ratio, RR_Low, RR_High), \(x) round(x, 3))) %>%
    mutate(P_Value = format.pval(P_Value, digits = 3, eps = 0.001))
}

results <- bind_rows(
  format_results(mod_exac_results, "Moderate Exacerbations"),
  format_results(sev_exac_results , "Severe Exacerbations")
)

cat("\n=== FINAL RESULTS ===\n")
print(results)
