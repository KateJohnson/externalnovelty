# ---------------- EPIC population test ----------------

## This is to confirm only patients with a confirmed diagnosis of COPD is captured and simulated patients have at least 1 year of follow up post-diagnosis to be eliglbe for matching

### Identify patients that fufil criteria
cand <- all_events2 %>%
  filter(diagnosis > 0, gold > 0, FEV1_baseline > 0) %>%
  group_by(id) %>%
  mutate(dx_time = min(local_time[diagnosis > 0], na.rm = TRUE)) %>%
  filter(local_time >= dx_time + 1) %>%
  slice_min(local_time, n = 1, with_ties = FALSE) %>%
  ungroup()

### Convert data-type in preparation for copula 
cand <- cand %>%
  mutate(
    smoking = suppressWarnings(as.integer(as.character(smoking_status))),
    pack_years         = suppressWarnings(as.numeric(pack_years)),
    gold               = as.integer(gold),
    medication_status  = as.integer(medication_status),
    exac_mod           = as.integer(exac_history_n_moderate),
    exac_sev           = as.integer(exac_history_n_severe_plus)
  ) %>%
  mutate(
    smoking_status = dplyr::case_when(
      smoking == 1L                                 ~ 1L,  # current
      smoking == 0L & !is.na(pack_years) &
        pack_years >  0.1                           ~ 2L,  # former
      smoking == 0L & !is.na(pack_years) &
        pack_years <= 0.1                           ~ 0L,  # never
      TRUE                                          ~ NA_integer_
    )
  )

### Test if the correct patient and time-point is identified
cand_test_id <- all_events2$id[
  which(all_events2$diagnosis == 1 & all_events2$gold > 0)[1]
]
cand_test_history <- all_events2 %>%
  filter(id == cand_test_id & diagnosis==1)

cand_id_EPIC<- cand %>%
  filter(id == cand_test_id)

# ---------------- Helper function tests ----------------

## Categorical variable transformation

pit_categorical_from_probs <- function(values, level_order, level_probs) {
  idx <- match(as.integer(values), as.integer(level_order))
  if (anyNA(idx)) stop("pit_categorical_from_probs: some values not in level_order.")
### Ensure probabilities are numeric and normalized
  p <- as.numeric(level_probs)
  if (any(!is.finite(p) | p < 0)) {
    stop("pit_categorical_from_probs: level_probs must be finite and non-negative.")
  }
  p <- p / sum(p)
#### Cumulative distribution over levels (in level_order)
  cum_all   <- cumsum(p)                 # e.g. [p1, p1+p2, p1+p2+p3]
  lower_all <- c(0, head(cum_all, -1))   # e.g. [0, p1, p1+p2]
### Map each observation to its [lower, upper] interval
  lower <- lower_all[idx]
  upper <- cum_all[idx]
  u <- stats::runif(length(values), min = lower, max = upper)
  clip01(u)
}

## Test on smoking variable
u_smoking <- pit_categorical_from_probs(cand$smoking_status, level_order = lev_smoking, level_probs = prob_smoking)

### Are values between 0-1
smoking_test <- all(u_smoking >= 0 & u_smoking <= 1)
cat("All values in [0,1]?:", in_range, "\n")

### Values must also span between the entirety of 0 and 1 (not for example 0-0.2)
cat("Range of u_smoking:", range(u_smoking), "\n\n")

### Expect to observe 3 clusters spaced out between 0-1 to account for the 3 categories (current, former, non-smoker)
km_smoking <- kmeans(u_smoking, centers = 3, nstart = 20)
cat("Cluster centers (sorted):\n")
print(sort(km_smoking$centers))



## Test on GOLD variable
u_gold <- pit_categorical_from_probs(cand$gold,level_order = lev_gold, level_probs = prob_gold)

### Are values between 0-1
gold_test <- all(u_gold >= 0 & u_gold <= 1)
cat("All values in [0,1]?:", gold_test, "\n")

### Values must also span between the entirety of 0 and 1 (not for example 0-0.2)
cat("Range of u_gold:", range(u_gold), "\n\n")

### Expect to observe 3 clusters spaced out between 0-1 to account for the 3 categories (GOLD 1-GOLD 4)
km_gold <- kmeans(u_gold, centers = 4, nstart = 20)

cat("Cluster centers (sorted):\n")
print(sort(km_gold$centers))



## Binary variable transformation

pit_bernoulli <- function(y01, p) {
  out <- rep(NA_real_, length(y01))
  i0 <- y01 == 0; i1 <- y01 == 1
  out[i0] <- runif(sum(i0, na.rm = TRUE), 0, 1 - p)
  out[i1] <- runif(sum(i1, na.rm = TRUE), 1 - p, 1)
  clip01(out)
}

## Test on sex variable
u_sex  <- pit_bernoulli(cand$female, p = p_sex)

### Are values between 0-1
sex_test <- all(u_sex >= 0 & u_sex <= 1)
cat("All values in [0,1]?:", sex_test, "\n")

### Values must also span between the entirety of 0 and 1 (not for example 0-0.2)
cat("Range of u_sex (non-NA):", range(u_sex), "\n\n")

### Expect to observe 2 clusters spaced out between 0-1 to account for the males and females
km_sex <- kmeans(u, centers = 2, nstart = 20)

cat("Cluster centers (sorted):\n")
print(round(sort(km_sex$centers), 2))



## Test on mmRC variable

u_mmrc<- pit_bernoulli(cand$dyspnea, p = p_mmrc)

### Are values between 0-1
mmrc_test <- all(u_sex >= 0 & u_sex <= 1)
cat("All values in [0,1]?:", mmrc_test, "\n")

### Values must also span between the entirety of 0 and 1 (not for example 0-0.2)
cat("Range of u_sex (non-NA):", range(u_sex), "\n\n")

### Expect to observe 2 clusters spaced out between 0-1 to account for the males and females
km_mmrc <- kmeans(u, centers = 2, nstart = 20)

cat("Cluster centers (sorted):\n")
print(round(sort(km$centers), 2))
