# NOVLETY Regression Model for Exacerbations

temp1 <- df_final %>% filter(grp01bl == "COPD") %>%
  mutate(Year = 1, usubjid = usubjid, LAMA = LAMA_y1, LABA = LABA_y1, ICS = ICS_y1, SABA = SABA_y1, FU_max1 = pmin(1, pmax(as.numeric(FU), 0)), rate_m = m1, rate_s = s1, rate_ms = m1s1)
temp2 <- df_final %>% filter(grp01bl == "COPD") %>%
  mutate(Year = 2, usubjid = usubjid, LAMA = LAMA_y2, LABA = LABA_y2, ICS = ICS_y2, SABA = SABA_y2, FU_max1 = pmin(1, pmax(as.numeric(FU)-1, 0)), rate_m = m2, rate_s = s2, rate_ms = m2s2)
temp3 <- df_final %>% filter(grp01bl == "COPD") %>%
  mutate(Year = 3, usubjid = usubjid, LAMA = LAMA_y3, LABA = LABA_y3, ICS = ICS_y3, SABA = SABA_y3, FU_max1 = pmin(1, pmax(as.numeric(FU)-2, 0)), rate_m = m3, rate_s = s3, rate_ms = m3s3)
df_reg <- rbind(temp1, temp2, temp3) %>% filter(FU_max1>0.01)

df_reg <- df_reg %>% #group_by(usubjid) %>%
  mutate(#FU_max3 = as.numeric(pmin(FU, 3)),
    #  rate_m = (m1+m2+m3)/FU_max3,
    #  rate_s = (s1+s2+s3)/FU_max3,
    #  rate_ms = (m1s1+m2s2+m3s3)/FU_max3,
    med = case_when(SABA == 0 & LAMA == 0 & LABA == 0 & ICS == 0 ~ 0,
                    SABA == 1 & LAMA == 0 & LABA == 0 & ICS == 0 ~ 1,
                    SABA == 0 & LAMA == 1 & LABA == 0 & ICS == 0 ~ 4,
                    SABA == 0 & LAMA == 1 & LABA == 1 & ICS == 0 ~ 6,
                    SABA == 0 & LAMA == 1 & LABA == 1 & ICS == 1 ~ 14,

                    SABA == 1 & LAMA == 1 & LABA == 1 & ICS == 1 ~ 14,
                    SABA == 0 & LAMA == 0 & LABA == 1 & ICS == 1 ~ 6,
                    SABA == 1 & LAMA == 0 & LABA == 0 & ICS == 1 ~ 4,
                    SABA == 1 & LAMA == 1 & LABA == 1 & ICS == 0 ~ 6,
                    SABA == 1 & LAMA == 0 & LABA == 1 & ICS == 1 ~ 6,
                    SABA == 1 & LAMA == 1 & LABA == 0 & ICS == 0 ~ 4,
                    SABA == 0 & LAMA == 0 & LABA == 0 & ICS == 1 ~ 4,
                    SABA == 1 & LAMA == 0 & LABA == 1 & ICS == 0 ~ 4,
                    SABA == 1 & LAMA == 0 & LABA == 1 & ICS == 0 ~ 4,
                    TRUE ~ NA_real_),
    med = ifelse(is.na(med), 20, med),
    med = factor(med, levels = c(0,1,4,6,14,20), ordered = FALSE)
  )

summary_rates <- lapply(df_reg[c("rate_m", "rate_s", "rate_ms")], summary)
summary_med <- table(df_reg$Year, df_reg$med)

library(glmmTMB)
# Full cohort
#model_m <- glmmTMB(rate_m ~  med + offset(log(as.numeric(FU_max1))) + (1|usubjid),data = df_reg[df_reg$med!="20",], family=nbinom2, zi = ~1 + med)
#model_s <- glmmTMB(rate_s ~  med + offset(log(as.numeric(FU_max1))) + (1|usubjid),data = df_reg[df_reg$med!="20",], family=nbinom1, zi = ~1 + med)
#model_ms <- glmmTMB(rate_ms ~  med + offset(log(as.numeric(FU_max1))) + (1|usubjid),data = df_reg[df_reg$med!="20",], family=nbinom1, zi = ~1 + med)
# Non-asthma cohort and both spirometry-confirmed cohorts
model_m <- glmmTMB(rate_m ~  med + offset(log(as.numeric(FU_max1))) + (1|usubjid),data = df_reg[df_reg$med!="20",], family=nbinom2)
model_s <- glmmTMB(rate_s ~  med + offset(log(as.numeric(FU_max1))) + (1|usubjid),data = df_reg[df_reg$med!="20",], family=nbinom2)
model_ms <- glmmTMB(rate_ms ~  med + offset(log(as.numeric(FU_max1))) + (1|usubjid),data = df_reg[df_reg$med!="20",], family=nbinom2)
