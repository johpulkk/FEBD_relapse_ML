## ===========================
## Time-to-event analysis 
## (the example uses the Finnish data)
## Author: M.D., Ph.D., Johannes Lieslehto
## email: johannes.lieslehto@ki.se
## ===========================


## ---- Required Packages ----
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, stringr, survival, survminer, ggplot2, intsurv)

###----Preprocessing the data----
# Keep original dates
study_end <- as.Date("2018-12-31")  # 2021-12-31 for the Swedish internal validation

# tulopvm = date of bipolar relapse
# kuolpv = date of death
# ced = start date of follow-up

fi <- fi_data2 %>%
  mutate(
    # times in days from cohort entry 'ced'
    t_relapse_raw = as.numeric(tulopvm - ced),                    # may be NA
    t_relapse     = as.numeric(pmin(tulopvm, study_end, na.rm=TRUE) - ced),  # relapse or admin censor
    t_death       = as.numeric(kuolpv - ced),                     # may be NA
    # analysis time = earliest of relapse/censor and death
    time_days     = pmin(t_relapse, t_death, na.rm = TRUE),
    # event = relapse occurs before death and before study end
    event_relapse = as.integer(!is.na(tulopvm) &
                                 (is.na(kuolpv) | tulopvm <= kuolpv) &
                                 tulopvm <= study_end),
    time_years    = time_days / 365.25
  )

# Basic sanity checks
stopifnot(all(is.finite(fi$time_days)))
stopifnot(all(fi$time_days >= 0))

# Quintiles based on cutpoints from the development (train):
qq <- quantile(train$rel_mod_cal, probs = seq(0, 1, by = 0.2), na.rm = TRUE)
qq <- unique(qq)  # guard against duplicates
qq[1] <- -Inf; qq[length(qq)] <- Inf

fi$quintile <- cut(fi$rel_mod_cal,
                   breaks = qq,
                   include.lowest = TRUE,
                   labels = paste0("Q", seq_len(length(qq) - 1)))

##----Plotting the risk-stratified KM curves----
km_fit <- survfit(Surv(time_years, event_relapse) ~ quintile, data = fi)

pdf("Cum_relapse_over_Quintiles.pdf", width = 8, height = 7.78)
ggsurvplot(
  km_fit,
  data = fi,
  fun = "event",              # 1 - S(t): note this is NOT a CIF with competing risks
  risk.table = TRUE,
  pval = FALSE,
  palette = c("#DF8F44FF", "#00A1D5FF", "#B24745FF", "#79AF97FF", "#6A6599FF"),
  conf.int = FALSE,
  legend.labs = levels(fi$quintile),
  surv.median.line = "hv",
  ylim = c(0, 1),
  xlim = c(0, 23),
  break.time.by = 3,
  censor = TRUE,
  risk.table.fontsize = 3.5,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = TRUE
)
dev.off()

### ----C-index (Harrell)----
set.seed(1)
c_index_value <- intsurv::cIndex(time = fi$time_days,
                                 event = fi$event_relapse,
                                 risk_score = fi$rel_mod_cal)
print(c_index_value)

# 95% CI by bootstrap
B <- 1000
ids <- seq_len(nrow(fi))
boot_c <- replicate(B, {
  ii <- sample(ids, replace = TRUE)
  intsurv::cIndex(time = fi$time_days[ii],
                  event = fi$event_relapse[ii],
                  risk_score = fi$rel_mod_cal[ii])["index"]
})
print(quantile(boot_c, c(.025, .975)))

print(summary(km_fit, times = seq(0, 20, 1)), digits = 6)

####----Cox Regression (cause-specific hazard)----
coxph3b <- coxph(Surv(time_years, event_relapse) ~ quintile, data = fi)
summary(coxph3b)
