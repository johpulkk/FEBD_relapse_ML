## ================================
## Fairness analysis: subgroup vs subgroup
## Author: M.D., Ph.D., Johannes Lieslehto
## email: johannes.lieslehto@ki.se
## ================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  dplyr, pROC, ggplot2,
)


## ---- Helper Functions: ----
clamp <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)

as01 <- function(y) {
  y2 <- suppressWarnings(as.integer(as.character(y)))
  if (any(is.na(y2))) stop("Outcome must be 0/1 or factor of '0'/'1'.")
  y2
}

# rel_mod_cal = calibrated model predictions
# outcome = 2-year hospital due to bipolar relapse

compute_metrics <- function(df) {
  stopifnot(all(c("outcome","rel_mod_cal") %in% names(df)))
  y <- as01(df$outcome)
  p <- clamp(df$rel_mod_cal)
  # Metrics
  brier <- mean((y - p)^2)
  lp <- qlogis(p)
  citl <- coef(glm(y ~ 1, family = binomial(), offset = lp))[1]
  slope <- coef(glm(y ~ lp, family = binomial()))["lp"]
  rocobj <- pROC::roc(y, p, quiet = TRUE)
  auc <- as.numeric(pROC::auc(rocobj))
  ci  <- pROC::ci.auc(rocobj)
  list(n = length(y), prev = mean(y), brier = brier,
       citl = citl, slope = slope, auc = auc, auc_ci = ci, roc = rocobj)
}

perm_test <- function(dfA, dfB, nIter = 1000, metric = c("brier","citl","slope","auc"), seed = 42) {
  metric <- match.arg(metric)
  pooled <- bind_rows(dfA %>% select(outcome, rel_mod_cal),
                      dfB %>% select(outcome, rel_mod_cal))
  nA <- nrow(dfA)
  mfun <- switch(metric,
                 brier = function(d) { y <- as01(d$outcome); p <- clamp(d$rel_mod_cal); mean((y - p)^2) },
                 citl  = function(d) { y <- as01(d$outcome); p <- clamp(d$rel_mod_cal); lp <- qlogis(p);
                 coef(glm(y ~ 1, family = binomial(), offset = lp))[1] },
                 slope = function(d) { y <- as01(d$outcome); p <- clamp(d$rel_mod_cal); lp <- qlogis(p);
                 coef(glm(y ~ lp, family = binomial()))[2] },
                 auc   = function(d) { y <- as01(d$outcome); p <- clamp(d$rel_mod_cal);
                 as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE))) }
  )
  obs <- mfun(dfA) - mfun(dfB)
  set.seed(seed)
  diffs <- numeric(nIter)
  for (i in seq_len(nIter)) {
    idx <- sample.int(nrow(pooled), size = nA, replace = FALSE)
    diffs[i] <- mfun(pooled[idx, ]) - mfun(pooled[-idx, ])
  }
  pval <- mean(abs(diffs) >= abs(obs))
  list(obs = obs, pval = pval, diffs = diffs)
}

plot_auc_compare <- function(mA, mB, labelA, labelB, title = NULL, file = NULL) {
  labA <- sprintf("%s (AUC=%.3f, 95%%CI=%.3f–%.3f)",
                  labelA, mA$auc, mA$auc_ci[1], mA$auc_ci[3])
  labB <- sprintf("%s (AUC=%.3f, 95%%CI=%.3f–%.3f)",
                  labelB, mB$auc, mB$auc_ci[1], mB$auc_ci[3])
  g <- pROC::ggroc(setNames(list(mA$roc, mB$roc), c(labA, labB)),
                   size = 0.7, legacy.axes = TRUE) +
    ggplot2::labs(title = title, x = "1 - Specificity", y = "Sensitivity") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 1, y = 0, yend = 1),
                          color = "darkgrey", linetype = "dashed") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position = c(0.7, 0.25),
                   legend.text = ggplot2::element_text(size = 6))
  if (!is.null(file)) {
    grDevices::pdf(file, width = 6.3, height = 6); print(g); grDevices::dev.off()
  }
  print(g)
}

fairness_analysis <- function(immi_df, muut_df, labelA, labelB, nIter = 1000,
                              plot_title = NULL, plot_file = NULL) {
  # keep only needed cols and ensure types
  A <- immi_df %>% select(outcome, rel_mod_cal)
  B <- muut_df %>% select(outcome, rel_mod_cal)
  
  mA <- compute_metrics(A)
  mB <- compute_metrics(B)
  
  # permutation tests for differences (two-sided)
  pt_auc   <- perm_test(A, B, nIter, "auc")
  pt_brier <- perm_test(A, B, nIter, "brier")
  pt_citl  <- perm_test(A, B, nIter, "citl")
  pt_slope <- perm_test(A, B, nIter, "slope")
  
  # DeLong test for AUC difference
  delong <- pROC::roc.test(mA$roc, mB$roc, paired = FALSE)
  
  cat("\n==============================\n")
  cat(sprintf("Fairness analysis: %s vs %s\n", labelA, labelB))
  cat("==============================\n")
  print(data.frame(
    group       = c(labelA, labelB),
    n           = c(mA$n, mB$n),
    prevalence  = round(c(mA$prev, mB$prev), 4),
    AUC         = round(c(mA$auc, mB$auc), 4),
    Brier       = signif(c(mA$brier, mB$brier), 4),
    CITL        = signif(c(mA$citl, mB$citl), 4),
    Slope       = signif(c(mA$slope, mB$slope), 4)
  ), row.names = FALSE)
  
  cat(sprintf("\nPermutation p-values (two-sided):\nΔAUC = %.4f, ΔBrier = %.4f, ΔCITL = %.4f, ΔSlope = %.4f\n",
              pt_auc$pval, pt_brier$pval, pt_citl$pval, pt_slope$pval))
  cat("\nDeLong test for AUC difference:\n"); print(delong)
  
  # ROC comparison plot
  plot_auc_compare(mA, mB, labelA, labelB, title = plot_title, file = plot_file)
  
  invisible(list(metrics_A = mA, metrics_B = mB,
                 pvals = list(auc_perm = pt_auc$pval, brier_perm = pt_brier$pval,
                              citl_perm = pt_citl$pval, slope_perm = pt_slope$pval),
                 delong = delong))
}

## ================================
## Subgroup comparisons
## Assumes df_2y has: outcome, rel_mod_cal, edu_cat.1, Kon, population_group
## ================================

# 1) Population group: Non–Swedish-born vs Swedish-born
# df_2y = held-out validation sample
immi <- subset(df_2y, df_2y$population_group != "Swedish-born")
muut <- subset(df_2y, df_2y$population_group == "Swedish-born")
fairness_analysis(immi, muut,
                  labelA = "Non–Swedish-born", labelB = "Swedish-born",
                  nIter = 1000,
                  plot_title = "ROC by Population Group",
                  plot_file  = "ROC_fairness_population_group.pdf")

#2) Education: edu_cat.1 == "1" vs others
immi <- subset(df_2y, df_2y$edu_cat.1 == "1")
muut <- subset(df_2y, df_2y$edu_cat.1 != "1")
fairness_analysis(immi, muut,
                   labelA = "Education cat 1", labelB = "Other education",
                   nIter = 1000,
                   plot_title = "ROC by Education",
                   plot_file  = "ROC_fairness_education.pdf")

#3) Sex (Kon): Kon != "1" vs Kon == "1"  (adjust labels to your coding)
immi <- subset(df_2y, df_2y$Kon != "1")
muut <- subset(df_2y, df_2y$Kon == "1")
fairness_analysis(immi, muut,
                   labelA = "Kon != 1", labelB = "Kon == 1",
                   nIter = 1000,
                   plot_title = "ROC by Sex",
                   plot_file  = "ROC_fairness_sex.pdf")
