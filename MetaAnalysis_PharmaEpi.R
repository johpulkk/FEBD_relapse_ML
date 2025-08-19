## ===========================
## Meta-analysis on pharmacoepidemiological
## results, stratified by the risk prediction:
## Author: M.D., Ph.D., Johannes Lieslehto
## email: johannes.lieslehto@ki.se



## ---- Required Packages ----
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  metafor, dplyr, ggplot2,
  forcats
)


## --- Helper: run fixed-effect meta per exposure within one risk group ---
fe_meta_by_exposure <- function(dat, group_value) {
  dat_g <- dat %>%
    filter(Machine_Learning_Prediction == group_value) %>%
    arrange(t1)
  
  # group by exposure index/name, run FE meta, return tidy rows
  res <- dat_g %>%
    group_by(t1, exp) %>%
    group_modify(function(.x, .y) {
      fit <- metafor::rma(yi = .x$TE, sei = .x$seTE, method = "FE")
      tibble::tibble(
        Exposure    = .y$exp,
        boxOdds     = exp(coef(fit)[1]),
        boxCILow    = exp(fit$ci.lb),
        boxCIHigh   = exp(fit$ci.ub),
        P_value     = fit$pval,
        k           = nrow(.x)
      )
    }) %>%
    ungroup() %>%
    arrange(t1) %>%
    mutate(
      Stimulus2 = Exposure,                 # left axis labels
      Stimulus  = Exposure,                 # used for the secondary axis below
      Lab       = paste0("k=", k),          # right/secondary axis labels
      `Machine Learning Prediction` = group_value,
      rowid = row_number()
    )
  
  res
}

# rr = Cox estimates from within-individual analyses per sample (stratified by the predicted risk).
# Data has to contain:
# TE = risk estimate (parameter estimate)
# seTE = standard error
# exp = medication name
# stud = cohort name (fin or swe)
# t1 = order to meta-analyse the medications
# Machine_Learning_Prediction = predicted risk group ("High risk" or "Low risk") 

## --- Compute results for High and Low risk groups ---
df_high <- fe_meta_by_exposure(rr, "High risk")
df_low  <- fe_meta_by_exposure(rr, "Low risk")

## --- Plot helper (mimics your styling) ---
plot_meta <- function(df, ylim_pair = c(0.25, 3), colors = c("black")) {
  # ensure the color mapping exists; if only one level, one color is fine
  levs <- unique(df$`Machine Learning Prediction`)
  if (length(colors) < length(levs)) colors <- rep(colors, length.out = length(levs))
  
  g <- ggplot(
    df,
    aes(
      x = rowid,
      y = boxOdds,
      colour = `Machine Learning Prediction`,
      group = `Machine Learning Prediction`
    )
  ) +
    geom_errorbar(aes(ymin = boxCILow, ymax = boxCIHigh), colour = "black", width = 0.1) +
    geom_point(size = 4) +
    labs(y = "HR (95% CI)", x = "") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.5) +
    coord_flip(ylim = ylim_pair, xlim = c(1, nrow(df))) +
    scale_color_manual(values = colors) +
    theme_bw() + theme(legend.position = "none") +
    scale_x_continuous(
      breaks = df$rowid,
      labels = forcats::fct_inorder(df$Stimulus2),
      sec.axis = sec_axis(
        ~ .,
        breaks = df$rowid,
        labels = df$Lab
      )
    ) +
    theme(axis.ticks.y = element_blank()) +
    scale_y_continuous(trans = "log10") +
    ggtitle("")
  g
}

## --- Make the plots ---
# High risk
g_high <- plot_meta(df_high, ylim_pair = c(0.18, 9), colors = c("darkred"))
print(g_high)

# Low risk 
g_low  <- plot_meta(df_low,  ylim_pair = c(0.25, 3), colors = c("#56B4E9"))
print(g_low)
