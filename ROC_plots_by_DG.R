## ===========================
## ROC curves stratified by relapse diagnosis
## Author: M.D., Ph.D., Johannes Lieslehto
## email: johannes.lieslehto@ki.se


## ---- Packages ----
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  dplyr, pROC, ggplot2,
  forcats)

## ---- Helper Functions: ----
# Make a positive-vs-all-negative subset for a set of ICD codes (exact code matches)
subset_pos_vs_neg <- function(dat, exact_codes = NULL) {
  dat <- dat %>% filter(!is.na(outcome), !is.na(rel_mod_cal))
  pos <- dat %>% filter(outcome == "1")
  if (!is.null(exact_codes)) {
    # exact matches only (e.g., "F314", not "F3141")
    rx <- paste0("^(", paste(exact_codes, collapse = "|"), ")$")
    pos <- pos %>% filter(grepl(rx, relapsedg))
  }
  neg <- dat %>% filter(outcome == "0")
  bind_rows(pos, neg)
}

# Safe ROC (returns NULL if not enough classes)
safe_roc <- function(df) {
  df <- df %>% filter(!is.na(outcome), !is.na(rel_mod_cal))
  if (n_distinct(df$outcome) < 2) return(NULL)
  pROC::roc(df$outcome, df$rel_mod_cal, levels = c("0", "1"), quiet = TRUE)
}

# Dynamic label: Name (AUC=0.XX, 95%CI=0.XX–0.XX)
fmt_label <- function(name, rocobj) {
  ci <- pROC::ci.auc(rocobj)
  sprintf("%s (AUROC=%.2f, 95%%CI=%.2f–%.2f)",
          name, as.numeric(pROC::auc(rocobj)), ci[1], ci[3])
}

# Build all ROC objects for one dataset and plot
plot_relapse_rocs <- function(dat, main_title, pdf_file,
                              y_limits = c(0.04, 0.97), x_limits = c(0.03, 0.97)) {
  
  N <- nrow(dat)
  
  ## All-cause (all positives vs all negatives)
  all_tt   <- dat %>% filter(!is.na(outcome), !is.na(rel_mod_cal))
  roc_all  <- safe_roc(all_tt)
  
  ## Cause-specific groups (exact codes)
  # Depressions (F313–F315), Mixed (F316)
  dep_tt   <- subset_pos_vs_neg(dat, c("F313","F314","F315"))
  mix_tt   <- subset_pos_vs_neg(dat, c("F316"))
  # Manias: All Mania codes, and then split by psychosis
  man_all  <- subset_pos_vs_neg(dat, c("F30","F300", "F310","F311","F312", "F3120","F31","F301","F302", "F3021", "F3121", "F308", "F309"))
  man_np   <- subset_pos_vs_neg(dat, c("F311","F301"))  # non-psychotic
  man_ps   <- subset_pos_vs_neg(dat, c("F312","F302", "F3120", "F3121", "F3020", "F3021"))  # psychotic
  # Non-psychotic / psychotic depression
  dep_np   <- subset_pos_vs_neg(dat, c("F314")) # Severe depression without psychosis
  dep_ps   <- subset_pos_vs_neg(dat, c("F315")) # Severe depression with psychosis
  
  # ROC objects
  roc_dep    <- safe_roc(dep_tt)
  roc_man    <- safe_roc(man_all)
  roc_dep_np <- safe_roc(dep_np)
  roc_dep_ps <- safe_roc(dep_ps)
  roc_mix    <- safe_roc(mix_tt)
  roc_man_np <- safe_roc(man_np)
  roc_man_ps <- safe_roc(man_ps)
  
  # Assemble named ROC list with dynamic labels (skip NULLs)
  roc_list <- list()
  add_curve <- function(obj, name) {
    if (!is.null(obj)) {
      nm <- fmt_label(name, obj)
      roc_list[[nm]] <<- obj
    }
  }
  
  add_curve(roc_all,    "All Bipolar Relapses")
  add_curve(roc_dep,    "All Bipolar Depressions")
  add_curve(roc_man,    "All Bipolar Manias")
  add_curve(roc_dep_np, "Non-Psychotic Depression")
  add_curve(roc_dep_ps, "Psychotic Depression")
  add_curve(roc_mix,    "Mixed Episode")
  add_curve(roc_man_np, "Non-Psychotic Mania")
  add_curve(roc_man_ps, "Psychotic Mania")
  
  # Colors (repeat/trim to match the number of curves)
  cols <- c("black", "#DF8F44FF", "#00A1D5FF", "#B24745FF",
            "#79AF97FF", "#6A6599FF", "#80796BFF", "#374E55FF")
  cols <- rep(cols, length.out = length(roc_list))
  
  # Plot
  g <- pROC::ggroc(roc_list, size = 0.7, legacy.axes = TRUE) +
    ggplot2::labs(
      title = sprintf("%s (N=%s)", main_title, format(N, big.mark = " ")),
      x = "1 - Specificity",
      y = "Sensitivity"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = 1, y = 0, yend = 1),
      color = "darkgrey", linetype = "dashed"
    ) +
    ggplot2::coord_cartesian(ylim = y_limits, xlim = x_limits) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = c(0.72, 0.25),
      legend.title    = ggplot2::element_blank(),
      legend.text     = ggplot2::element_text(size = 6.1)
    )
  
  # Save + print
  grDevices::pdf(pdf_file, width = 6.3, height = 6)
  print(g)
  grDevices::dev.off()
  print(g)
}

## ---- Run for Swedish internal & Finnish external ----

# Internal Validation: Swedish data
plot_relapse_rocs(
  dat       = test_dg,  # or sweden internal frame
  main_title = "Internal Validation: Swedish Data",
  pdf_file  = "ROC_kayra_se_bipo_inval.pdf",
  y_limits  = c(0.04, 0.97),
  x_limits  = c(0.03, 0.97)
)

# External Validation: Finnish data
plot_relapse_rocs(
  dat       = fin_dg,
  main_title = "External Validation: Finnish Data",
  pdf_file  = "ROC_kayra_fi_bipo_val.pdf",
  y_limits  = c(0.04, 0.97),
  x_limits  = c(0.03, 0.97)
)
