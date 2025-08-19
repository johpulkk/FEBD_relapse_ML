## ===========================
## Model development, internal & external 
## validation, and calibration analyses:
## Author: M.D., Ph.D., Johannes Lieslehto
## email: johannes.lieslehto@ki.se


## ---- Required Packages ----
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  mlr, ParamHelpers, xgboost,
  pROC, caret,
  dplyr, data.table,
  pmsampsize, SHAPforxgboost,
  parallel, parallelMap, rms, dcurves
)

set.seed(42)

## ---- Minimum sample size calculation ----
# Binary outcome: target C-statistic 0.7, max 15 parameters, 
# prevalence 6%, shrinkage 0.9
pmsampsize(type = "b", cstatistic = 0.7, parameters = 15, 
           prevalence = 0.06, shrinkage = 0.9)


## ---- Geographical split by county (Län in Swedish) ----
lan_levels <- sort(unique(sweden_df$Lan))
## Sample 10 counties for development
train_lan <- sample(lan_levels, 
                    size = min(10, length(lan_levels)), replace = FALSE)

dev_df  <- dplyr::filter(sweden_df, Lan %in% train_lan)
ival_df <- dplyr::filter(sweden_df, !Lan %in% train_lan)

## Drop the splitting variable (Lan) AFTER splitting
dev_df  <- dplyr::select(dev_df, -Lan)
ival_df <- dplyr::select(ival_df, -Lan)

## ---- Preprocessing ----
## 1) Predictor-wise: drop columns with >20% missingness (keep outcome!)
miss_prop <- function(df) sapply(dplyr::select(df, -outcome), 
                                 function(x) mean(is.na(x)))
bad_cols <- names(which(miss_prop(dev_df) > 0.20))
if (length(bad_cols)) {
  message("Dropping predictors with >20% missingness: ", 
          paste(bad_cols, collapse = ", "))
  dev_df  <- dplyr::select(dev_df,  -dplyr::all_of(bad_cols))
  ival_df <- dplyr::select(ival_df, -dplyr::any_of(bad_cols))
  finland_df <- dplyr::select(finland_df, -dplyr::any_of(bad_cols))
}


## 2) Row-wise: remove individuals with >40% missingness across ALL columns (including outcome)
drop_rows_by_na <- function(df, max_prop_na = 0.40) df[rowMeans(is.na(df)) <= max_prop_na, , drop = FALSE]
dev_df  <- drop_rows_by_na(dev_df,  0.40)
ival_df <- drop_rows_by_na(ival_df, 0.40)
finland_df <- drop_rows_by_na(finland_df, 0.40)

## 3) Ensure the outcome is present
dev_df$outcome  <- factor(dev_df$outcome,  levels = c("0","1"))
ival_df$outcome <- factor(ival_df$outcome, levels = c("0","1"))
finland_df$outcome <- factor(finland_df$outcome, levels = c("0","1"))


## ---- mlr tasks ----
train_task <- mlr::makeClassifTask(data = dev_df,  target = "outcome", positive = "1")
test_task  <- mlr::makeClassifTask(data = ival_df, target = "outcome", positive = "1")

## ---- XGBoost modeling ----
## Set a fixed scale_pos_weight to mitigate class imbalance:
tbl <- table(dev_df$outcome)
pos_weight <- as.numeric(tbl[names(tbl) == "1"]) / as.numeric(tbl[names(tbl) == "0"]) 
if (!is.finite(pos_weight)) pos_weight <- 1

xgb_learner <- mlr::makeLearner(
  "classif.xgboost",
  predict.type = "prob",
  par.vals = list(
    objective        = "binary:logistic",
    eval_metric      = "logloss",
    scale_pos_weight = pos_weight
  )
)

## Hyperparameter search space:
xgb_param_space <- ParamHelpers::makeParamSet(
  ParamHelpers::makeDiscreteParam("eta",              values = c(2^-8, 2^-6, 2^-4, 2^-2)),
  ParamHelpers::makeDiscreteParam("gamma",            values = c(2^-16, 2^-6, 2^2)),
  ParamHelpers::makeDiscreteParam("subsample",        values = c(0.5, 0.8)),
  ParamHelpers::makeDiscreteParam("max_depth",        values = c(3, 5, 7)),
  ParamHelpers::makeDiscreteParam("nrounds",          values = c(100, 300, 500)),
  ParamHelpers::makeDiscreteParam("min_child_weight", values = c(2^-16, 2^-6, 2^4, 2^8)),
  ParamHelpers::makeDiscreteParam("colsample_bytree", values = c(0.5, 0.8))
)

## Nested resampling
inner_desc <- mlr::makeResampleDesc("CV", iters = 3) # Inner cross-validation
outer_desc <- mlr::makeResampleDesc("CV", iters = 10) # Outer cross-validation

xgb_wrapper <- mlr::makeTuneWrapper(
  learner     = xgb_learner,
  resampling  = inner_desc,
  par.set     = xgb_param_space,
  control     = mlr::makeTuneControlGrid(),
  measures    = list(mlr::bac)  ## balanced accuracy
)

## ---- Run nested CV (parallelized) ----
parallelMap::parallelStartSocket(cpus = parallel::detectCores())
cv_res <- mlr::resample(
  learner     = xgb_wrapper,
  task        = train_task,
  resampling  = outer_desc,
  measures    = list(mlr::bac, mlr::auc),
  models      = TRUE,
  show.info   = TRUE
)
parallelMap::parallelStop()

print(cv_res$aggr)  ## mean BAC/AUC across outer folds

## ---- Variable importance across outer models ----
## Unwrap to the underlying xgboost model and average 'Gain' across CV folds
unwrap_to_xgb <- function(wrapped_model) mlr::getLearnerModel(wrapped_model, 
                                                              more.unwrap = TRUE)

imp_list <- lapply(cv_res$models, function(m) {
  xm <- try(unwrap_to_xgb(m), silent = TRUE)
  if (inherits(xm, "try-error")) return(NULL)
  xgboost::xgb.importance(model = xm)
})
imp_list <- Filter(Negate(is.null), imp_list)

if (length(imp_list) == 0) stop("No xgboost models found for importance extraction.")
imp_df <- data.table::rbindlist(imp_list, fill = TRUE)
imp_summary <- imp_df %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(mean_gain = mean(Gain, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(mean_gain))

top_k <- 15 # Selecting the top 15 features
top_features <- head(imp_summary$Feature, top_k)
print(head(imp_summary, top_k))

## ---- Sequential Forward Selection (SFS) on the top features ----
dev_top  <- dplyr::select(dev_df,  dplyr::any_of(c(top_features, "outcome")))
ival_top <- dplyr::select(ival_df, dplyr::any_of(c(top_features, "outcome")))
train_top_task <- mlr::makeClassifTask(data = dev_top, 
                                       target = "outcome", positive = "1")

feat_ctrl <- mlr::makeFeatSelControlSequential(method = "sfs", alpha = 5e-4)
parallelMap::parallelStartSocket(cpus = parallel::detectCores())
feat_sel <- mlr::selectFeatures(
  learner     = xgb_wrapper,
  task        = train_top_task,
  resampling  = inner_desc,
  control     = feat_ctrl,
  measures    = list(mlr::bac),
  show.info   = TRUE
)
parallelMap::parallelStop()

sel_feats <- feat_sel$x
message("Selected features (SFS): ", paste(sel_feats, collapse = ", "))

## ---- Final tuning on selected features + training ----
dev_final  <- dplyr::select(dev_df,  dplyr::all_of(c(sel_feats, "outcome")))
ival_final <- dplyr::select(ival_df, dplyr::all_of(c(sel_feats, "outcome")))
final_task <- mlr::makeClassifTask(data = dev_final, 
                                   target = "outcome", positive = "1")

parallelMap::parallelStartSocket(cpus = parallel::detectCores())
final_tune <- mlr::tuneParams(
  learner     = xgb_learner,
  task        = final_task,
  resampling  = inner_desc,
  par.set     = xgb_param_space,
  control     = mlr::makeTuneControlGrid(),
  measures    = list(mlr::bac)
)
parallelMap::parallelStop()

final_learner <- mlr::setHyperPars(xgb_learner, par.vals = final_tune$x)
final_model   <- mlr::train(final_learner, final_task)

## --- SHAP values & plots (using SHAPforxgboost) ---
# Example plot using the validation sample

# 1) Get the underlying xgboost booster from the mlr model
xgb_booster <- mlr::getLearnerModel(final_model, more.unwrap = TRUE)

# 2) Choose the data to explain (here: internal validation set)
#    Use exactly the features the model was trained on (no outcome column)
feat_names <- sel_feats     
X_val <- data.matrix(dplyr::select(ival_final, dplyr::all_of(feat_names)))

# 3) Compute SHAP values (long format)
shap_long <- SHAPforxgboost::shap.prep(
  xgb_model = xgb_booster,
  X_train   = X_val
)

# 4) Summary beeswarm plot
SHAPforxgboost::shap.plot.summary(shap_long, scientific = TRUE)


## ---- Evaluation on development + held-out internal validation ----
pred_dev  <- mlr::predict(final_model, newdata = dev_final)

pred_ival <- mlr::predict(final_model, newdata = ival_final)

cat("\nDevelopment set performance:\n")
print(mlr::performance(pred_dev,  measures = list(mlr::auc, mlr::bac)))

cat("\nInternal validation performance:\n")
print(mlr::performance(pred_ival, measures = list(mlr::auc, mlr::bac)))

# Plot ROC curve (internal validation sample):
plot(pROC::roc(pred_ival$data$truth, pred_ival$data$prob.1))

## Confusion matrix on the held-out internal validation:
ival_pred_label <- factor(ifelse(pred_ival$data$prob.1 > 0.5, "1", "0"),
                          levels = c("0", "1"))
print(caret::confusionMatrix(ival_pred_label, 
                             pred_ival$data$truth, positive = "1"))

## ---- External validation on Finnish external validation sample ---- 

# A similar approach was used for transdiagnostic validation

fi_final <- dplyr::select(finland_df, dplyr::any_of(c(sel_feats, "outcome")))
pred_fi  <- mlr::predict(final_model, newdata = fi_final)
cat("\nExternal validation (Finland):\n")
print(mlr::performance(pred_fi, measures = list(mlr::auc, mlr::bac)))
## pROC ROC curve / AUC if you want the classic object:
roc_fi <- pROC::roc(pred_fi$data$truth, pred_fi$data$prob.1)

# Plot ROC curve (Finnish external validation sample):
plot(roc_fi)

## Confusion matrix on the Finnish external validation sample:
eval_pred_label <- factor(ifelse(pred_fi$data$prob.1 > 0.5, "1", "0"),
                          levels = c("0", "1"))
print(caret::confusionMatrix(eval_pred_label, 
                             pred_fi$data$truth, positive = "1"))


## ---- Calibration with logistic regression (Platt scaling) ----
## fit on development; apply to others

to01 <- function(y) as.integer(as.character(y) == "1")

# 1) Fit Platt on development
y_dev <- to01(pred_dev$data$truth)
p_dev <- pred_dev$data$prob.1

cal_platt <- glm(y ~ p, family = binomial(),
                 data = data.frame(y = y_dev, p = p_dev))

# 2) Apply to development, internal and external validation samples
p_dev_cal  <- predict(cal_platt, 
                      newdata = data.frame(p = pred_dev$data$prob.1),  
                      type = "response")
p_ival_cal <- predict(cal_platt, 
                      newdata = data.frame(p = pred_ival$data$prob.1), 
                      type = "response")
p_fi_cal <- predict(cal_platt, 
                    newdata = data.frame(p = pred_fi$data$prob.1), 
                    type = "response")


# 3) Check calibration and draw calibration plots
cat("\nDevelopment (Platt-calibrated):\n")
val.prob(p_dev_cal,  y_dev)

cat("\nInternal validation (Platt-calibrated):\n")
val.prob(p_ival_cal, to01(pred_ival$data$truth))

cat("\nExternal Finland (Platt-calibrated):\n")
val.prob(p_fi_cal, to01(pred_fi$data$truth))

## ---- DCA Analysis ----

pred_ival$data$p_ival_cal<-p_ival_cal
pred_fi$data$p_fi_cal<-p_fi_cal

b<-dcurves::dca(to01(pred_ival$data$truth) ~ p_ival_cal, pred_ival, thresholds = seq(0, 0.5, by = 0.002), list(p_ival_cal="XGBoost Model")) %>%
  net_intervention_avoided()
plot(b, smooth = T, show_ggplot_code = T)+
  ggplot2::scale_color_manual(values = c("#00A1D5FF","black",  "#DF8F44FF"))

b<-dcurves::dca(to01(pred_fi$data$truth) ~ p_fi_cal, pred_fi, thresholds = seq(0, 0.5, by = 0.002), list(p_fi_cal="XGBoost Model")) %>%
  net_intervention_avoided()
plot(b, smooth = T, show_ggplot_code = T)+
  ggplot2::scale_color_manual(values = c("#00A1D5FF","black",  "#DF8F44FF"))




