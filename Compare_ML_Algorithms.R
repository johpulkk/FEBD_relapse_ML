## ===========================
## Compare ML Algorithms
## (SVM, RF, Elastic Net) vs XGB
## Author: M.D., Ph.D., Johannes Lieslehto
## email: johannes.lieslehto@ki.se
## ===========================

## ---- Required Packages ----
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  mlr, ParamHelpers, xgboost,
  pROC, caret,
  dplyr, data.table,
  pmsampsize, SHAPforxgboost,
  parallel, parallelMap, rms, dcurves
)


# Hyperparameter tuning CV:
inner_desc <- mlr::makeResampleDesc("CV", iters = 3) 

# 0) Tasks with ALL features (no SFS), same target & positive class
task_all_dev  <- mlr::makeClassifTask(data = dev_df,  target = "outcome", positive = "1")
task_all_ival <- mlr::makeClassifTask(data = ival_df, target = "outcome", positive = "1")

# helper for undersampling rate (majority shrink to match minority)
tbl <- table(dev_df$outcome)
us_rate <- as.numeric(tbl["1"] / tbl["0"])  # 0<rate<=1

# =============== SVM (linear kernel) ===============
svm_base <- mlr::makeLearner("classif.svm", predict.type = "prob")

# SVM-specific preprocessing pipeline

# impute + scale, then undersample majority
svm_pp  <- mlr::makePreprocWrapperCaret(svm_base, ppc.knnImpute = TRUE, ppc.nzv = TRUE,
                                        ppc.zv = TRUE, ppc.center = TRUE, ppc.scale = TRUE)
svm_lrn <- mlr::makeUndersampleWrapper(svm_pp, usw.rate = us_rate)

# Hyperparameter Space:
kernels <- c("linear")

svm_ps <- ParamHelpers::makeParamSet(
  ParamHelpers::makeDiscreteParam("kernel", values = kernels),
  ParamHelpers::makeDiscreteParam("cost", values = 2^seq(-8, 8, by = 2))
)

svm_tw <- mlr::makeTuneWrapper(learner    = svm_lrn,
                               resampling = inner_desc,
                               par.set    = svm_ps,
                               control    = mlr::makeTuneControlGrid(),
                               measures   = list(mlr::bac))

# =============== Random Forest ===============
rf_base <- mlr::makeLearner("classif.randomForest", predict.type = "prob")
# Preprocessing
rf_pp   <- mlr::makePreprocWrapperCaret(rf_base, ppc.knnImpute = TRUE, ppc.nzv = TRUE, ppc.zv = TRUE)
rf_lrn  <- mlr::makeUndersampleWrapper(rf_pp, usw.rate = us_rate)


# Hyperparameter Space: 
rf_ps <- ParamHelpers::makeParamSet(
  ParamHelpers::makeDiscreteParam("ntree", c(500)),
  ParamHelpers::makeDiscreteParam("mtry",     values = c(6,8,10,12)),
  ParamHelpers::makeDiscreteParam("nodesize", values = c(1, 3, 5)),
  ParamHelpers::makeDiscreteParam("maxnodes", values = c(5,10,15))
)

rf_tw <- mlr::makeTuneWrapper(learner    = rf_lrn,
                              resampling = inner_desc,
                              par.set    = rf_ps,
                              control    = mlr::makeTuneControlGrid(),
                              measures   = list(mlr::bac))

# =============== Elastic Net (glmnet) ===============
# cvglmnet does its own inner CV for lambda; we tune alpha and the "s" to use
enet_base <- mlr::makeLearner("classif.cvglmnet",id = "elastic", predict.type = "prob")
# Preprocessing
enet_lrn  <- mlr::makePreprocWrapperCaret(enet_base, ppc.knnImpute = TRUE, ppc.nzv = TRUE,
                                          ppc.zv = TRUE, ppc.center = TRUE, ppc.scale = TRUE)
enet_lrn  <- mlr::makeUndersampleWrapper(enet_lrn, usw.rate = us_rate)

# Hyperparameter Space:
enet_ps <- ParamHelpers::makeParamSet(
  ParamHelpers::makeDiscreteParam("alpha", values = c(0, 0.2, 0.4,0.5, 0.6, 0.8, 1)),   
  ParamHelpers::makeDiscreteParam("s",     values = c(0, 0.2, 0.5, 0.8, 1, 2, 5, 10))
)

enet_tw <- mlr::makeTuneWrapper(learner    = enet_lrn,
                                resampling = inner_desc,
                                par.set    = enet_ps,
                                control    = mlr::makeTuneControlGrid(),
                                measures   = list(mlr::bac))

# =============== Train tuned models on dev, predict on ival ===============
parallelMap::parallelStartSocket(cpus = parallel::detectCores())

svm_model  <- mlr::train(svm_tw,  task_all_dev)
rf_model   <- mlr::train(rf_tw,   task_all_dev)
enet_model <- mlr::train(enet_tw, task_all_dev)

pred_svm   <- predict(svm_model,  task = task_all_ival)
pred_rf    <- predict(rf_model,   task = task_all_ival)
pred_enet  <- predict(enet_model, task = task_all_ival)

parallelMap::parallelStop()

cat("\nAUC (internal validation):\n")
print(list(
  XGB  = mlr::performance(pred_ival, measures = mlr::auc),
  SVM  = mlr::performance(pred_svm,  measures = mlr::auc),
  RF   = mlr::performance(pred_rf,   measures = mlr::auc),
  ENet = mlr::performance(pred_enet, measures = mlr::auc)
))

# =============== Paired DeLong tests vs XGBoost ===============
xgb_roc  <- pROC::roc(pred_ival$data$truth, pred_ival$data$prob.1)
svm_roc  <- pROC::roc(pred_svm$data$truth,  pred_svm$data$prob.1)
rf_roc   <- pROC::roc(pred_rf$data$truth,   pred_rf$data$prob.1)
enet_roc <- pROC::roc(pred_enet$data$truth, pred_enet$data$prob.1)

cat("\nDeLong tests (paired) vs XGB:\n")
print(pROC::roc.test(svm_roc,  xgb_roc, paired = TRUE))
print(pROC::roc.test(rf_roc,   xgb_roc, paired = TRUE))
print(pROC::roc.test(enet_roc, xgb_roc, paired = TRUE))

# =============== Combined ROC plot ===============
# dynamic labels with AUC (mean CI)
fmt_lab <- function(name, rocobj) {
  ci <- pROC::ci.auc(rocobj)
  sprintf("%s (AUC=%.2f, 95%%CI=%.2f–%.2f)", name, as.numeric(pROC::auc(rocobj)), ci[1], ci[3])
}

lbl_xgb  <- fmt_lab("XGBoost",       xgb_roc)
lbl_svm  <- fmt_lab("SVM",  svm_roc)
lbl_rf   <- fmt_lab("Random Forest", rf_roc)
lbl_net  <- fmt_lab("Elastic Net",   enet_roc)

g1<-ggroc(list(`XGBoost with SFS (AUROC=0.71, 95%CI=0.69-0.72)`=xgb_roc,
               `Elastic Net  (AUROC=0.67, 95%CI=0.66-0.69)`=svm_roc,
               `SVM  (AUROC=0.68, 95%CI=0.67-0.70)`=rf_roc,
               `Random Forest  (AUROC=0.68, 95%CI=0.67-0.69)`=enet_roc), 
          size=0.7, legacy.axes = T) +ggtitle("Comparison of the ML Algorithms' Performances") 
g1+ ylab("Sensitivity") + xlab("1-Specificity") + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")+
  theme_bw() +coord_cartesian(ylim = c(0.04, 0.97), xlim = c(0.03, 0.97)) +
  scale_color_manual(values = c("#00A1D5FF", "#B24745FF", "#374E55FF", "#DF8F44FF"))+theme_bw()+
  theme(legend.position = c(0.69, 0.2))+
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=6.1))

