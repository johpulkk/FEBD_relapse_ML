# FEBD_relapse_ML
This repository provides the scripts to generate machine learning (ML) analyses to predict bipolar relapse in first-episode bipolar disorder (FEBD). 

Dev_Val_Cal.R = XGBoost ML model development, internal & external validation, calibration, and decision curve analyses.
Compare_ML_Algorithms.R = Comparison of XGBoost vs. other ML algorithms in predicting bipolar relapse. 
ROC_plots_by_DG.R = Comparison of discrimination performance (by AUROC), stratified by the specific cause of relapse.
Fairness_analyses.R = Fairness analyses: discrimination and calibration differences in different subgroups.
Time2event.R = Time to event analyses using the full follow-up data.
MetaAnalysis_PharmaEpi.R = Pharmacoepidemiological meta-analyses, stratified by predicted relapse risk, with plotting.
