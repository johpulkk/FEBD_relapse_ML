# FEBD_relapse_ML  
📊 Machine Learning to Predict Bipolar Relapse in First-Episode Bipolar Disorder

[![R](https://img.shields.io/badge/R-4.3+-blue)](https://www.r-project.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

---

<a href="https://johannes-lieslehto.shinyapps.io/biporacle/">
  <img src="logo.jpg" alt="BipOracle" width="200"/>
</a>

## Overview
This repository provides R scripts for machine learning analyses to predict the risk of bipolar relapse 
in **first-episode bipolar disorder (FEBD)**.  

The analyses include:
- Development and validation of ML models
- Comparison across algorithms
- Fairness analyses
- Time-to-event methods
- Pharmacoepidemiological meta-analyses

---

## Repository Contents

| Script                      | Purpose |
|------------------------------|---------|
| **Dev_Val_Cal.R**            | XGBoost ML model development, internal & external validation, calibration, and decision curve analyses |
| **Compare_ML_Algorithms.R**  | Comparison of XGBoost vs. other ML algorithms in predicting relapse |
| **ROC_plots_by_DG.R**        | Comparison of discrimination performance (AUROC), stratified by cause of relapse |
| **Fairness_analyses.R**      | Fairness analyses: discrimination and calibration across subgroups |
| **Time2event.R**             | Time-to-event analyses using the full follow-up data |
| **MetaAnalysis_PharmaEpi.R** | Pharmacoepidemiological meta-analyses stratified by predicted relapse risk, with plotting |

