# SmoothTensorDecompositionABPM

This is a repository for the implementation of the paper "Smooth Tensor Decomposition with Application to Ambulatory Blood Pressure Monitoring Data". The implementation depends on the [SmoothHOOI](https://github.com/IrinaStatsLab/SmoothHOOI) R package.

Data used in **HYPNOS Application** are confidential. Synthetically generated ABPM data that mimics the characteristics of data analyzed in the paper is presented in **Synthetic Example**.

## HYPNOS Application

This folder includes code for reproducing the results in Section 4 of the paper.

- `1a_Data Preprocessing.R`: data preprocessing, including low-quality data detection and removal
- `1b_Tensor Generation.R`: organization of original data into a tensor structure
- `2a_Autocorrelation.R`: check residual autocorrelation in the ABPM data
- `2b_Hyperparameter Tuning.R`: hyperparameter tuning, including rank reduction for parsimony
- `3_Algorithm Run.R`: implementation of SmoothHOOI algorithm, with optimal hyperparameter applied
- `4a_Result Visualization.R`: visualization of temporal components and estimated curves
- `4b_Chronotype Analysis.R`: validation of the interpretation of the third temporal component (plot of g score vs sleep times)
- `4c_Regression.R`: regression analysis
- `4d_Regression Interpretation_Effect_size.R`: visualization of effect sizes of all the variables
- `4e_Regression Interpretation_Profile.R`: estimation of DBP, SBP, and HR profiles for different groups of subjects
- `4f_CP_Analysis.R`: instability of CP decomposition in ABPM data
- `4g_MFPCA_Analysis.R`: MFPCA analysis for ABPM data 

## Simulation Studies

This folder includes code for reproducing the results in Section 3 of the paper. 

In `Case 1-Fixed ranks` and `Case 2-Flexible ranks` folders, the following abbreviations were used to name the files:
- `missing_rate`: random missingness
- `missing_struc`: structured missingness
- `noise_level`: noise level
- `p`: sample size

`sim_analysis.R` includes the code for plotting the figures related to simulation studies.

## Synthetic Example

This folder presents a synthetic example that shows the workflow of this study and the usage of the [SmoothHOOI](https://github.com/IrinaStatsLab/SmoothHOOI) R package.

- `synthetic_raw.Rda`: L, R, mean of G scores, covariance of G scores, and empirical residuals generated from HYPNOS Application, used to generate `synthetic_data.Rda`
- `synthetic_data.Rda`: synthetic data
- `Synthetic Example.Rmd`: code for this synthetic example
- `Synthetic-Example.pdf`: output for this synthetic example


