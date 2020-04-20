# Code to implement the method from the paper entitled "Detecting Participant Noncompliance Across Multiple Time Points By Modeling a Longitudinal Biomarker".

## Requires the R programming language with packages "gtools", "mvtnorm", "matrixStats", and "pROC".

## There are four files:
### 1) functions.R, an R file that defines a number of objects and functions required to implement the method. The true values of the parameters are defined here.
### 2) parameter_estimation_no_cens.Rmd, an RMD file that gives an example of how to estimate the parameters without censoring.
### 3) parameter_estimation_cens.Rmd, an RMD file that gives an example of how to estimate the parameters with a left-censored biomarker.
### 4) auc_calculation, an RMD file that uses a given set of parameter values to calculate the AUC values for the three compliance probabilities. A comparison with the method from Boatman et. al is included.

## Each file is defaulted to generate data as Compound Symmetry. Code in each file can be uncommented to generate data as AR(1). The models fit assume Compound Symmetry.

## Written by Ross Peterson, PhD Candidate in Biostatistics at the University of Minnesota. Email: pet00180@umn.edu