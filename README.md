# R Code for "Handling incomplete outcomes and covariates in cluster-randomized trials: doubly-robust estimation, efficiency considerations, and sensitivity analysis"

## Paper

The technical report is avaialble at [https://arxiv.org/abs/2401.11278](https://arxiv.org/abs/2312.03967).

## Data

### Abstract 

To demonstrate our theoretical results, we performed simulation studies and a data application on the electronic health records at the University of Michigan Health System

### Availability 
The EHR data are subject to the privacy rule at the University of Michigan Precision Health at https://precisionhealth.umich.edu/our-research/documents-for-researchers/.

## Code

### Abstract
The R code here can be used to reproduce all of the simulations and data applications.

### Description 

 - `TND-multiple-reasons.R` contains the code of simulations for no effect modification from X or H
 - `TND-multiple-reasons2.R` contains the code of simulations for effect modification from X and H
 - `data-application/data-preprocessing.R` contains the code for data preprocesssing
 - `data-application/ICD_analysis.R` contains the code for data application

### Reproducibility 
Numerical results in the main text and appendix can be reproduced using the R scripts.
