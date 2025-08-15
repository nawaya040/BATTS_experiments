# BATTS_experiments

This repository provides program code to reproduce the numerical experiments in Awaya, Xu, and Ma (2025).

The details are provided in 

Awaya, Xu, and Ma (2025). Two-sample comparison through additive tree models for density ratios. arXiv preprint arXiv:2508.03059.

To run the programs, the R package `BATTS` needs to be installed  using `devtools` as follows:

```
library(devtools)
install_github("nawaya040/BATTS")
```

The following files are included:

1. `experiment_2D.R`: density estimation for multi-variate data sets (Section 4.1)
2. `experiment_multi.R`: evaluation of the variable importance (Section 4.2)
3. `helpers`: programs to simulate data sets and run the other methods

