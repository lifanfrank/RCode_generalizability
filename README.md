# RCode_generalizability
Sample R Code for generalizing trial results to target populations in Li et al. (submitted)

You will need to change the directory to use the example code in script example_analysis.R
This folder includes the following functions:

1. IPSW_truePS.R : implement the IPSW1_tilde and IPSW2_tilde estimators for PATE and the associated sandwich variance estimators, the true treatment propensity score (e_i=1/2) is used;
2. OR_truePS.R : implement the OR_tilde estimator for PATE and the associated sandwich variance estimator, the true treatment propensity score (e_i=1/2) is used;
3. DR_truePS.R : implement the DR1_tilde and DR2_tilde estimators for PATE and the associated sandwich variance estimators, the true treatment propensity score (e_i=1/2) is used;
4. IPSW_estPS.R : implement the IPSW1_hat and IPSW2_hat estimators for PATE and the associated sandwich variance estimators, the estimated treatment propensity score (e_i_hat) is used;
5. OR_estPS.R : implement the OR_hat estimator for PATE and the associated sandwich variance estimator, the estimated treatment propensity score (e_i_hat) is used;
6. DR_estPS.R : implement the DR1_tilde and DR2_tilde estimators for PATE and the associated sandwich variance estimators, the estimated treatment propensity score (e_i_hat) is used;
7. simdata.R : sample R code to generate a simulated data set for illustrative analysis;
8. example_analysis.R : illustrative R code to perform the generalizability analysis using the above functions that implement point and variance estimators developed in Li et al. (submitted).
