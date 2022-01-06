5_gamma2_alpha2_smallTrial:

This folder contains simulation code to generate Web Table 15, when the trial sample size is 200 and the cohort sample size is 4000.

5_gamma2_alpha2_smallCohort: 

This folder contains simulation code to generate Web Table 16, when the trial sample size is 1000 and the cohort sample size is 800.

5_gamma2_alpha2_smallTrialCohort:

This folder contains simulation code to generate Web Table 16, when the trial sample size is 200 and the cohort sample size is 800.

In all three folders, the script files are

1_IPSW_a.R: IPSW1 and IPSW2 when the sampling score model is correctly specified;
2_IPSW_b.R: IPSW1 and IPSW2 when the sampling score model is misspecified;
3_OR_a.R:   REG when the outcome model is correctly specified;
4_OR_b.R:   REG when the outcome model is misspecified;
5_DR_a.R:   DR1 and DR2 when the sampling score model and outcome model are both correctly specified;
6_DR_ps.R:  DR1 and DR2 when the sampling score model is correctly specified but the outcome model is misspecified;
7_DR_po.R:  DR1 and DR2 when the sampling score model is misspecified but the outcome model is correctly specified;
8_DR_d.R:   DR1 and DR2 when the sampling score model and outcome model are both misspecified;


functions: 

This folder contains the specific functions that are called by the simulation scripts. The files are:

DR_est.R:  DR1 and DR2 estimators with the estimated propensity scores;
DR_true.R: DR1 and DR2 estimators with the true propensity scores;
IPSW_est.R: IPSW1 and IPSW2 estimators with the estimated propensity scores;
IPSW_true.R: IPSW1 and IPSW2 estimators with the true propensity scores;
OR_est.R:  REG estimators with the estimated propensity scores;
OR_true.R: REG estimators with the true propensity scores;
simdata.R: data generation process.




