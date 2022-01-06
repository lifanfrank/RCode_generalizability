This folder contains the simulation code to generate Table 3, Web Table 9 to 11. 

1_EstimatedPS: This folder contains the simulation code for generalizability estimators with the estimated propensity scores;

2_TruePS: This folder contains the simulation code for generalizability estimators with the true propensity scores;

In both folders, there are four sub-folders: 

1_gamma1_alpha1: moderate selection effect in trial participation and moderate effect modification;
2_gamma2_alpha1: strong selection effect in trial participation and moderate effect modification;
3_gamma1_alpha2: moderate selection effect in trial participation and strong effect modification;
4_gamma2_alpha2: strong selection effect in trial participation and strong effect modification;

In each of these sub-folders, the script files are:

1_IPSW_a.R: IPSW1 and IPSW2 when the sampling score model is correctly specified;
2_IPSW_b.R: IPSW1 and IPSW2 when the sampling score model is misspecified;
3_OR_a.R:   REG when the outcome model is correctly specified;
4_OR_b.R:   REG when the outcome model is misspecified;
5_DR_a.R:   DR1 and DR2 when the sampling score model and outcome model are both correctly specified;
6_DR_ps.R:  DR1 and DR2 when the sampling score model is correctly specified but the outcome model is misspecified;
7_DR_po.R:  DR1 and DR2 when the sampling score model is misspecified but the outcome model is correctly specified;
8_DR_d.R:   DR1 and DR2 when the sampling score model and outcome model are both misspecified;


The "functions" folder: 

This folder contains the specific functions that are called by the simulation scripts. The files are:

DR_est.R:  DR1 and DR2 estimators with the estimated propensity scores;
DR_true.R: DR1 and DR2 estimators with the true propensity scores;
IPSW_est.R: IPSW1 and IPSW2 estimators with the estimated propensity scores;
IPSW_true.R: IPSW1 and IPSW2 estimators with the true propensity scores;
OR_est.R:  REG estimators with the estimated propensity scores;
OR_true.R: REG estimators with the true propensity scores;
simdata.R: data generation process.




