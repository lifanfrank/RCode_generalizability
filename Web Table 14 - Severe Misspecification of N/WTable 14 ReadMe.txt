5_gamma2_alpha2_0.1million: 

This folder contains simulation code to generate Web Table 14, when the target population size is misspecified to be 0.1 million. 

5_gamma2_alpha2_0.05million:

This folder contains simulation code to generate Web Table 14, when the target population size is misspecified to be 0.05 million. 

In both folders, the files are

1_IPSW_a.R: IPSW1 and IPSW2 when the sampling score model is correctly specified;
2_IPSW_b.R: IPSW1 and IPSW2 when the sampling score model is misspecified;
3_OR_a.R:   REG when the outcome model is correctly specified;
4_OR_b.R:   REG when the outcome model is misspecified;
5_DR_a.R:   DR1 and DR2 when the sampling score model and outcome model are both correctly specified;
6_DR_ps.R:  DR1 and DR2 when the sampling score model is correctly specified but the outcome model is misspecified;
7_DR_po.R:  DR1 and DR2 when the sampling score model is misspecified but the outcome model is correctly specified;
8_DR_d.R:   DR1 and DR2 when the sampling score model and outcome model are both misspecified;
