###################################################################
# example use of function with simulated data
# need to change file directory to read in the functions
###################################################################

# set true coefficients
source("...\simdata.R")       ### change file directory here ###
gamma<-c(-7.374, 0.6,0.6)
alpha<-c(1,1)
# Delta<-2.4

set.seed(0601)
data<-simdata(gamma,alpha)    # generate data
N<-data$N                     # populatio size
trial<-data$trial             # trial portion of the population
cohort<-data$cohort           # simulated representative cohort
both<-data$both               # combined

Y<-trial$Y                    # (trial) outcome
S<-both$S                     # (both trial and cohort) participation indicator
Z<-cbind(1,both$Z1,both$Z2)   # (both trial and cohort) design matrix (including intercept)
X<-trial$X                    # (trial) treatment indicator

# analysis with true treatment propensity scores
source("...\IPSW_truePS.R")    ### change file directory here ###
source("...\OR_truePS.R")     ### change file directory here ###
source("...\DR_truePS.R")     ### change file directory here ###

IPSW1(Y=Y,S=S,Z=Z,X=X,N=N)
IPSW2(Y=Y,S=S,Z=Z,X=X,N=N)
OR1(Y=Y,S=S,Z=Z,X=X,N=N)
DR1(Y=Y,S=S,Zp=Z,Zo=Z,X=X,N=N)
DR2(Y=Y,S=S,Zp=Z,Zo=Z,X=X,N=N)

# analysis with estimated treatment propensity scores
source("...\IPSW_estPS.R")    ### change file directory here ###
source("...\OR_estPS.R")     ### change file directory here ###
source("...\DR_estPS.R")     ### change file directory here ###

IPSW1(Y=Y,S=S,Zp=Z,Ze=Z,X=X,N=N)
IPSW2(Y=Y,S=S,Zp=Z,Ze=Z,X=X,N=N)
OR1(Y=Y,S=S,Ze=Z,Zo=Z,X=X,N=N)
DR1(Y=Y,S=S,Zp=Z,Ze=Z,Zo=Z,X=X,N=N)
DR2(Y=Y,S=S,Zp=Z,Ze=Z,Zo=Z,X=X,N=N)




