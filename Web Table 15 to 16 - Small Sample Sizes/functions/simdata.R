##############################################
# function to simulate data
# Input: 
# gamma: coefficient in sampling score model
# alpha: coefficient of effect modifiers in 
#        potential outcome models
# m: size of the cohort
##############################################
simdata<-function(gamma,alpha,m){
  N=10^6 #set size of target
  tid<-1:N
  `%ni%` = Negate(`%in%`) #not in function
  
  p1<-0.4 #prev of binary covariate
  Z1<-rbinom(N, 1, p1)
  Z2<-rnorm(N)
  Z3<-Z1*Z2
  ps<-plogis(gamma[1]+gamma[2]*Z1+gamma[3]*Z2+gamma[4]*Z3)
  target <- data.frame(cbind(tid,Z1,Z2,Z3,ps)) 
  
  #Trial as a biased sample from target 
  target$S<-rbinom(N,1,target$ps)
  trial<-target[which(target$S==1),]
  n<-dim(trial)[1]
  
  #Cohort as SRS without replacement from target 
  #(less trial participants)
  trialtid<-trial$tid
  ntarget<-target[target$tid %ni% trialtid, ]
  nt<-(N-n)
  cohort<-target[sample(ntarget$tid,m,replace=F), ]
  cohort$S<-0
  
  #combine trial and cohort (S,Z) to estimate propensity scores
  #observed data (D=1)
  both<-rbind(cohort,trial) 
  
  #other variables for trial
  pt=0.5 #prob of treatment
  trial$X<-rbinom(n, 1, pt)
  trial$epsilon<-rnorm(n)
  
  #We need Y to depend on the covariates + treatment
  #outcome Y (EM by Z1, Z2, Z3)
  nu1<-nu2<-nu3<- -1
  xi<-2
  trial$Y<-nu1*trial$Z1 + nu2*trial$Z2 + nu3*trial$Z3 + 
    xi*trial$X + alpha[1]*trial$Z1*trial$X + 
    alpha[2]*trial$Z2*trial$X + alpha[3]*trial$Z3*trial$X + trial$epsilon 
  
  # Return objects
  return(list(N=N, trial=trial, cohort=cohort, both=both))
}
