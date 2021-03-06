##############################################
# IPSW estimators for generalizability
# IPSW1 and IPSW2
# true treatment PS is used
##############################################
# Input: 
# Y: (trial) outcome
# S: (both trial and cohort) participation indicator
# Z: (both trial and cohort) design matrix (including intercept) for sampling score
# X: (trial) treatment indicator
# N: target population size
 
IPSW1=function(Y,S,Z,X,N){
  
  # derive weights^-1 for sampling score model
  n <- sum(S); nt <- (N-n); m <- sum(1-S)
  pw <- rep(1, length(S))
  pw[S==0] <- m/nt
  
  #estimate selection propensity scores using weighted logistic regression
  mylogit <- glm(S ~  -1+Z , family = "binomial",weights=pw^(-1))
  gamma.h <- as.numeric(coef(mylogit))
  
  #selection propensity score w_i
  p <- plogis(c(Z%*%gamma.h))
  p1 <- p[S==1]
  
  # true propensity scores
  e <- rep(0.5, length(X))
  
  # fit outcome models
  Y1 <- Y[X==1]
  Y0 <- Y[X==0]
  Z1 <- Z[S==1,]      # in the trial
  Z11 <- Z1[X==1,]    # in the trial and treated
  Z10 <- Z1[X==0,]    # in the trial but untreated
  
  # point estimation
  mu1.h<-sum((X*Y)/(p1*e))/N
  mu0.h<-sum(((1-X)*Y)/(p1*(1-e)))/N
  delta<-mu1.h-mu0.h
  
  # get the dimensions right first
  pp <- ncol(model.matrix(mylogit))
  
  # functions to compute estimating equations
  # here we output a (2+pp) x N matrix
  phi<-function(theta){
    mu1<-theta[1]
    mu0<-theta[2]
    gamma<-theta[3:(2+pp)]
    
    # sampling scores
    p<-plogis(c(Z%*%gamma))
    p1<-p[S==1]
    
    # components
    f1<-f2<-rep(0,n+m)
    f1[S==1]<-(X*Y)/(p1*e)
    f2[S==1]<-((1-X)*Y)/(p1*(1-e))
    f1t<-c(f1,rep(0,N-(n+m)))-rep(mu1,N)
    f2t<-c(f2,rep(0,N-(n+m)))-rep(mu0,N)
    
    f4<-Z*(S-p)*(pw)^(-1)
    f4t<-rbind(f4, matrix(0,N-(n+m),pp))
    
    f <- rbind(f1t,f2t,t(f4t))
    return(f)
  }
  
  # score operator
  # mphi<-function(theta){
  #   as.numeric(rowMeans(phi(theta)))
  # }
  
  # covariance operator
  Omega<-function(theta){
    phis<-phi(theta)
    return(tcrossprod(phis)/N)
  }
  
  # differentiation operator
  Gradient<-function(theta){
    mu1<-theta[1]
    mu0<-theta[2]
    gamma<-theta[3:(2+pp)]
    
    # sampling scores
    p<-plogis(c(Z%*%gamma))
    p1<-p[S==1]
    
    # components
    r1gamma<- colSums(-X*Y*Z1*(1-p1)/(p1*e))/N
    r2gamma<- colSums(-(1-X)*Y*Z1*(1-p1)/(p1*(1-e)))/N
    r1<-c(-1,0,c(r1gamma))
    r2<-c(0,-1,c(r2gamma))
    
    Zstar<-Z*c(sqrt(p*(1-p)*pw^(-1)))
    r4gamma<- -crossprod(Zstar)/N
    r4<- cbind(matrix(0,pp,2), r4gamma)
    
    r <- rbind(r1,r2,r4)
    return(r)
  }
  
  theta.h<-c(mu1.h, mu0.h, gamma.h)
  
  # compute variance
  # Atheta<-jacobian(mphi,theta.h)
  Atheta<-Gradient(theta.h)
  Btheta<-Omega(theta.h)
  invAtheta<-solve(Atheta)
  Vtheta<-invAtheta%*%Btheta%*%t(invAtheta)/N
  
  a<-c(1,-1,rep(0,pp))
  vdelta<-t(a)%*%Vtheta%*%a
  se<-sqrt(vdelta)
  lcl<-delta-qnorm(0.975)*se
  ucl<-delta+qnorm(0.975)*se
  
  # output to screen
  return(list(delta=delta,se=se,lcl=lcl,ucl=ucl))
}

IPSW2=function(Y,S,Z,X,N){
  
  # derive weights^-1 for sampling score model
  n <- sum(S); nt <- (N-n); m <- sum(1-S)
  pw <- rep(1, length(S))
  pw[S==0] <- m/nt
  
  #estimate selection propensity scores using weighted logistic regression
  mylogit <- glm(S ~  -1+Z , family = "binomial",weights=pw^(-1))
  gamma.h <- as.numeric(coef(mylogit))
  
  #selection propensity score w_i
  p <- plogis(c(Z%*%gamma.h))
  p1 <- p[S==1]
  
  # true propensity scores
  e <- rep(0.5, length(X))
  
  # fit outcome models
  Y1 <- Y[X==1]
  Y0 <- Y[X==0]
  Z1 <- Z[S==1,]      # in the trial
  Z11 <- Z1[X==1,]    # in the trial and treated
  Z10 <- Z1[X==0,]    # in the trial but untreated
  
  # point estimation
  mu1.h<-sum((X*Y)/(p1*e))/sum(X/(p1*e))
  mu0.h<-sum(((1-X)*Y)/(p1*(1-e)))/sum((1-X)/(p1*(1-e)))
  delta<-mu1.h-mu0.h
  
  # get the dimensions right first
  pp <- ncol(model.matrix(mylogit))
  
  # functions to compute estimating equations
  # here we output a (2+pp) x N matrix
  phi<-function(theta){
    mu1<-theta[1]
    mu0<-theta[2]
    gamma<-theta[3:(2+pp)]
    
    # sampling scores
    p<-plogis(c(Z%*%gamma))
    p1<-p[S==1]
    
    # components
    f1<-f2<-rep(0,n+m)
    f1[S==1]<-X*(Y-mu1)/(p1*e)
    f2[S==1]<-(1-X)*(Y-mu0)/(p1*(1-e))
    f1t<-c(f1,rep(0,N-(n+m)))
    f2t<-c(f2,rep(0,N-(n+m)))
    
    f4<-Z*(S-p)*(pw)^(-1)
    f4t<-rbind(f4, matrix(0,N-(n+m),pp))
    
    f <- rbind(f1t,f2t,t(f4t))
    return(f)
  }
  
  # score operator
  # mphi<-function(theta){
  #   as.numeric(rowMeans(phi(theta)))
  # }
  
  # covariance operator
  Omega<-function(theta){
    phis<-phi(theta)
    return(tcrossprod(phis)/N)
  }
  
  # differentiation operator
  Gradient<-function(theta){
    mu1<-theta[1]
    mu0<-theta[2]
    gamma<-theta[3:(2+pp)]
    
    # sampling scores
    p<-plogis(c(Z%*%gamma))
    p1<-p[S==1]
    
    # components
    r1gamma<- colSums(-X*(Y-mu1)*Z1*(1-p1)/(p1*e))/N
    r2gamma<- colSums(-(1-X)*(Y-mu0)*Z1*(1-p1)/(p1*(1-e)))/N
    r1<-c(sum(-X/(p1*e))/N, 0, c(r1gamma))
    r2<-c(0, sum(-(1-X)/(p1*(1-e)))/N, c(r2gamma))
    
    Zstar<-Z*c(sqrt(p*(1-p)*pw^(-1)))
    r4gamma<- -crossprod(Zstar)/N
    r4<- cbind(matrix(0,pp,2), r4gamma)
    
    r <- rbind(r1,r2,r4)
    return(r)
  }

  theta.h<-c(mu1.h, mu0.h, gamma.h)
  
  # compute variance
  # Atheta<-jacobian(mphi,theta.h)
  Atheta<-Gradient(theta.h)
  Btheta<-Omega(theta.h)
  invAtheta<-solve(Atheta)
  Vtheta<-invAtheta%*%Btheta%*%t(invAtheta)/N
  
  a<-c(1,-1,rep(0,pp))
  vdelta<-t(a)%*%Vtheta%*%a
  se<-sqrt(vdelta)
  lcl<-delta-qnorm(0.975)*se
  ucl<-delta+qnorm(0.975)*se
  
  # output to screen
  return(list(delta=delta,se=se,lcl=lcl,ucl=ucl))
}
