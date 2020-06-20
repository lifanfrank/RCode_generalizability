##############################################
# IPW estimators for generalizability
# IPW1 and IPW2
# estimated treatment PS is used
##############################################
# Y: (trial) outcome
# S: (both trial and cohort) participation indicator
# Zp: (both) design matrix (including intercept) for trial participation
# Ze: (trial) design matrix (including intercept) for propensity 
# X: (trial) treatment indicator
# N: target population size

IPW1=function(Y,S,Zp,Ze,X,N){
  
  # derive weights^-1 for sampling score model
  n <- sum(S); nt <- (N-n); m <- sum(1-S)
  pw <- rep(1, length(S))
  pw[S==0] <- m/nt
  
  #estimate selection propensity scores using weighted logistic regression
  mylogit <- glm(S ~  -1+Zp , family = "binomial",weights=pw^(-1))
  gamma.h <- as.numeric(coef(mylogit))
  
  #selection propensity score w_i
  p <- plogis(c(Zp%*%gamma.h))
  p1 <- p[S==1]
  
  # fit outcome models
  Y1 <- Y[X==1]
  Y0 <- Y[X==0]
  Zp1 <- Zp[S==1,]        # in the trial
  Ze1 <- Ze[S==1,]        # in the trial
  Zp11 <- Zp1[X==1,]      # in the trial and treated
  Ze11 <- Ze1[X==1,]      # in the trial and treated
  Zp10 <- Zp1[X==0,]      # in the trial but untreated
  Ze10 <- Ze1[X==0,]      # in the trial but untreated
  
  #estimate treatment propensity scores using logistic regression
  mylogit2 <- glm(X ~  -1+Ze1 , family = "binomial")
  beta.h <- as.numeric(coef(mylogit2))
  e<-plogis(c(Ze1%*%beta.h))
  e1<-e[X==1]
  e0<-e[X==0]
  
  # point estimation
  mu1.h<-sum((X*Y)/(p1*e))/N
  mu0.h<-sum(((1-X)*Y)/(p1*(1-e)))/N
  delta<-mu1.h-mu0.h
  
  # get the dimensions right first
  pp<-ncol(model.matrix(mylogit))
  qq<-ncol(model.matrix(mylogit2))
  
  # functions to compute estimating equations
  # here we output a (2+pp+qq) x N matrix
  phi<-function(theta){
    mu1<-theta[1]
    mu0<-theta[2]
    gamma<-theta[3:(2+pp)]
    beta<-theta[(3+pp):(2+pp+qq)]
    
    # scores
    p <- plogis(c(Zp%*%gamma))
    p1 <- p[S==1]
    e<-plogis(c(Ze1%*%beta))
    
    # components
    f1<-f2<-rep(0,n+m)
    f1[S==1]<-X*Y/(p1*e)
    f2[S==1]<-(1-X)*Y/(p1*(1-e))
    f1t<-c(f1,rep(0,N-(n+m)))-rep(mu1,N)
    f2t<-c(f2,rep(0,N-(n+m)))-rep(mu0,N)
    
    f4<-Zp*(S-p)*(pw)^(-1)
    f4t<-rbind(f4, matrix(0,N-(n+m),pp))
    
    fe<-matrix(0,n+m,qq)
    fe[S==1,]<-Ze1*(X-e)
    fet<-rbind(fe, matrix(0,N-(n+m),qq))
    
    f <- rbind(f1t,f2t,t(f4t),t(fet))
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
    beta<-theta[(3+pp):(2+pp+qq)]
    
    # scores
    p<-plogis(c(Zp%*%gamma))
    p1<-p[S==1]
    e<-plogis(c(Ze1%*%beta))
    
    # components
    r1gamma<- colSums(-X*Y*Zp1*(1-p1)/(p1*e))/N
    r2gamma<- colSums(-(1-X)*Y*Zp1*(1-p1)/(p1*(1-e)))/N
    r1beta<- colSums(-X*Y*Ze1*(1-e)/(p1*e))/N
    r2beta<- colSums((1-X)*Y*Ze1*e/(p1*(1-e)))/N
    r1<-c(-1,0,c(r1gamma,r1beta))
    r2<-c(0,-1,c(r2gamma,r2beta))
    
    Zstar<-Zp*c(sqrt(p*(1-p)*pw^(-1)))
    r4gamma<- -crossprod(Zstar)/N
    r4<- cbind(matrix(0,pp,2), r4gamma, matrix(0,pp,qq))
    
    Ztilde<-Ze1*c(sqrt(e*(1-e)))
    rebeta<- -crossprod(Ztilde)/N
    re<- cbind(matrix(0,qq,2+pp), rebeta)
    
    r <- rbind(r1,r2,r4,re)
    return(r)
  }
  
  theta.h<-c(mu1.h, mu0.h, gamma.h, beta.h)
  
  # compute variance
  #Atheta<-jacobian(mphi,theta.h)
  Atheta<-Gradient(theta.h)
  Btheta<-Omega(theta.h)
  invAtheta<-solve(Atheta)
  Vtheta<-invAtheta%*%Btheta%*%t(invAtheta)/N
  
  a<-c(1,-1,rep(0,pp+qq))
  vdelta<-t(a)%*%Vtheta%*%a
  se<-sqrt(vdelta)
  lcl<-delta-qnorm(0.975)*se
  ucl<-delta+qnorm(0.975)*se
  
  # output to screen
  return(list(delta=delta,se=se,lcl=lcl,ucl=ucl))
}

IPW2=function(Y,S,Zp,Ze,X,N){
  
  # derive weights^-1 for sampling score model
  n <- sum(S); nt <- (N-n); m <- sum(1-S)
  pw <- rep(1, length(S))
  pw[S==0] <- m/nt
  
  #estimate selection propensity scores using weighted logistic regression
  mylogit <- glm(S ~  -1+Zp , family = "binomial",weights=pw^(-1))
  gamma.h <- as.numeric(coef(mylogit))
  
  #selection propensity score w_i
  p <- plogis(c(Zp%*%gamma.h))
  p1 <- p[S==1]
  
  # fit outcome models
  Y1 <- Y[X==1]
  Y0 <- Y[X==0]
  Zp1 <- Zp[S==1,]        # in the trial
  Ze1 <- Ze[S==1,]        # in the trial
  Zp11 <- Zp1[X==1,]      # in the trial and treated
  Ze11 <- Ze1[X==1,]      # in the trial and treated
  Zp10 <- Zp1[X==0,]      # in the trial but untreated
  Ze10 <- Ze1[X==0,]      # in the trial but untreated
  
  #estimate treatment propensity scores using logistic regression
  mylogit2 <- glm(X ~  -1+Ze1 , family = "binomial")
  beta.h <- as.numeric(coef(mylogit2))
  e<-plogis(c(Ze1%*%beta.h))
  e1<-e[X==1]
  e0<-e[X==0]
  
  # point estimation
  mu1.h<-sum((X*Y)/(p1*e))/sum(X/(p1*e))
  mu0.h<-sum(((1-X)*Y)/(p1*(1-e)))/sum((1-X)/(p1*(1-e)))
  delta<-mu1.h-mu0.h
  
  # get the dimensions right first
  pp <- ncol(model.matrix(mylogit))
  qq<-ncol(model.matrix(mylogit2))
  
  # functions to compute estimating equations
  # here we output a (2+pp) x N matrix
  phi<-function(theta){
    mu1<-theta[1]
    mu0<-theta[2]
    gamma<-theta[3:(2+pp)]
    beta<-theta[(3+pp):(2+pp+qq)]
    
    # scores
    p<-plogis(c(Zp%*%gamma))
    p1<-p[S==1]
    e<-plogis(c(Ze1%*%beta))
    
    # components
    f1<-f2<-rep(0,n+m)
    f1[S==1]<-X*(Y-mu1)/(p1*e)
    f2[S==1]<-(1-X)*(Y-mu0)/(p1*(1-e))
    f1t<-c(f1,rep(0,N-(n+m)))
    f2t<-c(f2,rep(0,N-(n+m)))
    
    f4<-Zp*(S-p)*(pw)^(-1)
    f4t<-rbind(f4, matrix(0,N-(n+m),pp))
    
    fe<-matrix(0,n+m,qq)
    fe[S==1,]<-Ze1*(X-e)
    fet<-rbind(fe, matrix(0,N-(n+m),qq))
    
    f <- rbind(f1t,f2t,t(f4t),t(fet))
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
    beta<-theta[(3+pp):(2+pp+qq)]
    
    # scores
    p<-plogis(c(Zp%*%gamma))
    p1<-p[S==1]
    e<-plogis(c(Ze1%*%beta))
    
    # components
    r1gamma<- colSums(-X*(Y-mu1)*Zp1*(1-p1)/(p1*e))/N
    r2gamma<- colSums(-(1-X)*(Y-mu0)*Zp1*(1-p1)/(p1*(1-e)))/N
    r1beta<- colSums(-X*(Y-mu1)*Ze1*(1-e)/(p1*e))/N
    r2beta<- colSums((1-X)*(Y-mu0)*Ze1*e/(p1*(1-e)))/N
    r1<-c(sum(-X/(p1*e))/N, 0, c(r1gamma,r1beta))
    r2<-c(0, sum(-(1-X)/(p1*(1-e)))/N, c(r2gamma,r2beta))
    
    Zstar<-Zp*c(sqrt(p*(1-p)*pw^(-1)))
    r4gamma<- -crossprod(Zstar)/N
    r4<- cbind(matrix(0,pp,2), r4gamma, matrix(0,pp,qq))
    
    Ztilde<-Ze1*c(sqrt(e*(1-e)))
    rebeta<- -crossprod(Ztilde)/N
    re<- cbind(matrix(0,qq,2+pp), rebeta)
    
    r <- rbind(r1,r2,r4,re)
    return(r)
  }
  
  theta.h<-c(mu1.h, mu0.h, gamma.h, beta.h)
  
  # compute variance
  # Atheta<-jacobian(mphi,theta.h)
  Atheta<-Gradient(theta.h)
  Btheta<-Omega(theta.h)
  invAtheta<-solve(Atheta)
  Vtheta<-invAtheta%*%Btheta%*%t(invAtheta)/N
  
  a<-c(1,-1,rep(0,pp+qq))
  vdelta<-t(a)%*%Vtheta%*%a
  se<-sqrt(vdelta)
  lcl<-delta-qnorm(0.975)*se
  ucl<-delta+qnorm(0.975)*se
  
  # output to screen
  return(list(delta=delta,se=se,lcl=lcl,ucl=ucl))
}
