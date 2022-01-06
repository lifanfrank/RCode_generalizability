#########################################################
# Outcome regression estimators for generalizability
# estimated treatment PS is used
#########################################################
# Input: 
# Y: (trial) outcome
# S: (both trial and cohort) participation indicator
# Ze: (trial) design matrix (including intercept) for propensity 
# Zo: (both trial and cohort) design matrix (including intercept) for outcome
# X: (trial) treatment indicator
# N: target population size

OR1=function(Y,S,Ze,Zo,X,N){
  
  # derive weights^-1 for sampling score model
  n <- sum(S); nt <- (N-n); m <- sum(1-S)
  pw <- rep(1, length(S))
  pw[S==0] <- m/nt
  
  # fit outcome models
  Y1 <- Y[X==1]
  Y0 <- Y[X==0]
  Ze1 <- Ze[S==1,]        # in the trial
  Zo1 <- Zo[S==1,]        # in the trial
  Ze11 <- Ze1[X==1,]      # in the trial and treated
  Zo11 <- Zo1[X==1,]      # in the trial and treated
  Ze10 <- Ze1[X==0,]      # in the trial but untreated
  Zo10 <- Zo1[X==0,]      # in the trial but untreated
  
  #estimate treatment propensity scores using logistic regression
  mylogit2 <- glm(X ~  -1+Ze1 , family = "binomial")
  beta.h <- as.numeric(coef(mylogit2))
  e<-plogis(c(Ze1%*%beta.h))
  e1<-e[X==1]
  e0<-e[X==0]
  
  om.fit1 <- lm(Y1 ~ -1+Zo11, weights=e1^(-1))
  om.fit0 <- lm(Y0 ~ -1+Zo10, weights=(1-e0)^(-1))
  
  # predicted potential outcomes
  alpha1.h <- as.numeric(coef(om.fit1))
  alpha0.h <- as.numeric(coef(om.fit0))
  
  m1.h <- c(Zo %*% alpha1.h)
  m0.h <- c(Zo %*% alpha0.h)
  m1.h1 <- m1.h[S==1]   # predicted Y^1 in the trial
  m0.h1 <- m0.h[S==1]   # predicted Y^0 in the trial
  
  # point estimation
  c.h<-(S+(1-S)/pw)
  delta <- sum(c.h*(m1.h-m0.h)) / N
  
  # get the dimensions right first
  rr1 <- ncol(model.matrix(om.fit1))
  rr0 <- ncol(model.matrix(om.fit0))
  qq<-ncol(model.matrix(mylogit2))
  
  # functions to compute estimating equations
  # here we output a (1+rr1+rr0+qq) x N matrix
  phi<-function(theta){
    nu3<-theta[1]
    alpha1<-theta[2:(1+rr1)]
    alpha0<-theta[(2+rr1):(1+rr1+rr0)]
    beta<-theta[(2+rr1+rr0):(1+rr1+rr0+qq)]
    
    # scores
    e<-plogis(c(Ze1%*%beta))
    
    # outcomes
    m1.h <- c(Zo %*% alpha1)
    m0.h <- c(Zo %*% alpha0)
    m1.h1 <- m1.h[S==1]   # predicted Y^1 in the trial
    m0.h1 <- m0.h[S==1]   # predicted Y^0 in the trial
    
    # components
    f3 <- c.h*(m1.h-m0.h)
    f3t <- c(f3,rep(0,N-(n+m)))-rep(nu3,N)
    
    f5 <- matrix(0,n+m,rr1)
    f6 <- matrix(0,n+m,rr0)
    f5[S==1,] <- Zo1*(Y-m1.h1)*X/e
    f6[S==1,] <- Zo1*(Y-m0.h1)*(1-X)/(1-e)
    f5t <- rbind(f5, matrix(0,N-(n+m),rr1))
    f6t <- rbind(f6, matrix(0,N-(n+m),rr0))
    
    fe<-matrix(0,n+m,qq)
    fe[S==1,]<-Ze1*(X-e)
    fet<-rbind(fe, matrix(0,N-(n+m),qq))
    
    f <- rbind(f3t, t(f5t), t(f6t), t(fet))
    return(f)
  }
  
  # score operator
  mphi<-function(theta){
    as.numeric(rowMeans(phi(theta)))
  }
  
  # covariance operator
  Omega<-function(theta){
    phis<-phi(theta)
    return(tcrossprod(phis)/N)
  }
  
  # differentiation operator
  Gradient<-function(theta){
    nu3<-theta[1]
    alpha1<-theta[2:(1+rr1)]
    alpha0<-theta[(2+rr1):(1+rr1+rr0)]
    beta<-theta[(2+rr1+rr0):(1+rr1+rr0+qq)]
    
    # scores
    e<-plogis(c(Ze1%*%beta))
    
    # outcomes
    m1.h <- c(Zo %*% alpha1)
    m0.h <- c(Zo %*% alpha0)
    m1.h1 <- m1.h[S==1]   # predicted Y^1 in the trial
    m0.h1 <- m0.h[S==1]   # predicted Y^0 in the trial
    
    # components
    r3alpha1<- colSums(c.h*Zo)/N
    r3alpha0<- colSums(-c.h*Zo)/N
    r3<-c(-1,c(r3alpha1, r3alpha0),rep(0,qq))
    
    r5alpha1<- -crossprod(X*Zo1/sqrt(e))/N
    r6alpha0<- -crossprod((1-X)*Zo1/sqrt(1-e))/N
    
    # Zcheck1<- Zo1*(X*(1-e)/e)
    # Zcheck0<- Zo1*((1-X)*e/(1-e))
    # r5beta<- -crossprod(Zcheck1, Ze1*(Y-m1.h1))/N
    # r6beta<- crossprod(Zcheck0, Ze1*(Y-m0.h1))/N
    
    Zcheck1<- Zo1*(Y-m1.h1)
    Zcheck0<- Zo1*(Y-m0.h1)
    r5beta<- -crossprod(Zcheck1, Ze1*(X*(1-e)/e))/N
    r6beta<- crossprod(Zcheck0, Ze1*((1-X)*e/(1-e)))/N
    
    r5<- cbind(matrix(0,rr1,1), r5alpha1, matrix(0,rr1,rr0),r5beta)
    r6<- cbind(matrix(0,rr0,1+rr1), r6alpha0,r6beta)
    
    Ztilde<-Ze1*c(sqrt(e*(1-e)))
    rebeta<- -crossprod(Ztilde)/N
    re<- cbind(matrix(0,qq,1+rr1+rr0), rebeta)
    
    r <- rbind(r3,r5,r6,re)
    return(r)
  }
  
  theta.h<-c(delta, alpha1.h, alpha0.h, beta.h)
  
  # compute variance
  # Atheta<-jacobian(mphi,theta.h)
  Atheta<-Gradient(theta.h)
  Btheta<-Omega(theta.h)
  invAtheta<-solve(Atheta)
  Vtheta<-invAtheta%*%Btheta%*%t(invAtheta)/N
  
  a<-c(1,rep(0,rr1+rr0+qq))
  vdelta<-t(a)%*%Vtheta%*%a
  se<-sqrt(vdelta)
  lcl<-delta-qnorm(0.975)*se
  ucl<-delta+qnorm(0.975)*se
  
  # output to screen
  return(list(delta=delta,se=se,lcl=lcl,ucl=ucl))
}
