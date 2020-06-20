##############################################
# Doubly-Robust estimators
# DR1 and DR2 for generalizability
# estimated treatment PS is used 
##############################################
# Y: (trial) outcome
# S: (both trial and cohort) participation indicator
# Zp: (both trial and cohort) design matrix (including intercept) for participation
# Ze: (trial) design matrix (including intercept) for treatment propensity 
# Zo: (both trial and cohort) design matrix (including intercept) for outcome models
# X: (trial) treatment indicator
# N: target population size

DR1=function(Y,S,Zp,Ze,Zo,X,N){
  
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
  Zo1 <- Zo[S==1,]        # in the trial
  Zp11 <- Zp1[X==1,]      # in the trial and treated
  Ze11 <- Ze1[X==1,]      # in the trial and treated
  Zo11 <- Zo1[X==1,]      # in the trial and treated
  Zp10 <- Zp1[X==0,]      # in the trial but untreated
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
  nu1.h <- sum(X*(Y-m1.h1)/(p1*e)) / N
  nu2.h <- sum((1-X)*(Y-m0.h1)/(p1*(1-e))) / N
  nu3.h <- sum(c.h*(m1.h-m0.h)) / N
  delta <- nu1.h - nu2.h + nu3.h
  
  # get the dimensions right first
  pp <- ncol(model.matrix(mylogit))
  rr1 <- ncol(model.matrix(om.fit1))
  rr0 <- ncol(model.matrix(om.fit0))
  qq<-ncol(model.matrix(mylogit2))
  
  # functions to compute estimating equations
  # here we output a (3+pp+rr1+rr0+qq) x N matrix
  phi <- function(theta){
    nu1 <- theta[1]
    nu2 <- theta[2]
    nu3 <- theta[3]
    gamma <- theta[4:(3+pp)]
    alpha1 <- theta[(4+pp):(3+pp+rr1)]
    alpha0 <- theta[(4+pp+rr1):(3+pp+rr1+rr0)]
    beta<-theta[(4+pp+rr1+rr0):(3+pp+rr1+rr0+qq)]
    
    # scores
    p <- plogis(c(Zp%*%gamma))
    p1 <- p[S==1]
    e<-plogis(c(Ze1%*%beta))
    
    m1.h <- c(Zo %*% alpha1)
    m0.h <- c(Zo %*% alpha0)
    m1.h1 <- m1.h[S==1]   # predicted Y^1 in the trial
    m0.h1 <- m0.h[S==1]   # predicted Y^0 in the trial
    
    # components
    f1 <- f2 <- rep(0,n+m)
    f1[S==1] <- X*(Y-m1.h1)/(p1*e)
    f2[S==1] <- (1-X)*(Y-m0.h1)/(p1*(1-e))
    f3 <- c.h*(m1.h-m0.h)
    f1t <- c(f1,rep(0,N-(n+m)))-rep(nu1,N)
    f2t <- c(f2,rep(0,N-(n+m)))-rep(nu2,N)
    f3t <- c(f3,rep(0,N-(n+m)))-rep(nu3,N)
    
    f4 <- Zp*(S-p)*(pw)^(-1)
    f4t <- rbind(f4, matrix(0,N-(n+m),pp))
    
    f5 <- matrix(0,n+m,rr1)
    f6 <- matrix(0,n+m,rr0)
    f5[S==1,] <- Zo1*(Y-m1.h1)*X/e
    f6[S==1,] <- Zo1*(Y-m0.h1)*(1-X)/(1-e)
    f5t <- rbind(f5, matrix(0,N-(n+m),rr1))
    f6t <- rbind(f6, matrix(0,N-(n+m),rr0))
    
    fe<-matrix(0,n+m,qq)
    fe[S==1,]<-Ze1*(X-e)
    fet<-rbind(fe, matrix(0,N-(n+m),qq))
    
    f <- rbind(f1t, f2t, f3t, t(f4t), t(f5t), t(f6t),t(fet))
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
    nu1 <- theta[1]
    nu2 <- theta[2]
    nu3 <- theta[3]
    gamma <- theta[4:(3+pp)]
    alpha1 <- theta[(4+pp):(3+pp+rr1)]
    alpha0 <- theta[(4+pp+rr1):(3+pp+rr1+rr0)]
    beta<-theta[(4+pp+rr1+rr0):(3+pp+rr1+rr0+qq)]
    
    # scores
    p <- plogis(c(Zp%*%gamma))
    p1 <- p[S==1]
    e<-plogis(c(Ze1%*%beta))
    
    m1.h <- c(Zo %*% alpha1)
    m0.h <- c(Zo %*% alpha0)
    m1.h1 <- m1.h[S==1]   # predicted Y^1 in the trial
    m0.h1 <- m0.h[S==1]   # predicted Y^0 in the trial
    
    # components
    r1gamma<- colSums(-X*(Y-m1.h1)*Zp1*(1-p1)/(p1*e))/N
    r2gamma<- colSums(-(1-X)*(Y-m0.h1)*Zp1*(1-p1)/(p1*(1-e)))/N
    r1alpha1<- colSums(-X*Zo1/(p1*e))/N
    r2alpha0<- colSums(-(1-X)*Zo1/(p1*(1-e)))/N
    r3alpha1<- colSums(c.h*Zo)/N
    r3alpha0<-colSums(-c.h*Zo)/N
    
    r1beta<- colSums(-X*(Y-m1.h1)*Ze1*(1-e)/(p1*e))/N
    r2beta<- colSums((1-X)*(Y-m0.h1)*Ze1*e/(p1*(1-e)))/N
    
    r1<-c(-1,0,0,c(r1gamma,r1alpha1,rep(0,rr0),r1beta))
    r2<-c(0,-1,0,c(r2gamma,rep(0,rr1),r2alpha0,r2beta))
    r3<-c(0,0,-1,c(rep(0,pp),r3alpha1,r3alpha0,rep(0,qq)))
    
    Zstar<-Zp*c(sqrt(p*(1-p)*pw^(-1)))
    r4gamma<- -crossprod(Zstar)/N
    r4<- cbind(matrix(0,pp,3), r4gamma, matrix(0,pp,rr1+rr0+qq))
    
    r5alpha1<- -crossprod(X*Zo1/sqrt(e))/N
    r6alpha0<- -crossprod((1-X)*Zo1/sqrt(1-e))/N
    Zcheck1<- Zo1*(Y-m1.h1)
    Zcheck0<- Zo1*(Y-m0.h1)
    r5beta<- -crossprod(Zcheck1, Ze1*(X*(1-e)/e))/N
    r6beta<- crossprod(Zcheck0, Ze1*((1-X)*e/(1-e)))/N
    
    r5<- cbind(matrix(0,rr1,3+pp), r5alpha1, matrix(0,rr1,rr0),r5beta)
    r6<- cbind(matrix(0,rr0,3+pp+rr1), r6alpha0,r6beta)
    
    Ztilde<-Ze1*c(sqrt(e*(1-e)))
    rebeta<- -crossprod(Ztilde)/N
    re<- cbind(matrix(0,qq,3+pp+rr1+rr0), rebeta)
    
    r <- rbind(r1,r2,r3,r4,r5,r6,re)
    return(r)
  }
  
  theta.h <- c(nu1.h, nu2.h, nu3.h, gamma.h, alpha1.h, alpha0.h, beta.h)
  
  # compute variance
  # Atheta<-jacobian(mphi,theta.h)
  Atheta<-Gradient(theta.h)
  Btheta<-Omega(theta.h)
  invAtheta<-solve(Atheta)
  Vtheta<-invAtheta%*%Btheta%*%t(invAtheta)/N
  
  a<-c(1,-1,1,rep(0,pp+rr1+rr0+qq))
  vdelta<-t(a)%*%Vtheta%*%a
  se<-sqrt(vdelta)
  lcl<-delta-qnorm(0.975)*se
  ucl<-delta+qnorm(0.975)*se
  
  # output to screen
  return(list(delta=delta,se=se,lcl=lcl,ucl=ucl))
}



DR2=function(Y,S,Zp,Ze,Zo,X,N){
  
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
  
  # true propensity scores
  e <- rep(0.5, length(X))
  
  # fit outcome models
  Y1 <- Y[X==1]
  Y0 <- Y[X==0]
  Zp1 <- Zp[S==1,]        # in the trial
  Ze1 <- Ze[S==1,]        # in the trial
  Zo1 <- Zo[S==1,]        # in the trial
  Zp11 <- Zp1[X==1,]      # in the trial and treated
  Ze11 <- Ze1[X==1,]      # in the trial and treated
  Zo11 <- Zo1[X==1,]      # in the trial and treated
  Zp10 <- Zp1[X==0,]      # in the trial but untreated
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
  nu1.h <- sum(X*(Y-m1.h1)/(p1*e)) / sum(X/(p1*e))
  nu2.h <- sum((1-X)*(Y-m0.h1)/(p1*(1-e))) / sum((1-X)/(p1*(1-e)))
  nu3.h <- sum(c.h*(m1.h-m0.h)) / N
  delta <- nu1.h - nu2.h + nu3.h
  
  # get the dimensions right first
  pp <- ncol(model.matrix(mylogit))
  rr1 <- ncol(model.matrix(om.fit1))
  rr0 <- ncol(model.matrix(om.fit0))
  qq<-ncol(model.matrix(mylogit2))
  
  # functions to compute estimating equations
  # here we output a (3+pp+rr1+rr0+qq) x N matrix
  phi <- function(theta){
    nu1 <- theta[1]
    nu2 <- theta[2]
    nu3 <- theta[3]
    gamma <- theta[4:(3+pp)]
    alpha1 <- theta[(4+pp):(3+pp+rr1)]
    alpha0 <- theta[(4+pp+rr1):(3+pp+rr1+rr0)]
    beta<-theta[(4+pp+rr1+rr0):(3+pp+rr1+rr0+qq)]
    
    # scores
    p <- plogis(c(Zp%*%gamma))
    p1 <- p[S==1]
    e<-plogis(c(Ze1%*%beta))
    
    m1.h <- c(Zo %*% alpha1)
    m0.h <- c(Zo %*% alpha0)
    m1.h1 <- m1.h[S==1]   # predicted Y^1 in the trial
    m0.h1 <- m0.h[S==1]   # predicted Y^0 in the trial
    
    # components
    f1 <- f2 <- rep(0,n+m)
    f1[S==1] <- X*(Y-m1.h1-nu1)/(p1*e)
    f2[S==1] <- (1-X)*(Y-m0.h1-nu2)/(p1*(1-e))
    f3 <- c.h*(m1.h-m0.h)
    f1t <- c(f1,rep(0,N-(n+m)))
    f2t <- c(f2,rep(0,N-(n+m)))
    f3t <- c(f3,rep(0,N-(n+m)))-rep(nu3,N)
    
    f4 <- Zp*(S-p)*(pw)^(-1)
    f4t <- rbind(f4, matrix(0,N-(n+m),pp))
    
    f5 <- matrix(0,n+m,rr1)
    f6 <- matrix(0,n+m,rr0)
    f5[S==1,] <- Zo1*(Y-m1.h1)*X/e
    f6[S==1,] <- Zo1*(Y-m0.h1)*(1-X)/(1-e)
    f5t <- rbind(f5, matrix(0,N-(n+m),rr1))
    f6t <- rbind(f6, matrix(0,N-(n+m),rr0))
    
    fe<-matrix(0,n+m,qq)
    fe[S==1,]<-Ze1*(X-e)
    fet<-rbind(fe, matrix(0,N-(n+m),qq))
    
    f <- rbind(f1t, f2t, f3t, t(f4t), t(f5t), t(f6t),t(fet))
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
    nu1 <- theta[1]
    nu2 <- theta[2]
    nu3 <- theta[3]
    gamma <- theta[4:(3+pp)]
    alpha1 <- theta[(4+pp):(3+pp+rr1)]
    alpha0 <- theta[(4+pp+rr1):(3+pp+rr1+rr0)]
    beta<-theta[(4+pp+rr1+rr0):(3+pp+rr1+rr0+qq)]
    
    # scores
    p <- plogis(c(Zp%*%gamma))
    p1 <- p[S==1]
    e<-plogis(c(Ze1%*%beta))
    
    m1.h <- c(Zo %*% alpha1)
    m0.h <- c(Zo %*% alpha0)
    m1.h1 <- m1.h[S==1]   # predicted Y^1 in the trial
    m0.h1 <- m0.h[S==1]   # predicted Y^0 in the trial
    
    # components
    r1gamma<- colSums(-X*(Y-m1.h1-nu1)*Zp1*(1-p1)/(p1*e))/N
    r2gamma<- colSums(-(1-X)*(Y-m0.h1-nu2)*Zp1*(1-p1)/(p1*(1-e)))/N
    r1alpha1<- colSums(-X*Zo1/(p1*e))/N
    r2alpha0<- colSums(-(1-X)*Zo1/(p1*(1-e)))/N
    r3alpha1<- colSums(c.h*Zo)/N
    r3alpha0<-colSums(-c.h*Zo)/N
    
    r1beta<- colSums(-X*(Y-m1.h1-nu1)*Ze1*(1-e)/(p1*e))/N
    r2beta<- colSums((1-X)*(Y-m0.h1-nu2)*Ze1*e/(p1*(1-e)))/N
    
    r1<-c(sum(-X/(p1*e))/N,0,0,c(r1gamma,r1alpha1,rep(0,rr0),r1beta))
    r2<-c(0,sum(-(1-X)/(p1*(1-e)))/N,0,c(r2gamma,rep(0,rr1),r2alpha0,r2beta))
    r3<-c(0,0,-1,c(rep(0,pp),r3alpha1,r3alpha0,rep(0,qq)))
    
    Zstar<-Zp*c(sqrt(p*(1-p)*pw^(-1)))
    r4gamma<- -crossprod(Zstar)/N
    r4<- cbind(matrix(0,pp,3), r4gamma, matrix(0,pp,rr1+rr0+qq))
    
    r5alpha1<- -crossprod(X*Zo1/sqrt(e))/N
    r6alpha0<- -crossprod((1-X)*Zo1/sqrt(1-e))/N
    Zcheck1<- Zo1*(Y-m1.h1)
    Zcheck0<- Zo1*(Y-m0.h1)
    r5beta<- -crossprod(Zcheck1, Ze1*(X*(1-e)/e))/N
    r6beta<- crossprod(Zcheck0, Ze1*((1-X)*e/(1-e)))/N
    
    r5<- cbind(matrix(0,rr1,3+pp), r5alpha1, matrix(0,rr1,rr0),r5beta)
    r6<- cbind(matrix(0,rr0,3+pp+rr1), r6alpha0,r6beta)
    
    Ztilde<-Ze1*c(sqrt(e*(1-e)))
    rebeta<- -crossprod(Ztilde)/N
    re<- cbind(matrix(0,qq,3+pp+rr1+rr0), rebeta)
    
    r <- rbind(r1,r2,r3,r4,r5,r6,re)
    return(r)
  }
  
  theta.h <- c(nu1.h, nu2.h, nu3.h, gamma.h, alpha1.h, alpha0.h, beta.h)
  
  # compute variance
  # Atheta<-jacobian(mphi,theta.h)
  Atheta<-Gradient(theta.h)
  Btheta<-Omega(theta.h)
  invAtheta<-solve(Atheta)
  Vtheta<-invAtheta%*%Btheta%*%t(invAtheta)/N
  
  a<-c(1,-1,1,rep(0,pp+rr1+rr0+qq))
  vdelta<-t(a)%*%Vtheta%*%a
  se<-sqrt(vdelta)
  lcl<-delta-qnorm(0.975)*se
  ucl<-delta+qnorm(0.975)*se
  
  # output to screen
  return(list(delta=delta,se=se,lcl=lcl,ucl=ucl))
}


