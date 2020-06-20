##############################################
# Outcome regression estimators for generalizability
# true treatment PS is used
##############################################
# Y: (trial) outcome
# S: (both trial and cohort) participation indicator
# Z: (both trial and cohort) design matrix (including intercept) for outcome model
# X: (trial) treatment indicator
# N: target population size

OR1=function(Y,S,Z,X,N){
  
  # derive weights^-1 for sampling score model
  n <- sum(S); nt <- (N-n); m <- sum(1-S)
  pw <- rep(1, length(S))
  pw[S==0] <- m/nt
  
  # fit outcome models
  Y1 <- Y[X==1]
  Y0 <- Y[X==0]
  Z1 <- Z[S==1,]      # in the trial
  Z11 <- Z1[X==1,]    # in the trial and treated
  Z10 <- Z1[X==0,]    # in the trial but untreated
  
  om.fit1 <- lm(Y1 ~ -1+Z11)
  om.fit0 <- lm(Y0 ~ -1+Z10)
  
  # predicted potential outcomes
  alpha1.h <- as.numeric(coef(om.fit1))
  alpha0.h <- as.numeric(coef(om.fit0))
  
  m1.h <- c(Z %*% alpha1.h)
  m0.h <- c(Z %*% alpha0.h)
  m1.h1 <- m1.h[S==1]   # predicted Y^1 in the trial
  m0.h1 <- m0.h[S==1]   # predicted Y^0 in the trial
  
  # point estimation
  c.h<-(S+(1-S)/pw)
  delta <- sum(c.h*(m1.h-m0.h)) / N
  
  # get the dimensions right first
  rr1 <- ncol(model.matrix(om.fit1))
  rr0 <- ncol(model.matrix(om.fit0))
  
  # functions to compute estimating equations
  # here we output a (1+rr1+rr0) x N matrix
  phi<-function(theta){
    nu3<-theta[1]
    alpha1<-theta[2:(1+rr1)]
    alpha0<-theta[(2+rr1):(1+rr1+rr0)]
    
    # outcomes
    m1.h <- c(Z %*% alpha1)
    m0.h <- c(Z %*% alpha0)
    m1.h1 <- m1.h[S==1]   # predicted Y^1 in the trial
    m0.h1 <- m0.h[S==1]   # predicted Y^0 in the trial
    
    # components
    f3 <- c.h*(m1.h-m0.h)
    f3t <- c(f3,rep(0,N-(n+m)))-rep(nu3,N)
    
    f5 <- matrix(0,n+m,rr1)
    f6 <- matrix(0,n+m,rr0)
    f5[S==1,] <- Z1*(Y-m1.h1)*X
    f6[S==1,] <- Z1*(Y-m0.h1)*(1-X)
    f5t <- rbind(f5, matrix(0,N-(n+m),rr1))
    f6t <- rbind(f6, matrix(0,N-(n+m),rr0))
    
    f <- rbind(f3t, t(f5t), t(f6t))
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
    
    # outcomes
    m1.h <- c(Z %*% alpha1)
    m0.h <- c(Z %*% alpha0)
    m1.h1 <- m1.h[S==1]   # predicted Y^1 in the trial
    m0.h1 <- m0.h[S==1]   # predicted Y^0 in the trial
    
    # components
    r3alpha1<- colSums(c.h*Z)/N
    r3alpha0<- colSums(-c.h*Z)/N
    r3<-c(-1,c(r3alpha1, r3alpha0))
     
    r5alpha1<- -crossprod(X*Z1)/N
    r6alpha0<- -crossprod((1-X)*Z1)/N
    r5<- cbind(matrix(0,rr1,1), r5alpha1, matrix(0,rr1,rr0))
    r6<- cbind(matrix(0,rr0,1+rr1), r6alpha0)
    
    r <- rbind(r3,r5,r6)
    return(r)
  }
  
  theta.h<-c(delta, alpha1.h, alpha0.h)
  
  # compute variance
  # Atheta<-jacobian(mphi,theta.h)
  Atheta<-Gradient(theta.h)
  Btheta<-Omega(theta.h)
  invAtheta<-solve(Atheta)
  Vtheta<-invAtheta%*%Btheta%*%t(invAtheta)/N
  
  a<-c(1,rep(0,rr1+rr0))
  vdelta<-t(a)%*%Vtheta%*%a
  se<-sqrt(vdelta)
  lcl<-delta-qnorm(0.975)*se
  ucl<-delta+qnorm(0.975)*se
  
  # output to screen
  return(list(delta=delta,se=se,lcl=lcl,ucl=ucl))
}
