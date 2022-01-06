#############################
# DR - (d) incorrect ps 
#          incorrect or
#############################
rm(list = ls())

# source functions
dir<-"D:/A1 Methods Research/Causal Inference/DR Generalizability/R1_Revision/functions/"
source(paste0(dir,"simdata.R"))
source(paste0(dir,"DR_true.R"))

# logistics
gammatype=1
alphatype=2

if(gammatype==1){
  gamma<-c(-7.148, 0.3,0.3,0.3)
} 
if(gammatype==2){
  gamma<-c(-7.698, 0.6,0.6,0.6)
} 
if(gammatype==3){
  gamma<-c(-8.642, 0.9,0.9,0.9)
} 

if(alphatype==1){
  alpha<-c(1,1,1)
  Delta<-2.4
} else{
  alpha<-c(2,2,2)
  Delta<-2.8
}


nsim<-5000
out1<-matrix(NA,nsim,4)
out2<-matrix(NA,nsim,4)
colnames(out1)<-colnames(out2)<-c("EST","SE","LCL","UCL")

for(i in 1:nsim){
  set.seed(100+i)
  data<-simdata(gamma,alpha)
  N<-data$N
  trial<-data$trial
  cohort<-data$cohort
  both<-data$both
  
  Y<-trial$Y                    # (trial) outcome
  S<-both$S                     # (both) sampling indicator
  Z<-cbind(1,both$Z1,both$Z2)   # (both) covariate vector for sampling
  X<-trial$X                    # (trial) treatment indicator
  
  Zp<-cbind(1,both$Z1,both$Z2)
  Zo<-cbind(1,both$Z1,both$Z2)
  
  res1<-DR1(Y,S,Zp,Zo,X,N)
  res2<-DR2(Y,S,Zp,Zo,X,N)
  
  # output
  out1[i,]<-c(res1$delta,res1$se,res1$lcl,res1$ucl)
  out2[i,]<-c(res2$delta,res2$se,res2$lcl,res2$ucl)
  
  # Loop control
  if(i%%50 == 0) print(i)
}

# summary
bias1<-mean(out1[,1]-Delta, na.rm = T)
bias2<-mean(out2[,1]-Delta, na.rm = T)

mcse1<-sd(out1[,1], na.rm = T)
mcse2<-sd(out2[,1], na.rm = T)

avgse1<-mean(out1[,2], na.rm = T)
avgse2<-mean(out2[,2], na.rm = T)

cvg1<-mean((Delta>=out1[,3]) & (Delta<=out1[,4]), na.rm = T)
cvg2<-mean((Delta>=out2[,3]) & (Delta<=out2[,4]), na.rm = T)

row1<-c(bias1,mcse1,avgse1,cvg1)
row2<-c(bias2,mcse2,avgse2,cvg2)
tab<-rbind(row1,row2)
colnames(tab)<-c("Bias","Empirical SE", "Average SE", "Coverage")
print(tab)
print(tab,digits=3)
