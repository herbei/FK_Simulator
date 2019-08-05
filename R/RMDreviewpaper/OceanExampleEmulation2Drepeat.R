## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,cache=TRUE)


## ----packages,message=FALSE----------------------------------------------
library(ggplot2)
library(colorRamps)
library(gridExtra)
library(DiceKriging)
library(DiceDesign)
library(hetGP)
library(MCMCpack)
library(lhs)
library(mvtnorm)


## ----fksettings----------------------------------------------------------
source("SimulatorAndFunctions.R")


## ----seed----------------------------------------------------------------
seed=1234
set.seed(seed)


Nrep = 100
vseed = floor(runif(Nrep,1,1e6)) # for lhsDesign function


library(parallel)
RESrep = mclapply(1:Nrep,
function(k){
  ## ----space filling design------------------------------------------------
n = 50
X0 = lhsDesign(n,2,seed=vseed[k])
X0 = maximinSA_LHS(X0$design)$design
X <- rbind(X0, X0, X0, X0, X0, X0, X0, X0, X0, X0)
X <- rbind(X, X)
Z <- apply(X, 1, simulator)


## ----normalization-------------------------------------------------------
Zm <- mean(Z)
Zv <- var(Z)
Z <- (Z - Zm)/sqrt(Zv)


## ----GP settings---------------------------------------------------------
# settings for GPs
lower <- rep(0.01, 2)
upper <- rep(30, 2)
covtype <- "Matern5_2"
nc <- list(g_min=1e-6, g_bounds=c(1e-6, 1), # k_theta_g_bounds=c(1,10), 
           lowerDelta=log(1e-6))
settings <- list(linkThetas="none", logN=TRUE, initStrategy="smoothed", 
                 checkHom=TRUE, penalty=TRUE, trace=0, return.matrices=TRUE, return.hom=FALSE)
control <- list(tol_dist=1e-4, tol_diff=1e-4, multi.start=30)


## ----GP fit--------------------------------------------------------------
## Homoskedastic
Ghom <- mleHomGP(X, Z, lower=lower, upper=upper, covtype=covtype,
                known=list(beta0=0), maxit=1000) 

## Heteroskedastic 
Ghet <- mleHetGP(X, Z, lower=lower, upper=upper, covtype=covtype, 
                noiseControl=nc, settings=settings, known=list(beta0=0), maxit=1000)



## ----seqdesign,results="hide"--------------------------------------------
ninitseq = 10
n=nrow(X0)
Xseq <- X[1:(ninitseq*n),]
Y <- Z[1:(ninitseq*n)]
mod <- mleHetGP(Xseq, Y, lower=lower, upper=upper, covtype=covtype, 
                noiseControl=nc, settings=settings, known=list(beta0=0), maxit=1000)

nadd = nrow(X)-nrow(Xseq)  

h <- rep(NA, nadd)
## acquisitions
for(i in 1:nadd) { 
  
  ## choose lookahead horizon and solve IMSPE
  h[i] <- horizon(mod)
  opt <- IMSPE_optim(mod, h[i], control=control)
  cat("i=", i, ", h=", h[i], "\n", sep="")
  
  ## evaluate the simulator
  ynew <- (simulator(opt$par) -Zm)/sqrt(Zv)
  
  ## update the fit
  mod <- update(mod, Xnew=opt$par, Znew=ynew, ginit=mod$g*1.01)
  if(i %% 25 == 0){
    mod2 <- mleHetGP(list(X0=mod$X0, Z0=mod$Z0, mult=mod$mult),
                     Z=mod$Z, lower=lower, upper=upper, covtype=covtype, 
                     noiseControl=nc, settings=settings, known=list(beta0=0), 
                     maxit=1000)
    if(mod2$ll > mod$ll) mod <- mod2  
  }
}
seqGhet=mod


## ----load data test2-----------------------------------------------------
load("testdesign2D.Rdata")
#vmNorm = (vm-Zm)/sqrt(Zv)
library(hetGP)
#Hom GP
predhom = predict(x = testdesign, object = Ghom)
#Het GP
predhet = predict(x = testdesign, object = Ghet)
#Het GP seq
predHetseq = predict(x = testdesign, object = mod)
# RMSE
Ztest.true = (vm - Zm)/sqrt(Zv)
msehom = mean(((predhom$mean-Ztest.true))^2)
msehet = mean(((predhet$mean-Ztest.true))^2)
msehetseq = mean(((predHetseq$mean-Ztest.true))^2)
# scores
schom = mean(-(Ztest.true-predhom$mean)^2/(predhom$sd2) -log(predhom$sd2))
schet = mean(-(Ztest.true-predhet$mean)^2/(predhet$sd2) -log(predhet$sd2))
schetseq = mean(-(Ztest.true-predHetseq$mean)^2/(predHetseq$sd2) -log(predHetseq$sd2))
# scores 2
Ztest = (apply(testdesign,1,simulator)-Zm)/sqrt(Zv) # for computing scores
schom2 = mean(-(Ztest-predhom$mean)^2/(predhom$sd2+predhom$nugs) -log(predhom$sd2+predhom$nugs))
schet2 = mean(-(Ztest-predhet$mean)^2/(predhet$sd2+predhet$nugs) -log(predhet$sd2+predhet$nugs))
schetseq2 = mean(-(Ztest-predHetseq$mean)^2/(predHetseq$sd2+predHetseq$nugs) -log(predHetseq$sd2+predHetseq$nugs))


## ----results-------------------------------------------------------------
return(c(rmsehom=msehom,rmsehet=msehet,rmsehetse=msehetseq,scorehom=schom,scorehet=schet,scorehetse=schetseq,scorehom=schom2,scorehet=schet2,scorehetse=schetseq2))
},mc.cores=10)


#save(RESrep,seed,file="RepetitonEmulation.Rdata")


library(ggplot2)
library(plyr)
library(reshape2)


RES = Reduce(rbind,RESrep)
dim(RES)

#MSE
RMSE = melt(RES[,1:3])
names(RMSE) = c("v","GP","MSE")
levels(RMSE$GP)
RMSE$GP =revalue(RMSE$GP, c("rmsehom"="homGP", "rmsehet"="hetGP","rmsehetse"="seqhetGP"))
p1=ggplot(RMSE,aes(x=GP,y=MSE))+theme_bw()+geom_boxplot()


#score Mean
ScoreMean = melt(RES[,4:6])
names(ScoreMean) = c("v","GP","ScoreMean")
levels(ScoreMean$GP)
ScoreMean$GP =revalue(ScoreMean$GP, c("scorehom"="homGP", "scorehet"="hetGP","scorehetse"="seqhetGP"))
p2=ggplot(ScoreMean,aes(x=GP,y=ScoreMean))+theme_bw()+geom_boxplot()

#Score
Score = melt(RES[,7:9])
names(Score) = c("v","GP","Score")
levels(Score$GP)
Score$GP =revalue(Score$GP, c("scorehom"="homGP", "scorehet"="hetGP","scorehetse"="seqhetGP"))
p3=ggplot(Score,aes(x=GP,y=Score))+theme_bw()+geom_boxplot()

pdf(file="BoxplotEmulation.pdf")
grid.arrange(p1,p2,p3,nrow=1)
dev.off()


dfmprob = as.data.frame(mprob)
names(dfmprob)=c("t","Ig","Id","Te","Tp","Te.or.Tp")
meltData <- melt(dfmprob)
names(meltData)=c("input","prob")
p <- ggplot(meltData, aes(input,prob)) 
pdf("PVbox.pdf")
p + geom_boxplot() +theme_bw()+ylab("prob of activeness")
dev.off()