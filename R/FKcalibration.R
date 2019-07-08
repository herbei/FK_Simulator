
# packages ----------------------------------------------------------------
library(ggplot2)
library(gridExtra)
library(DiceKriging)
library(DiceDesign)
library(hetGP)
library(MCMCpack)
library(lhs)


# loading FK setting and simulator function -------------------------------

rm(list=ls())
set.seed(123)

load("fksettings.Rdata")
N = 37
M = 19

delx=(domain[2]-domain[1])/(N-1)
dely=(domain[4]-domain[3])/(M-1)

fk_simulator  = function(domain, sites, U, V, KX, KY, OXY, NPATHS){
  
  M = dim(OXY)[1]
  N = dim(OXY)[2]
  NOBS = dim(sites)[1]
  HH    = 1e7
  lam   = 1e-11
  
  OUT = matrix(0, nrow = NOBS, ncol = NPATHS)
  
  for(i in 1:NOBS){
    #print(i)
    for(kk in 1:NPATHS){
      
      
      xc = sites[i,1]
      yc = sites[i,2]
      
      tau=0;
      while( ((xc - domain[1])*(xc-domain[2])<0) & ( (yc-domain[3])*(yc-domain[4])<0)    ){
        
        #--------
        # find indices
        jx = ceiling( (xc - domain[1])/delx );
        ix = ceiling( (yc - domain[3])/dely );
        if (jx<=0){
          jx=1
        }
        else{
          if (jx>N){
            jx=N
          }
        }
        
        if (ix<=0){
          ix=1
        }
        else{
          if (ix>M) ix=M
        }
        #------------
        
        uc = U[ix,jx]
        vc = V[ix,jx]
        
        # Euler step
        xn = xc + HH * uc + sqrt(HH)*sqrt(2*KX)*rnorm(1);
        yn = yc + HH * vc + sqrt(HH)*sqrt(2*KY)*rnorm(1);    
        
        xc=xn;
        yc=yn;
        tau = tau + HH;
        
        
      }
      JJ = ceiling( (xc - domain[1])/delx );
      II = ceiling( (yc - domain[3])/dely );
      
      if (JJ<=0){
        JJ=1
      }
      else{
        if (JJ>N) JJ=N
      }
      
      if (II<=0){
        II=1
      }
      else{
        if (II>M) II=M
      }
      
      OUT[i,kk] = OXY[II,JJ]*exp(-lam*tau)
      
    }
  }
  
  return(OUT)
}



# function to be called by IMSPE_optim
# v has dim 4: 2 first dimensions are long, lat, then KX and KY diffusion coefficients
simulator = function(v,NPATHS=6) 
{
  #unnormalize
  long = v[1] * (domain[2]-domain[1]) + domain[1]
  lat = v[2] * (domain[4]-domain[3]) + domain[3]
  KX = v[3] * 900 + 100
  KY = v[4] * 900 + 100
  #KX = 700
  #KY = 200
  
  sim = fk_simulator(domain, matrix(c(long,lat),nrow=1), U, V, KX,KY, OXY, NPATHS)
  return(mean(sim))
}



# Simulation for field data -----------------------------------------------

# echantillonnage selong un LHS
n = 150
X0 = lhsDesign(n,2)
X0 = maximinSA_LHS(X0$design)$design
u1 = (700 -100)/900  
u2 = (200- 100)/900
Xfield = cbind(X0,u1,u2)

library(parallel)
Yfield = mclapply(1:n,function(i)
{
  mean(sapply(1:200,function(k) simulator(Xfield[i,])))
},mc.cores = 10)
Yfield = Reduce(c,Yfield)

vareps = 2^2 
Yfieldnoise = Yfield + rnorm(length(Yfield),0,sd=sqrt(vareps)) 




# DOE in 4d for GP --------------------------------------------------------

source("function_LHD_optim_grid.R")
b_inf_theta=c(0,0)
b_sup_theta=c(1,1)
size_design = 500

out = GRID_design(Xfield[,1:2],b_inf_theta,b_sup_theta,size_design)
out_opt=maximinSA_LHS_grid(out,dim_xf=2,c=0.95,dim_theta=2,Imax=2000)
#dim(out_opt)
#head(out_opt)
plot(out_opt[,c(1,2)])
points(Xfield,col=2)

#10 replicates
Xsim = matrix(rep(t(out_opt),10),ncol=4,byrow = T)
head(Xsim)
out_opt[1:6,]
Xsim[c(1,501,1001),]
#ok

Ysim = apply(Xsim,1,simulator)






# GPhom and GPhet ---------------------------------------------------------
covtype <- "Matern5_2"
noiseControl <- list(g_min=1e-6, g_bounds=c(1e-6, 1), lowerDelta=log(1e-6))
settings <- list(linkThetas="none", initStrategy="smoothed", return.hom=TRUE)
lower <- c(0.01, 0.01, 0.001, 0.001) 
upper <- c(30, 30, 100, 100) 

het <- mleHetGP(Xsim, Ysim, lower=lower, upper=upper, covtype=covtype, noiseControl=noiseControl, 
                settings=settings, maxit=10000)
hom <- het$modHom



# Calibration -------------------------------------------------------------

library(mvtnorm)
lpost.invert <- function(theta, XF, yF, GP)
{
  ## input processing and checking
  if(length(theta) != ncol(GP$X0) - ncol(XF) + 1) 
    stop("length(theta), ncol(XF), ncol(GP$X0) mismatch")
  u <- theta[-length(theta)]
  s2 <- theta[length(theta)]
  
  ## prior checking  
  if(any(u < 0 | u > 1)) return (-Inf)
  if(s2 < 0) return(-Inf)
  
  ## derive predictive distribution for XF paired with u
  XFU <- cbind(XF, matrix(rep(u, nrow(XF)), ncol=length(u), byrow=TRUE)) 
  p <- predict(GP, XFU, xprime=XFU)
  C <- s2*diag(nrow(p$cov)) + (p$cov + t(p$cov))/2 
  
  ## gaussian log density evaluation for yF under that predictive
  return(dmvnorm(yF, p$mean, C, log=TRUE) - log(s2))
}


V <- diag(c(0.01, 0.005, 0.05))
calHom <- MCMCmetrop1R(lpost.invert, XF=Xfield[,1:2], yF=Yfieldnoise, GP=hom, logfun=TRUE, theta.init=c(0.5,0.5,1), 
                       burnin=1000, mcmc=2e4, V=V)

calHet <- MCMCmetrop1R(lpost.invert, XF=Xfield[,1:2], yF=Yfieldnoise, GP=het, logfun=TRUE, theta.init=c(0.5,0.5,1), 
                       burnin=1000, mcmc=2e4, V=V)


par(mfrow=c(4,3))
plot(as.numeric(calHom[,1]), type="l", ylab="u1", main="hom u1 trace")
plot(as.numeric(calHom[,2]), type="l", ylab="u2", main="hom u2 trace")
plot(as.numeric(calHom[,3]), type="l", ylab="s2", main="hom s2 trace")
hist(calHom[,1], xlim=c(0,1), xlab="u1", main="hom u1 posterior")
abline(v=2/3, lwd=2, col=2)
hist(calHom[,2], xlim=c(0,1), xlab="u2", main="hom u2 posterior")
abline(v=1/9, lwd=2, col=2)
hist(calHom[,3], xlab="s2", main="hom s2 posterior")
plot(as.numeric(calHet[,1]), type="l", ylab="u1", main="het u1 trace")
plot(as.numeric(calHet[,2]), type="l", ylab="u2", main="het u2 trace")
plot(as.numeric(calHet[,3]), type="l", ylab="s2", main="het s2 trace")
hist(calHet[,1], xlim=c(0,1), xlab="u1", main="het u1 posterior")
abline(v=2/3, lwd=2, col=2)
hist(calHet[,2], xlim=c(0,1), xlab="u2", main="het u2 posterior")
abline(v=1/9, lwd=2, col=2)
hist(calHet[,3], xlab="s2", main="het s2 posterior")

save.image("FKcaliboutput2.Rdata")

dfcalhom = as.data.frame(calHom)
names(dfcalhom) = c("u1","u2","s2")
dfcalhom[,1:2] = dfcalhom[,1:2] * 900 +100
hom1 = ggplot(dfcalhom,aes(x=u1,stat(density))) + geom_histogram(bins= 50) + theme_bw() + xlim(50,1000) + geom_vline(aes(xintercept=700),color="red")
hom2 = ggplot(dfcalhom,aes(x=u2,stat(density))) + geom_histogram(bins = 50) + theme_bw() +xlim(50,1000)  + geom_vline(aes(xintercept=200),color="red")
hom3 = ggplot(dfcalhom,aes(x=s2,stat(density))) + geom_histogram(bins = 50) + theme_bw() + xlim(1,7) + geom_vline(aes(xintercept=4),color="red")
dfcalhet = as.data.frame(calHet)
names(dfcalhet) = c("u1","u2","s2")
dfcalhet[,1:2] = dfcalhet[,1:2] * 900 +100
het1 = ggplot(dfcalhet,aes(x=u1,stat(density))) + geom_histogram(bins = 50) + theme_bw() + xlim(50,1000) + geom_vline(aes(xintercept=700),color="red")
het2 = ggplot(dfcalhet,aes(x=u2,stat(density))) + geom_histogram(bins = 50) + theme_bw() + xlim(50,1000) + geom_vline(aes(xintercept=200),color="red")
het3 = ggplot(dfcalhet,aes(x=s2,stat(density))) + geom_histogram(bins=50) + theme_bw() + xlim(1,7) + geom_vline(aes(xintercept=4),color="red")

grid.arrange(hom1,hom2,hom3,het1,het2,het3,nrow=2,ncol=3)
