rm(list=ls())

set.seed(1234)

## fct definitions and loading some parameters
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
simulator = function(v,NPATHS=20) 
{
  #unnormalize
  long = v[1] * (domain[2]-domain[1]) + domain[1]
  lat = v[2] * (domain[4]-domain[3]) + domain[3]
  #KX = v[3] * 900 + 100
  #KY = v[4] * 900 + 100
  KX = 700
  KY = 200
    
  sim = fk_simulator(domain, matrix(c(long,lat),nrow=1), U, V, KX,KY, OXY, NPATHS)
  return(mean(sim))
}


# space filling design
library(DiceDesign)
n = 50
X = lhsDesign(n,2)
X = maximinSA_LHS(X$design)$design
plot(X)

# running simulator with ten replicates

Nrep = 20
out = matrix(NA,n,Nrep)
for (i in 1:n)
  for (j in 1:Nrep)
    out[i,j] = simulator(X[i,])
    
# out = sapply(rep(1:nrow(X),each=Nrep), function(i)
#   {simulator(X[i,])})
# out = t(out)
# dim(out)

# normalization
load("testdesign2D.Rdata")
n=nrow(out)
Mnorm = mean(vm)
SDnorm = sd(vm)
out=(((out)-Mnorm)/SDnorm)
dim(out)

meanOUT = apply(out,1,mean)

library(hetGP)
dimX = 2

# Bobby's settings 
lower <- rep(0.01, 2)
upper <- rep(30, 2)
covtype <- "Matern5_2"
nc <- list(g_min=1e-6, g_bounds=c(1e-6, 1), # k_theta_g_bounds=c(1,10), 
           lowerDelta=log(1e-6))
settings <- list(linkThetas="none", logN=TRUE, initStrategy="smoothed", 
                 checkHom=TRUE, penalty=TRUE, trace=0, return.matrices=TRUE, return.hom=FALSE)
control <- list(tol_dist=1e-4, tol_diff=1e-4, multi.start=30)


Z= matrix((out),Nrep*nrow(X),1,byrow = T)
Xrep = X
for (i in 1:(Nrep-1))
{Xrep=rbind(Xrep,X)}

# Homoskedastic
#Ghom = mleHomGP(list(X0=X,Z0=meanOUT,mult=rep(Nrep,nrow(X))),Z= matrix(t(out),Nrep*nrow(X),1,byrow = T),lower = rep(0.01, dimX), upper = rep(10, dimX),covtype = "Gaussian",maxit=1000) 
Ghom = mleHomGP(Xrep, Z, lower=lower, upper=upper, covtype=covtype,
         known=list(beta0=0), maxit=1000) 

# Heteroskedastic 
#Ghet = mleHetGP(list(X0=X,Z0=meanOUT,mult=rep(Nrep,nrow(X))),Z= matrix(t(out),Nrep*nrow(X),1,byrow = T),lower = rep(0.01, dimX), upper = rep(10, dimX),covtype = "Gaussian",maxit=1000)
Ghet =  mleHetGP(Xrep, Z, lower=lower, upper=upper, covtype=covtype, 
                 noiseControl=nc, settings=settings, known=list(beta0=0), maxit=1000)


# plots 
de = as.matrix(expand.grid(seq(0,1,.01),seq(0,1,.01)))
predHomde = predict(x = de, object = Ghom)
gridHom = as.data.frame(cbind(de[,1:2],predHomde$mean,predHomde$sd2,predHomde$nugs))
names(gridHom) = c("long","lat","mean","varmean","nug")
predhetde = predict(x = de, object = Ghet)
gridhet = as.data.frame(cbind(de[,1:2],predhetde$mean,predhetde$sd2,predhetde$nugs))
names(gridhet) = c("long","lat","mean","varmean","nug")

library(ggplot2)
library(colorRamps)

ggplot(gridHom, aes(long, lat)) + geom_raster(aes(fill = mean), interpolate = TRUE) + scale_fill_gradientn(colours=matlab.like(10),limits=c(min(gridHom$mean,gridhet$mean),max(gridHom$mean,gridhet$mean)))
#scale_fill_gradient(low="blue", high="red",limits=c(min(gridHom$mean,gridhet$mean),max(gridHom$mean,gridhet$mean)))
ggplot(gridHom, aes(long, lat)) + geom_raster(aes(fill = varmean), interpolate = TRUE) + scale_fill_gradientn(colours=matlab.like(10),limits=c(min(gridHom$varmean,gridhet$varmean),max(gridHom$varmean,gridhet$varmean)))

ggplot(gridhet, aes(long, lat)) + geom_raster(aes(fill = mean), interpolate = TRUE) +  scale_fill_gradientn(colours=matlab.like(10),limits=c(min(gridHom$mean,gridhet$mean),max(gridHom$mean,gridhet$mean)))
ggplot(gridhet, aes(long, lat)) + geom_raster(aes(fill = varmean), interpolate = TRUE) + scale_fill_gradientn(colours=matlab.like(10),limits=c(min(gridHom$varmean,gridhet$varmean),max(gridHom$varmean,gridhet$varmean)))
ggplot(gridhet, aes(long, lat)) + geom_raster(aes(fill = nug), interpolate = TRUE) + scale_fill_gradientn(colours=matlab.like(10))



islice = which(gridHom$long==.5)
ggplot(gridHom[islice,],aes(x=lat,y=mean)) +theme_bw()+
  geom_ribbon(aes(ymin = mean-2*sqrt(varmean+nug), ymax = mean+2*sqrt(varmean+nug),fill="prediction"))+geom_ribbon(aes(ymin = mean-2*sqrt(varmean), ymax = mean+2*sqrt(varmean),fill="confidence"))+scale_fill_manual("",values=c("steelblue","lightpink"))+geom_line(aes(x=lat,y=mean,colour="mean"))+scale_colour_manual("",values="black")

ggplot(gridhet[islice,],aes(x=lat,y=mean)) +theme_bw()+
  geom_ribbon(aes(ymin = mean-2*sqrt(varmean+nug), ymax = mean+2*sqrt(varmean+nug),fill="prediction"))+geom_ribbon(aes(ymin = mean-2*sqrt(varmean), ymax = mean+2*sqrt(varmean),fill="confidence"))+scale_fill_manual("",values=c("steelblue","lightpink"))+geom_line(aes(x=lat,y=mean,colour="mean"))+scale_colour_manual("",values="black")

Xinit = X

# seq design
nadd=500
#mod = Ghet
Y = matrix((out),Nrep*nrow(X),1,byrow = T)
h <- rep(NA, nadd)
#dim(Y)

Xseq <- Xrep[1:(10*n),]
Y <- Z[1:(10*n)]
mod <- mleHetGP(Xseq, Y, lower=lower, upper=upper, covtype=covtype, 
                noiseControl=nc, settings=settings, known=list(beta0=0), maxit=1000)


for(i in 1:nadd) { 
  h[i] <- horizon(mod)
  opt <- IMSPE_optim(mod, h[i])
  cat("i=", i, ", h=", h[i], "\n", sep="")
  Xseq <- rbind(Xseq, opt$par)
  Ynew <- (simulator(opt$par)-Mnorm)/SDnorm
  Y <- c(Y, Ynew)
  mod <- update(mod, Xnew = opt$par, Znew = Ynew, ginit = mod$g * 1.01)
  if(i %% 25 == 0){
    mod2 <- mleHetGP(X = list(X0 = mod$X0, Z0 = mod$Z0, mult = mod$mult),Z = mod$Z, lower=lower, upper=upper, covtype=covtype, 
                     noiseControl=nc, settings=settings, known=list(beta0=0), 
                     maxit=1000)
    if(mod2$ll > mod$ll) mod <- mod2
  }
  #print(g)
}
#save.image(file="Designseq2Dbisafterseq500.Rdata")

#save.image(file="designseq2D10paths2.Rdata")



designseq = as.data.frame(cbind(mod$X0,mod$mult))
names(designseq) = c("long","lat","rep")
predhetseqde = predict(x = de, object = mod)
gridhetseq = as.data.frame(cbind(de[,1:2],predhetseqde$mean,predhetseqde$sd2,predhetseqde$nugs, sqrt(predhetseqde$sd2 + predhetseqde$nugs)))
names(gridhetseq) = c("long","lat","mean","varmean","nug","psd")

#save.image(file="Designseq2Dds.Rdata")


g <- ggplot(gridhetseq, aes(long, lat)) + geom_raster(aes(fill = psd), interpolate = TRUE) +scale_fill_gradientn(colours=matlab.like(10)) +geom_point(data=designseq,aes(x=long,y=lat,colour=rep))+scale_color_gradient(low="grey", high="black")
print(g)


nrow(mod$X0)

mod$mult

designseq = as.data.frame(cbind(mod$X0,mod$mult))
names(designseq) = c("long","lat","rep")

predhetde = predict(x = de, object = mod)
gridhet = as.data.frame(cbind(de[,1:2],predhetde$mean,predhetde$sd2,predhetde$nugs, sqrt(predhetde$sd2 + predhetde$nugs)))
names(gridhet) = c("long","lat","mean","varmean","nug","psd")



ggplot(gridhet, aes(long, lat)) + geom_raster(aes(fill = psd), interpolate = TRUE) + scale_fill_gradient(low="blue", high="red")+geom_point(data=designseq,aes(x=long,y=lat,colour=rep))+scale_color_gradient(low="grey", high="black")


