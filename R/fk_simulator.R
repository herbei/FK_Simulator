

setwd("~/Dropbox/0Work/SAMSI/MUMS/FK/R/")
rm(list=ls())
xy=c(-34.0, 11.0, -32.0, -3.0)

lat0=-17.5
domain=rep(0,4)
domain[1]=xy[1]*111e3*cos(lat0*pi/180.0)
domain[2]=xy[2]*111e3*cos(lat0*pi/180.0)
domain[3]=xy[3]*111e3
domain[4]=xy[4]*111e3


N = 37
M = 19

delx=(domain[2]-domain[1])/(N-1)
dely=(domain[4]-domain[3])/(M-1)

# velocity vectors
U = read.table('gridu.txt')
U =-U
V = read.table('gridv.txt')
V = -V

# difussion coefficients
KX = 700
KY=200



NOBS = 200

# generate some spatial locations where simulator will run
sites= matrix(0, nrow=NOBS, ncol=2)
for(i in 1:NOBS){
  sites[i,1] = runif(1, min = domain[1], max = domain[2])
  sites[i,2] = runif(1, min = domain[3], max = domain[4])
}


#this is needed only for the boundary values
OXY = as.matrix(read.table('oxygengrid.txt'))

NPATHS = 200


#run FK simulator
start_time <- Sys.time()
OUT = fk_simulator(domain, sites, U, V, KX, KY, OXY, NPATHS)
Sys.time() - start_time


OXYFK = matrix(0, nrow = NOBS, ncol = 3)
for(i in 1:NOBS){
  OXYFK[i,1:2] = sites[i,]
  OXYFK[i,3] = mean(OUT[i,])
}

library(MASS)
write.matrix(OXYFK, "oxygen_fk.txt", sep = ",")

# obtaining plot
OXYFKdf = as.data.frame(OXYFK)
names(OXYFKdf) = c("long","lat","oxy")
library(ggplot2)
plotfk = ggplot(OXYFKdf, aes(x=long, y=lat, color=oxy)) + geom_point()+scale_color_gradient(low="blue", high="red")
plotfk


## - thanks Pierre.
OXYFKdf = as.data.frame(OXYFK)
names(OXYFKdf) = c("long","lat","oxy")
library(ggplot2)
plotfk = ggplot(OXYFKdf, aes(x=long, y=lat, color=oxy)) + geom_point()+scale_color_gradient(low="blue", high="red")
plotfk




####################################################################
fk_simulator  = function(domain, sites, U, V, KX, KY, OXY, NPATHS){
  
  
  
  M = dim(OXY)[1]
  N = dim(OXY)[2]
  NOBS = dim(sites)[1]
  HH    = 1e7
  lam   = 1e-11
  

  OUT = matrix(0, nrow = NOBS, ncol = NPATHS)
    
  for(i in 1:NOBS){
    print(i)
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
  

