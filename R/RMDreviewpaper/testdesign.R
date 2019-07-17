rm(list=ls())

set.seed(123)



chemingal="/projet/extern/save/pbarbillon/Radu/"

setwd(chemingal)


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
simulator = function(v,NPATHS=1) 
{
  #unnormalize
  long = v[1] * (domain[2]-domain[1]) + domain[1]
  lat = v[2] * (domain[4]-domain[3]) + domain[3]
  #KX = v[3] * 900 + 100
  #KY = v[4] * 900 + 100
  KX = 700
  KY = 200
  
  fk_simulator(domain, matrix(c(long,lat),nrow=1), U, V, KX,KY, OXY, NPATHS)
}

library(parallel)
simulatorpar = function(v)
{
  Reduce(c,mclapply(1:400, function(i) {return(simulator(v))},mc.cores = 10))
}

#simulatorpar(c(.5,.5))

testdesign = lhs::randomLHS(500,2)
testdesign <- maximinSA_LHS(testdesign)$design
plot(testdesign)
vm = numeric(nrow(testdesign))
vs = numeric(nrow(testdesign))
vn = numeric(nrow(testdesign))



for (i in 1:nrow(testdesign))
{
   res = as.vector(simulatorpar(testdesign[i,]))
   m = mean(res)
   s = sd(res)  
   n = length(res)
   while( s/sqrt(n) >.1 )
   {
     res = c(res,as.vector(simulatorpar(testdesign[i,])))
     s = sd(res)
     n = length(res)
     m = mean(res)
   }
   vm[i] = m
   vs[i] = s
   vn[i] = length(res)
   print(i)
}



save.image(file = "testdesign2D.Rdata")


