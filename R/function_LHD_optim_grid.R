###DESIGN OPTIMIZED LHS RESTRICTED THE GRID OF FIELD VALUES###

maximinSA_LHS_grid=function (design,dim_xf,dim_theta, T0 = 10, c = 0.95, it = 2000, p = 50, profile = "GEOM", 
    Imax = 100) 
{
    crit <- NULL
    temp <- NULL
    proba <- NULL
    if (profile == "GEOM") {
        m <- design
        i <- 0
        T <- T0
        fi_p <- phiP(m, p)
        crit <- fi_p
        while (T > 0 & i < it) {
            G <- perturbationLHS(m,dim_xf,dim_theta)
            fi_p_ep <- phiP(G, p)
            diff <- min(exp((fi_p - fi_p_ep)/T), 1)
            if (diff == 1) {
                m <- G
                fi_p <- fi_p_ep
            }
            else {
                Bernoulli <- rbinom(1, 1, diff)
                if (Bernoulli == 1) {
                  m <- G
                  fi_p <- fi_p_ep
                }
            }
            i <- i + 1
            crit <- c(crit, fi_p)
            temp <- c(temp, T)
            proba <- c(proba, diff)
            T <- (c^i) * (T0)
      }
}
    return(m)
    
}



perturbationLHS=function(M,dim_xf,dim_theta) 
{  
#M un LHS
## On genère un autre LHS à partir de M, en le perturbant de façon minimale
n=nrow(M)
deb=dim_xf+1
fin=dim_xf+dim_theta
G=M
u=trunc(runif(2,0,n))  ##choix ligne
u=u+1
u1=trunc(runif(1,0,2)) ##choix bloc x/theta
u1=u1+1                ##1 ou 2
if(u1==1)
 {x=G[u[1],1:dim_xf]
  G[u[1],1:dim_xf]=G[u[2],1:dim_xf]
  G[u[2],1:dim_xf]=x}
else {x=G[u[1],deb:fin]
      G[u[1],deb:fin]=G[u[2],deb:fin]
      G[u[2],deb:fin]=x}
return(G)
} 


GRID_design=function(Xf,b_inf_theta,b_sup_theta,size_design)

 #Xf=input field data (matrix or vector)
 #b_inf_theta=lower bound of the parameter to be calibrated
 #b_sup_theta=upper bound of the parameter to be calibrated
 #design_point=number of point in the design 

{ Xf=as.matrix(Xf)
  dim_theta=length(b_inf_theta)
  N=nrow(Xf) 
  number_bloc=size_design%/%N
  rest=size_design%%N
  Xff=Xf
  if(number_bloc>=2)
   {for(i in 1:(number_bloc-1))
      {Xff=rbind(Xff,Xf)
      }}
  if(rest>0)
   {if(ncol(Xf)>1)
     {Xf=Xf[1:rest,]
      Xff=rbind(Xff,Xf)}
    else {Xf=as.matrix(Xf[1:rest,])
          Xff=rbind(Xff,Xf)}
   }
  
  Theta=c()
  for(j in 1:dim_theta)
  {M=randomLHS(size_design,1)*(b_sup_theta[j]-b_inf_theta[j])+b_inf_theta[j]
  Theta=cbind(Theta,M)}
  #print(Xff)
  #print(Theta)
  out=cbind(Xff,Theta)
  return(out)

}
  
# library(lhs)
# library(DiceDesign)
# ##EXEMPLE##
# #Xf=matrix(rep(c(0.1,0.3,0.8),2),ncol=2)
# Xf = matrix(runif(30),ncol=2)
# Xf = matrix(1:20,ncol=2)
# #b_inf_theta=c(0,0,0)
# #b_sup_theta=c(1,1,1)
# #size_design=62
# 
# #Xf=seq(0.1,0.9,by=0.1)
# #b_inf_theta=0
# #b_sup_theta=1
# #size_design=30
# 
# #out=GRID_design(Xf,b_inf_theta,b_sup_theta,size_design)
# #out_opt=maximinSA_LHS_grid(out,dim_xf=1,c=0.95,dim_theta=1,Imax=2000)
# plot(out_opt[,c(1,4)])
# table(out_opt[,c(1)],out_opt[,c(2)])
# 
# 
# 
# ## sites locations
# load("RUN_700_200bis.Rdata")
# # sites deja a l'echelle
# # 0 1 pour KX et KY a mettre a l'echelle pour runs
# b_inf_theta=c(0,0)
# b_sup_theta=c(1,1)
# size_design = 500
# 
# out = GRID_design(sites,b_inf_theta,b_sup_theta,size_design)
# out_opt=maximinSA_LHS_grid(out,dim_xf=2,c=0.95,dim_theta=2,Imax=2000)
# dim(out_opt)
# head(out_opt)
# plot(out_opt[,c(3,4)])
# 
# design_restricted = out_opt
#save(design_restricted,file="designrestricted.Rdata")
#save(design_restricted,file="designrestricted2.Rdata") # for testsing
