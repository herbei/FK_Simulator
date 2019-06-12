load("testdesign2D.Rdata")


summary(vn)
summary(vs)
summary(vm)



# prediction
load("Designseq2D.Rdata")
load("Designseq2Dbisafterseq.Rdata")
load("Designseq2Dbisafterseq500.Rdata")
library(hetGP)
load("designseq2D10paths2.Rdata")


#Hom GP
predhom = predict(x = testdesign, object = Ghom)

#Het GP
predhet = predict(x = testdesign, object = Ghet)

#Het GP seq
predHetseq = predict(x = testdesign, object = mod)


vmNorm = (vm-Mnorm)/SDnorm

# RMSE
msehom = mean(((predhom$mean-vmNorm)^2))
msehet = mean(((predhet$mean-vmNorm)^2))
msehetseq = mean(((predHetseq$mean-vmNorm)^2))

dfmse = as.data.frame(cbind((predhom$mean-vmNorm)^2,(predhet$mean-vmNorm)^2,(predHetseq$mean-vmNorm)^2))
names(dfmse) = c("hom","het","seqhet")
boxplot(dfmse,main="MSE",ylim=c(0,.2))

which(dfmse$seqhet>1000)
dfmse[383,]
vmNorm[383]


# scores
schom = mean(-(vmNorm-predhom$mean)^2/(predhom$sd2) -log(predhom$sd2))
schet = mean(-(vmNorm-predhet$mean)^2/(predhet$sd2) -log(predhet$sd2))
schetseq = mean(-(vmNorm-predHetseq$mean)^2/(predHetseq$sd2) -log(predHetseq$sd2))


dfscore = as.data.frame(cbind(-(vm-predhom$mean)^2/(predhom$sd2) -log(predhom$sd2),-(vm-predhet$mean)^2/(predhet$sd2) -log(predhet$sd2),-(vm-predHetseq$mean)^2/(predHetseq$sd2) -log(predHetseq$sd2)))
names(dfscore) = c("hom","het","seqhet")
boxplot(dfscore,ylim=c(-20,2),main="scores")





