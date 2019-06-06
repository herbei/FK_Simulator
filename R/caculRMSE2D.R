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



# RMSE
msehom = mean(((predhom$mean-vm)/vm)^2)
msehet = mean(((predhet$mean-vm)/vm)^2)
msehetseq = mean(((predHetseq$mean-vm)/vm)^2)


boxplot(cbind(((predhom$mean-vm)/vm)^2,((predhet$mean-vm)/vm)^2,((predHetseq$mean-vm)/vm)^2),ylim=c(0,5e-4))


# scores
schom = mean(-(vm-predhom$mean)^2/(predhom$sd2) -log(predhom$sd2))
schet = mean(-(vm-predhet$mean)^2/(predhet$sd2) -log(predhet$sd2))
schetseq = mean(-(vm-predHetseq$mean)^2/(predHetseq$sd2) -log(predHetseq$sd2))

boxplot(cbind(-(vm-predhom$mean)^2/(predhom$sd2) -log(predhom$sd2),-(vm-predhet$mean)^2/(predhet$sd2) -log(predhet$sd2),-(vm-predHetseq$mean)^2/(predHetseq$sd2) -log(predHetseq$sd2)),ylim=c(-20,2))





