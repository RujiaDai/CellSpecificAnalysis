setwd('/mnt/comp/Groups/LiuLab/User/Rdai/data/deconvolution')

##bMIND deconvolution with CAM proportions
library(bMIND)
bulk<-read.csv('datExpr_AllRegressed.csv',header=T,sep=',',row.names=1)
prop<-read.csv('prop_group_m5.csv',header=T,sep=',',row.names=1)
deconv = bMIND(bulk, frac = frac, y = rbinom(n = nrow(frac), size = 1, prob = .5), ncore = 16)

save(deconv,file='braingvex_bmind_unsup.RData')




