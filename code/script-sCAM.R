# swCAM

#input: 
#Xn: observation mixture matrix (sample * gene)
#Aest: estimated proporiton matrix from CAM (sample * subtype)
#Sest: estimated subtype expression matrix from CAM (subtype * gene)
#eta: eta parameter
#iteradmm: max iteration number in ADMM; larger value will cause longer running time 


Aest <-Aest
Sest <-Sest
eta <- 1000
iteradmm <- 1000
Xn<-t(Xn)


source('sCAMfastNonNeg.R')


rsCAM <- sCAMfastNonNeg(Xn, Aest, Sest, eta = eta, iteradmm=iteradmm, silent = T)

# eta parameter decides the penalty term of nuclear norm minimization.
# It can be decided by cross-validation scheme.
# According to my experience, large sample size need larger eta parameter.
# eta setting in Chuan's data:
# window 1: 500
# window 2: 350
# Window 3: 150
# window 4: 40
# window 5: 90
# Window 6: 80


Swest <- rsCAM$S #sample-wise expression
sample_names <- rownames(Xn)
gene_names <- colnames(Xn)
save(Swest, sample_names, gene_names, file = "sample-wise S.RData")