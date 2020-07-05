##cell network preservation test
load('ssCAM_networks.RData')

colors<-lapply(network,function(x){x$colors})
datExpr1<-t(datExpr[[1]])
colors1<-colors[[1]]

datExpr2<-t(datExpr[[3]])
colors2<-colors[[3]]



multiExpr = list(astrocyte = list(data = datExpr1), upneuron = list(data = datExpr2));
multiColor = list(astrocyte = colors1, upneuron = colors2);
nSets = 2

system.time( {
mp = modulePreservation(multiExpr, multiColor,
referenceNetworks = c(1:2),
networkType = "signed",
nPermutations = 100,
randomSeed = 1,
verbose = 3)
} );


save(mp,file='mp13.RData')



datExpr1<-t(datExpr[[1]])
colors1<-colors[[1]]

datExpr2<-t(datExpr[[4]])
colors2<-colors[[4]]



multiExpr = list(astrocyte = list(data = datExpr1), microglia = list(data = datExpr2));
multiColor = list(astrocyte = colors1, microglia = colors2);
nSets = 2

system.time( {
mp = modulePreservation(multiExpr, multiColor,
referenceNetworks = c(1:2),
networkType = "signed",
nPermutations = 100,
randomSeed = 1,
verbose = 3)
} );


save(mp,file='mp14.RData')




datExpr1<-t(datExpr[[1]])
colors1<-colors[[1]]

datExpr2<-t(datExpr[[5]])
colors2<-colors[[5]]



multiExpr = list(astrocyte = list(data = datExpr1), oligodendrocyte = list(data = datExpr2));
multiColor = list(astrocyte = colors1, oligodendrocyte = colors2);
nSets = 2

system.time( {
mp = modulePreservation(multiExpr, multiColor,
referenceNetworks = c(1:2),
networkType = "signed",
nPermutations = 100,
randomSeed = 1,
verbose = 3)
} );


save(mp,file='mp15.RData')


datExpr1<-t(datExpr[[2]])
colors1<-colors[[2]]

datExpr2<-t(datExpr[[4]])
colors2<-colors[[4]]



multiExpr = list(deepneuron = list(data = datExpr1), microglia = list(data = datExpr2));
multiColor = list(deepneuron = colors1, microglia = colors2);
nSets = 2

system.time( {
mp = modulePreservation(multiExpr, multiColor,
referenceNetworks = c(1:2),
networkType = "signed",
nPermutations = 100,
randomSeed = 1,
verbose = 3)
} );


save(mp,file='mp24.RData')



datExpr1<-t(datExpr[[2]])
colors1<-colors[[2]]

datExpr2<-t(datExpr[[5]])
colors2<-colors[[5]]



multiExpr = list(deepneuron = list(data = datExpr1), oligodendrocyte= list(data = datExpr2));
multiColor = list(deepneuron = colors1, oligodendrocyte = colors2);
nSets = 2

system.time( {
mp = modulePreservation(multiExpr, multiColor,
referenceNetworks = c(1:2),
networkType = "signed",
nPermutations = 100,
randomSeed = 1,
verbose = 3)
} );


save(mp,file='mp25.RData')


datExpr1<-t(datExpr[[3]])
colors1<-colors[[3]]

datExpr2<-t(datExpr[[4]])
colors2<-colors[[4]]



multiExpr = list(upneuron = list(data = datExpr1), microglia= list(data = datExpr2));
multiColor = list(upneuron = colors1, microglia = colors2);
nSets = 2

system.time( {
mp = modulePreservation(multiExpr, multiColor,
referenceNetworks = c(1:2),
networkType = "signed",
nPermutations = 100,
randomSeed = 1,
verbose = 3)
} );


save(mp,file='mp34.RData')



datExpr1<-t(datExpr[[3]])
colors1<-colors[[3]]

datExpr2<-t(datExpr[[5]])
colors2<-colors[[5]]



multiExpr = list(upneuron = list(data = datExpr1), oligodendrocyte= list(data = datExpr2));
multiColor = list(upneuron = colors1, oligodendrocyte = colors2);
nSets = 2

system.time( {
mp = modulePreservation(multiExpr, multiColor,
referenceNetworks = c(1:2),
networkType = "signed",
nPermutations = 100,
randomSeed = 1,
verbose = 3)
} );


save(mp,file='mp35.RData')


datExpr1<-t(datExpr[[4]])
colors1<-colors[[4]]

datExpr2<-t(datExpr[[5]])
colors2<-colors[[5]]



multiExpr = list(microglia = list(data = datExpr1), oligodendrocyte= list(data = datExpr2));
multiColor = list(microglia = colors1, oligodendrocyte = colors2);
nSets = 2

system.time( {
mp = modulePreservation(multiExpr, multiColor,
referenceNetworks = c(1:2),
networkType = "signed",
nPermutations = 100,
randomSeed = 1,
verbose = 3)
} );


save(mp,file='mp45.RData')



datExpr1<-t(ssdata[[1]])
colors1<-sscolors[[1]]

datExpr2<-t(repdata[[2]])
colors2<-repcolors[[2]]



multiExpr = list(bgx = list(data = datExpr1),rep= list(data = datExpr2));
multiColor = list(bgx = colors1, rep = colors2);
nSets = 2

system.time( {
mp = modulePreservation(multiExpr, multiColor,
referenceNetworks = 1,
networkType = "signed",
nPermutations = 100,
randomSeed = 1,
verbose = 3)
} );


save(mp,file='mp_rep_v1.RData')


datExpr1<-t(ssdata[[2]])
colors1<-sscolors[[2]]

datExpr2<-t(repdata[[1]])
colors2<-repcolors[[1]]



multiExpr = list(bgx = list(data = datExpr1),rep= list(data = datExpr2));
multiColor = list(bgx = colors1, rep = colors2);
nSets = 2

system.time( {
mp = modulePreservation(multiExpr, multiColor,
referenceNetworks = 1,
networkType = "signed",
nPermutations = 100,
randomSeed = 1,
verbose = 3)
} );


save(mp,file='mp_rep_v2.RData')


datExpr1<-t(ssdata[[3]])
colors1<-sscolors[[3]]

datExpr2<-t(repdata[[3]])
colors2<-repcolors[[3]]



multiExpr = list(bgx = list(data = datExpr1),rep= list(data = datExpr2));
multiColor = list(bgx = colors1, rep = colors2);
nSets = 2

system.time( {
mp = modulePreservation(multiExpr, multiColor,
referenceNetworks = 1,
networkType = "signed",
nPermutations = 100,
randomSeed = 1,
verbose = 3)
} );


save(mp,file='mp_rep_v3.RData')


datExpr1<-t(ssdata[[4]])
colors1<-sscolors[[4]]

datExpr2<-t(repdata[[4]])
colors2<-repcolors[[4]]



multiExpr = list(bgx = list(data = datExpr1),rep= list(data = datExpr2));
multiColor = list(bgx = colors1, rep = colors2);
nSets = 2

system.time( {
mp = modulePreservation(multiExpr, multiColor,
referenceNetworks = 1,
networkType = "signed",
nPermutations = 100,
randomSeed = 1,
verbose = 3)
} );


save(mp,file='mp_rep_v4.RData')


datExpr1<-t(ssdata[[5]])
colors1<-sscolors[[5]]

datExpr2<-t(repdata[[5]])
colors2<-repcolors[[5]]



multiExpr = list(bgx = list(data = datExpr1),rep= list(data = datExpr2));
multiColor = list(bgx = colors1, rep = colors2);
nSets = 2

system.time( {
mp = modulePreservation(multiExpr, multiColor,
referenceNetworks = c(1,2),
networkType = "signed",
nPermutations = 100,
randomSeed = 1,
verbose = 3)
} );


save(mp,file='mp_rep_v5.RData')




pdf('GSE12649_bulk_cell_preservation.pdf',width=12,height=6)


ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];















pdf('testunisgned_refcsu_preservation.pdf',width=12,height=6)
plotMods = !(modColors %in% c("grey", "gold"));
text = "";
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
mains = c("Preservation Median rank", "Preservation Zsummary");
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
min = min(plotData[, p], na.rm = TRUE);
max = max(plotData[, p], na.rm = TRUE);
# Adjust ploting ranges appropriately
if (p==2)
{
if (min > -max/10)
min = -max/10
ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
}
else
ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
main = mains[p],
cex = 2.4,
ylab = mains[p], xlab = "Module size", log = "x",
ylim = ylim,
xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.2)
labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs =
0.08);
if (p==2)
{
abline(h=0)
abline(h=2, col = "blue", lty = 2)
abline(h=10, col = "darkgreen", lty = 2)
}
}


dev.off()





pdf('GSE12649_neuron1_preservation.pdf',width=12,height=6)
to_plot=TRUE
testid<-c(1,3,4,5,6)
names<-c('bulk','astro','neuron2','micro','oligo')
for(i in 1:5){
ref = 2
test = testid[[i]]
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];



plotMods = !(modColors %in% c("grey", "gold"));
text = "";
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
mains = c("Preservation Median rank", "Preservation Zsummary");
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
min = min(plotData[, p], na.rm = TRUE);
max = max(plotData[, p], na.rm = TRUE);
# Adjust ploting ranges appropriately
if (p==2)
{
if (min > -max/10)
min = -max/10
ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
}
else
ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
main = mains[p],
cex = 2.4,
ylab = mains[p], xlab = "Module size", log = "x",
ylim = ylim,
xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.2)
labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs =
0.08);
if (p==2)
{
abline(h=0)
abline(h=2, col = "blue", lty = 2)
abline(h=10, col = "darkgreen", lty = 2)
}
}

write.csv(( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) ),file=paste0('neuron1_zsummary',names[[i]],'.csv'))
}
dev.off()

p=2
plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
main = mains[p],
cex = 2.4,
ylab = mains[p], xlab = "Module size", log = "x",
ylim = ylim,
xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.2)
labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs =
0.08);
if (p==2)
{
abline(h=0)
abline(h=2, col = "blue", lty = 2)
abline(h=10, col = "darkgreen", lty = 2)
}
