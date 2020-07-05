load('sample-specific S.RData')
v1<-log2(Si[1,,]+1)
v2<-log2(Si[2,,]+1)
v3<-log2(Si[3,,]+1)
v4<-log2(Si[4,,]+1)
v5<-log2(Si[5,,]+1)

library(WGCNA)
enableWGCNAThreads()
powers<- c(c(1:10), seq(from = 12, to=20, by=2))
sft<-pickSoftThreshold(data= t(v1), networkType = "signed", corFnc="bicor",verbose=5,powerVector=powers)
   # Power SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
# 1      1 0.493000 15.7000          0.931 8090.00  8080.000 9070.0
# 2      2 0.337000  5.9100          0.940 4450.00  4430.000 5610.0
# 3      3 0.190000  2.7500          0.960 2480.00  2450.000 3510.0
# 4      4 0.061500  1.1100          0.974 1390.00  1370.000 2250.0
# 5      5 0.000609 -0.0879          0.975  793.00   769.000 1470.0
# 6      6 0.076100 -0.8580          0.980  457.00   436.000  977.0
# 7      7 0.217000 -1.3300          0.978  267.00   249.000  662.0
# 8      8 0.412000 -1.7700          0.972  157.00   143.000  455.0
# 9      9 0.610000 -2.1300          0.978   94.00    83.200  320.0
# 10    10 0.750000 -2.4000          0.987   56.90    48.700  230.0
# 11    12 0.862000 -2.6000          0.994   21.70    17.000  123.0
# 12    14 0.911000 -2.5900          0.997    8.69     6.170   69.4
# 13    16 0.929000 -2.4900          0.998    3.68     2.300   40.5
# 14    18 0.936000 -2.3700          0.997    1.64     0.882   24.6
# 15    20 0.931000 -2.2900          0.994    0.77     0.348   15.5

						 
astro = blockwiseModules(datExpr=t(v1), networkType="signed",corType="bicor",  
                         power = 12, mergeCutHeight= 0.1, minModuleSize= 30, pamStage=FALSE, reassignThreshold=1e-6, 
                         saveTOMs=FALSE, 
                         verbose = Inf, deepSplit=4)

						
powers<- c(c(1:10), seq(from = 12, to=20, by=2))
sft<-pickSoftThreshold(data= t(v2), networkType = "signed", corFnc="bicor",verbose=5,powerVector=powers)				 
   # Power SFT.R.sq  slope truncated.R.sq  mean.k. median.k. max.k.
# 1      1  0.41600 13.500          0.813 8060.000  8060.000 8970.0
# 2      2  0.27000  5.040          0.823 4410.000  4400.000 5480.0
# 3      3  0.15800  2.420          0.872 2440.000  2430.000 3390.0
# 4      4  0.05370  1.010          0.917 1370.000  1350.000 2130.0
# 5      5  0.00139 -0.130          0.965  774.000   753.000 1370.0
# 6      6  0.05460 -0.705          0.987  443.000   424.000  891.0
# 7      7  0.19100 -1.210          0.982  256.000   241.000  590.0
# 8      8  0.41800 -1.730          0.980  150.000   138.000  397.0
# 9      9  0.59900 -2.090          0.987   88.800    79.500  275.0
# 10    10  0.71500 -2.340          0.990   53.200    46.200  194.0
# 11    12  0.82700 -2.570          0.998   19.900    15.900  101.0
# 12    14  0.87600 -2.590          0.997    7.790     5.640   54.7
# 13    16  0.90100 -2.550          0.997    3.210     2.060   31.0
# 14    18  0.91800 -2.450          0.998    1.400     0.778   18.2
# 15    20  0.92900 -2.380          0.999    0.637     0.301   11.2


nueron1 = blockwiseModules(datExpr=t(v2), networkType="signed",corType="bicor",  
                         power = 12, mergeCutHeight= 0.1, minModuleSize= 30, pamStage=FALSE, reassignThreshold=1e-6, 
                         saveTOMs=FALSE, 
                         verbose = Inf, deepSplit=4)
						 
powers<- c(c(1:10), seq(from = 12, to=20, by=2))
sft<-pickSoftThreshold(data= t(v3), networkType = "signed", corFnc="bicor",verbose=5,powerVector=powers)						 
   # Power SFT.R.sq  slope truncated.R.sq  mean.k. median.k.  max.k.
# 1      1  0.37300 13.900          0.922 8090.000   8100.00 9010.00
# 2      2  0.24400  5.260          0.937 4450.000   4450.00 5530.00
# 3      3  0.11600  2.250          0.948 2480.000   2470.00 3430.00
# 4      4  0.02900  0.800          0.962 1390.000   1380.00 2150.00
# 5      5  0.00111 -0.125          0.976  790.000    776.00 1370.00
# 6      6  0.07330 -0.898          0.984  454.000    440.00  891.00
# 7      7  0.21300 -1.420          0.985  263.000    251.00  588.00
# 8      8  0.37100 -1.880          0.986  155.000    145.00  398.00
# 9      9  0.51300 -2.290          0.990   91.800     84.10  275.00
# 10    10  0.62400 -2.550          0.994   55.100     49.30  193.00
# 11    12  0.75500 -2.910          0.995   20.600     17.30   98.20
# 12    14  0.82400 -3.010          0.998    8.030      6.25   52.30
# 13    16  0.85500 -3.050          0.998    3.280      2.33   29.00
# 14    18  0.87900 -2.920          0.995    1.410      0.89   16.60
# 15    20  0.88900 -2.810          0.987    0.633      0.35    9.72
nueron2 = blockwiseModules(datExpr=t(v3), networkType="signed",corType="bicor",  
                         power = 14, mergeCutHeight= 0.1, minModuleSize= 30, pamStage=FALSE, reassignThreshold=1e-6, 
                         saveTOMs=FALSE, 
                         verbose = Inf, deepSplit=4)
						 
powers<- c(c(1:10), seq(from = 12, to=20, by=2))
sft<-pickSoftThreshold(data= t(v4), networkType = "signed", corFnc="bicor",verbose=5,powerVector=powers)
   # Power SFT.R.sq  slope truncated.R.sq  mean.k. median.k.  max.k.
# 1      1   0.4150 12.000          0.835 8110.000  8120.000 9120.00
# 2      2   0.2400  4.360          0.845 4490.000  4490.000 5670.00
# 3      3   0.0963  1.800          0.883 2510.000  2510.000 3590.00
# 4      4   0.0115  0.470          0.920 1430.000  1420.000 2300.00
# 5      5   0.0135 -0.425          0.950  822.000   808.000 1500.00
# 6      6   0.0993 -1.040          0.972  479.000   466.000  988.00
# 7      7   0.2150 -1.470          0.981  282.000   271.000  660.00
# 8      8   0.3480 -1.810          0.986  169.000   159.000  446.00
# 9      9   0.4670 -2.050          0.989  102.000    94.300  305.00
# 10    10   0.5680 -2.230          0.992   62.300    56.500  210.00
# 11    12   0.7210 -2.440          0.998   24.100    20.800  103.00
# 12    14   0.8020 -2.700          0.999    9.770     7.890   54.80
# 13    16   0.8520 -2.770          0.998    4.130     3.080   30.20
# 14    18   0.8730 -2.810          0.998    1.820     1.260   17.10
# 15    20   0.8820 -2.820          0.997    0.839     0.524    9.99
micro = blockwiseModules(datExpr=t(v4), networkType="signed",corType="bicor",  
                         power = 14, mergeCutHeight= 0.1, minModuleSize= 30, pamStage=FALSE, reassignThreshold=1e-6, 
                         saveTOMs=FALSE, 
                         verbose = Inf, deepSplit=4)

sft<-pickSoftThreshold(data= t(v5), networkType = "signed", corFnc="bicor",verbose=5,powerVector=powers)
   # Power SFT.R.sq slope truncated.R.sq  mean.k. median.k.  max.k.
# 1      1  0.36800 12.30          0.870 8060.000  8070.000 9000.00
# 2      2  0.18300  4.14          0.884 4420.000  4420.000 5550.00
# 3      3  0.03570  1.21          0.894 2450.000  2440.000 3480.00
# 4      4  0.00247 -0.25          0.928 1370.000  1360.000 2210.00
# 5      5  0.06850 -1.14          0.962  777.000   761.000 1430.00
# 6      6  0.19600 -1.76          0.983  445.000   430.000  940.00
# 7      7  0.34600 -2.27          0.988  257.000   245.000  625.00
# 8      8  0.47400 -2.55          0.991  150.000   141.000  421.00
# 9      9  0.59300 -2.81          0.994   88.700    81.400  288.00
# 10    10  0.68600 -2.96          0.993   52.900    47.500  199.00
# 11    12  0.79400 -3.07          0.996   19.500    16.500   97.80
# 12    14  0.85300 -3.02          0.995    7.480     5.890   50.10
# 13    16  0.89000 -2.91          0.993    3.000     2.170   26.60
# 14    18  0.91400 -2.78          0.995    1.250     0.822   14.50
# 15    20  0.92600 -2.64          0.994    0.549     0.319    8.16

oligo = blockwiseModules(datExpr=t(v5), networkType="signed",corType="bicor",  
                         power = 14, mergeCutHeight= 0.1, minModuleSize= 30, pamStage=FALSE, reassignThreshold=1e-6, 
                         saveTOMs=FALSE, 
                         verbose = Inf, deepSplit=4)

astro<-astro
neuron1<-nueron1
neuron2<-nueron2 
network<-list(astro,neuron1,neuron2,micro,oligo)

power<-c(12,12,14,14,14)
adjset<-list()
tomset<-list()
for(i in 1:5){
adjset[[i]]<-adjacency(datExpr=t(datExpr[[i]]),type='signed',power=power[[i]],corFnc='bicor')
tomset[[i]]<-TOMsimilarity(adjset[[i]], TOMType = "signed")
print(i)
}
consensusTOM = pmin(TOM[1, , ], TOM[2, , ]);
consTree = hclust(as.dist(1-consensusTOM), method = "average")
minModuleSize = 30;

co1<-network[[1]]$colors
co2<-network[[2]]$colors
co3<-network[[3]]$colors
co4<-network[[4]]$colors
co5<-network[[5]]$colors
pdf('sscam_dendroplot_merged_bulk.pdf')
plotDendroAndColors(consTree, cbind(bco,co1,co2,co3,co4,co5),
c('bulk','astrocyte','neuron_1','neuron_2','microglia','oligodendrocyte'),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

##trait association
for(i in 1:length(network)){
+     module<-network[[i]]$colors
+ MEs = network[[i]]$MEs
+     netname<-names(network)[i]
+ modTrait = data.frame()
+ 
+ for(i in 1:ncol(MEs)) {
+   me = MEs[,i]
+   moduleColor = gsub("ME", "", colnames(MEs)[i])
+ 
+ 
+   s = summary(lm(me ~ Group, data = datM))$coefficients
+   for(grp in c("BD", "SCZ")) {
+ rowID = paste0("Group", grp)
+ modTrait = rbind(modTrait, 
+    data.frame(Module=moduleColor, Group=grp, 
+   beta = s[rowID, "Estimate"], SE = s[rowID, "Std. Error"], t=s[rowID, "t value"], p=s[rowID, "Pr(>|t|)"]))
+   }
+ }
+ 
+ modTrait$fdr = p.adjust(modTrait$p, method = "fdr")
+ modTrait$signedLog10fdr = -log10(modTrait$fdr) * sign(modTrait$beta)
+ modTrait$signedLog10fdr[modTrait$fdr > .05] = 0
+ modTrait$text = signif(modTrait$beta, 1)
+ modTrait$text[modTrait$fdr > 0.05] = ""
+ modTrait$Label = ""
+ modTrait$Label[modTrait$fdr<.05] = "*"
+     
+     write.csv(modTrait,file=paste0(netname,'modTrait.csv'))
+ 
+ }
kME<-list()
for(i in 1:length(network)){
   MEs = network[[i]]$MEs
   kMEtable[[i]]=as.data.frame(cor(t(datExpr[[i]]),MEs,use="p"))
}




##GWAS enrichment
d = dir("/zs32/home/rjdai/PsychENCODE/combined/12/MAGMA/",pattern=".genes.out")

modGWAS = data.frame()
for(i in 1:length(d)) {
  gwas = read.delim(paste0("/zs32/home/rjdai/PsychENCODE/combined/12/MAGMA/", d[[i]]), sep="")
  gwas = gwas[match(rownames(kMEtable[[1]]), gwas$GENE),]
  gwasName = gsub(".genes.out", "", d[[i]])
  gwas$geneLength = gwas$STOP - gwas$START
  gwas$logGeneLen = log10(gwas$geneLength)
  gwas$logP = -log10(gwas$P)
  
  for(j in 1:ncol(kMEtable[[1]])){
    moduleColor  = gsub("ME", "", colnames(kMEtable)[j])
    #Spearman's Correlation between module membership and GWAS significance
    cor = cor.test(kMEtable[,j], gwas$ZSTAT, method="spearman", use="pairwise.complete.obs")
    modGWAS = rbind(modGWAS, 
                    data.frame(Module=moduleColor, 
                               GWAS=gwasName, rho=as.numeric(cor$estimate), p=as.numeric(cor$p.value))}
}

modGWAS$fdr = p.adjust(modGWAS$p, method = "fdr")
 modGWAS$fdr[modGWAS$rho<0] =1 
modGWAS$log10fdr = -log10(modGWAS$fdr) 
modGWAS$text = signif(modGWAS$log10fdr,2)
modGWAS$text[modGWAS$fdr > 0.1] = ""

#hub-network
modules<-network[[2]]$colors
for(i in 1:ncol(kMEtable[[2]])) {
  moduleColor = gsub("ME", "", colnames(kMEtable[[2]])[i])
  moduleGenes = annot$gene_id[modules==moduleColor]
  
  hubGenes = moduleGenes[order(kMEtable[[2]][moduleGenes,i], decreasing = T)[1:20]]
  hubGene.symbols = annot$gene_name[match(hubGenes, annot$gene_id)]
  adjMat = adjacency(t(datExpr[[2]][hubGenes,]),type = "signed",corFnc = "bicor", power=12)
  adjMat[adjMat < quantile(adjMat,0.1)]=0
  graph <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=T,diag=F)
  plotcord= data.frame(layout_with_fr(graph))
  colnames(plotcord) = c("X1","X2")
  edgelist <- get.edgelist(graph,names = F)
  edges <- data.frame(plotcord[edgelist[,1],], plotcord[edgelist[,2],])
  colnames(edges) <- c("X1","Y1","X2","Y2")
  plotcord = cbind(plotcord, data.frame(gene=hubGene.symbols))
  ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = 0.5, colour="grey") + 
    geom_point(aes(X1, X2), data=plotcord,color=moduleColor,size=4) + geom_text(aes(x=X1, y=X2+.2, label=gene),data=plotcord) +
    theme_classic() + labs(x="", y="") + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  ggsave( file = paste0("./nueron1", moduleColor, ".png"), width = 12, height = 6, type = "cairo", dpi = 2000)
}

##cell type enrichment
zhang.datExpr = read.csv("/zs32/home/rjdai/PsychENCODE/combined/12/enrichment/datExpr.zhangHuman.avgForPSI.GSE21653.csv",row.names=1)
pSI.zhang = specificity.index(pSI.in=zhang.datExpr,bts=100,p_max=.1, e_min=0.3);


modCellType = data.frame()
modules<-unique(network[[1]]$colors)
for(i in 1:length(unique(modules))) {
  print(i)
  me<-modules[[i]]
  moduleColor =network[[1]]$colors
  genesInMod = as.character(annot$gene_name[me==moduleColor])
  
  f.zhang = fisher.iteration(pSI.zhang, genesInMod,p.adjust = F)

  modCellType = rbind(modCellType, 
                      data.frame(Dataset="Zhang", Module=me,  CellType=rownames(f.zhang), p=f.zhang[,1]))
                      
}


modCellType$fdr = p.adjust(modCellType$p, method = "fdr")
modCellType$log10fdr = -log10(modCellType$fdr) 
modCellType$text = signif(modCellType$log10fdr,2)
modCellType$text[modCellType$fdr > 0.05] = ""
modCellType$CellType = gsub("Oligodendrocyte", "Oligo", modCellType$CellType)
write.csv(modCellType,'astromodCellType.csv')


##Go term
me1<-network[[1]]$MEs
me2<-network[[2]]$MEs
me3<-network[[3]]$MEs
me4<-network[[4]]$MEs
me5<-network[[5]]$MEs

names(me1)<-gsub('ME','astro ',names(me1))
names(me2)<-gsub('ME','neuron1 ',names(me2))
names(me3)<-gsub('ME','neuron2 ',names(me3))
names(me4)<-gsub('ME','micro ',names(me4))
names(me5)<-gsub('ME','oligo ',names(me5))


for(i in 1:length(oligo)){

name<-oligo[i]
ggene<-df[,1][df$V5==name]
go = gprofiler(query=as.vector(ggene),  correction_method = "fdr",custom_bg = as.vector(df[,1]),ordered_query = F)
go = go[order(go$p.value)[1:min(10,nrow(go))],]
ggplot(go, aes(x=reorder(term.name, -log10(p.value)), y=-log10(p.value))) + geom_bar(stat="identity", fill="royalblue") +theme_light(base_size=15)+ coord_flip() + xlab("") + geom_hline(yintercept=-log10(0.05), lty=2, color="red")

ggsave( file = paste0("./oligo",name, ".png"), width = 10, height = 6, type = "cairo", dpi = 2000)
}


##top20 network
load('ssCAM_network.RData')
annot<-read.csv('annot.csv',header=T,sep=",")
datMeta<-read.csv('datMeta_QC.csv',header=T,sep=",")
kMEtable<-read.csv('BGX_COE_kME.csv',header=T,sep=",")

mkme<-subset(kMEtable,V2=='violet')
moduleColor<-'violet'
hubGenes = mkme$X[order(mkme$kme2.value,decreasing=T)][1:20]
hubGene.symbols = annot$gene_id[match(hubGenes, annot$gene_name)]
adjMat = adjacency(t(datExpr[[2]][hubGene.symbols,]),type = "signed",corFnc = "bicor", power=10)
adjMat[!is.na(adjMat)]<-0.01
graph <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=T,diag=F)
plotcord= data.frame(layout_with_fr(graph))
colnames(plotcord) = c("X1","X2")
edgelist <- get.edgelist(graph,names = F)
edges <- data.frame(plotcord[edgelist[,1],], plotcord[edgelist[,2],])
colnames(edges) <- c("X1","Y1","X2","Y2")
plotcord = cbind(plotcord, data.frame(gene=hubGenes))
ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = 0.5, colour="grey") + 
  geom_point(aes(X1, X2), data=plotcord,color=moduleColor,size=3) + geom_text(aes(x=X1, y=X2+.2, label=gene),data=plotcord,size=5) +
  theme_classic() + labs(x="", y="") + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),legend.position='none')


##trait and GWAS plot
library(gridExtra)
setwd('G:/ssCAM/COE/trait/')
files<-list.files(getwd(),full.name=TRUE)
trait<-lapply(files,function(x){read.csv(x,header=T,sep=",")})
setwd('G:/ssCAM/COE/GWAS/')
files<-list.files(getwd(),full.name=TRUE)
files<-files[c(2,6,10,14,18)]
gwas<-lapply(files,function(x){read.csv(x,header=T,sep=",")})

modTrait<-subset(trait[[5]],Module=='cyan')
modGWAS<-subset(gwas[[5]],Module=='cyan')
g1<-ggplot(modTrait,aes(x=Group,y=beta, Label=Label, fill=Group)) + geom_bar(stat="identity", position=position_dodge(width=.75),width=.75) +
  geom_errorbar(aes(x=Group,ymin=beta-SE, ymax=beta+SE),width=.25,position=position_dodge(width = .75))  + 
  geom_text(aes(label=Label, y=beta + sign(beta)*SE + sign(beta)*.0005), position=position_dodge(width=.75),color="red",size=8) + 
  theme_classic(base_size=15)+
  theme(axis.text.x = element_text(angle=-60, hjust=0),legend.position='none') + labs(y="log2FC", x="")  
  
g2<-ggplot(modGWAS, aes(x=reorder(GWAS, log10fdr),y=log10fdr)) +geom_bar(stat='identity',position='dodge')+ coord_flip()+theme_classic(base_size=15)+labs(x='GWAS')+ geom_hline(yintercept=-log10(0.05), lty=2, color="red")
grid.arrange(g1,g2)


##go plot
setwd('F:/oligo/')
kme<-read.csv('BGX_COE_kME.csv',header=T,sep=",")
library(gProfileR)
kme2<-kme[order(kme$kme2.value,decreasing=T),]
go = gprofiler(query=as.vector(subset(kme2,V2=='violet')), max_set_size = 1000, correction_method = "fdr",hier_filtering = "strong", custom_bg = as.vector(kme2[,1]), src_filter = c("GO", "KEGG"),ordered_query = T)
go2 = go[order(go$p.value)[1:min(10,nrow(go))],]
ggplot(go2, aes(x=reorder(term.name, -log10(p.value)), y=-log10(p.value))) + geom_bar(stat="identity", fill="royalblue") + coord_flip() + xlab("") + geom_hline(yintercept=-log10(0.05), lty=2, color="red")+theme(axis.text.y = element_text(size=15))
  