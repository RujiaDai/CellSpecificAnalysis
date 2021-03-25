setwd('/mnt/comp/Groups/LiuLab/User/Rdai/data/Habib2017/')
load('getx_ewce.RData')


library(EWCE)
library(sctransform)
scT = sctransform::vst(as.matrix(raw2), return_cell_attr = TRUE)
exp_scT = correct_counts(scT, as.matrix(raw2))

# Generate celltype data for just the cortex/hippocampus data
exp_CortexOnly_DROPPED = drop.uninformative.genes(exp=exp_scT,level2annot = annot2$level2class)
annotLevels = list(level1class=annot2$level1class,level2class=annot2$level2class)
fNames_CortexOnly = generate.celltype.data(exp=exp_CortexOnly_DROPPED,annotLevels=annotLevels,groupName="kiCortexOnly")

load(fNames_CortexOnly[1])
load('/mnt/Comp3/Groups/LiuLab/User/Rdai/data/allenbrain/swcam_list.RData')##ctDEGs 
reps=10000 # <- Use 100 bootstrap lists so it runs quickly, for publishable analysis use >10000
level=1 # <- Use level 1 annotations (i.e. Interneurons)


library(gridExtra)
glist<-list()
full_results1<-list()
for(i in 1:5){
full_results1[[i]] = bootstrap.enrichment.test(sct_data=ctd,hits=as.character(clist[[i]]),bg=union(gene,rownames(raw3)),
                                reps=reps,annotLevel=level,genelistSpecies="human",sctSpecies="human")
glist[[i]]<-ewce.plot(full_results1[[i]]$results,mtc_method="BH")$plain								
								
}
save(full_results1,glist,file='ewce_result_habib2017.RData')
pdf('ewce_habib2017.pdf')
grid.arrange(glist[[1]],glist[[2]],glist[[3]],glist[[4]],glist[[5]])
dev.off()						





