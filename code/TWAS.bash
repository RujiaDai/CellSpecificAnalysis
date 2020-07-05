##Make expression weights files
bash TWAS_weights.bash expression.matrix GENO covariant.txt

#Association analysis
for i in {1..22}; do Rscript ~/softwares/fusion_twas-master/FUSION.assoc_test.R --sumstats gwas_z.txt --weights gene_pos.txt --weights_dir ~/softwares/fusion_twas-master/ --ref_ld_chr ~/database/TWAS/LDREF/1000G.EUR. --chr $i --out scz."$i".dat; done