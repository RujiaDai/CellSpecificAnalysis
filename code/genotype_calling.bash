# Next Generation Sequencing (NGS) technology generates high throughput data in FASTQ format. The genotype calling procedure can be divided into the following steps:

## 1.	QC of FASTQ files (trim_galore, another software is FASTQC)
## https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
trim_galore SAMPLE_1.fq.gz SAMPLE_2.fq.gz -o OUTDIR --phred64 -q 20 -a GATCGGAAGAGCACACGTCT -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -e 0.1 --length 50 --stringency 1 --fastqc --paired --dont_gzip

##2.	map reads to reference genome (BWA)
##http://bio-bwa.sourceforge.net
bwa mem -M -t 2 -R '@RG\tID:97-51\tLB:97-51\tPL:ILLUMINA\tSM:97-51\tDT:2019-09-05T22:52:17' 97-51_1.fq.gz 97-51_2.fq.gz 2>./logs/97-51.log | samtools view -bS -@ 2 -o ./bamfiles/97-51.tmp.bam - 2>>./logs/97-51.log 
samtools sort -@ 2 -o ./bamfiles/97-51.sorted.bam ./bamfiles/97-51.tmp.bam 2>>./logs/97-51.log 
rm -rf ./bamfiles/97-51.tmp.bam  
samtools index ./bamfiles/97-51.sorted.bam

##3.	Mark PCR duplicates (GATK)
##https://gatk.broadinstitute.org/hc/en-us/articles/360037224712--Tool-Documentation-Index
gatk-4.1.0.0/gatk --java-options "-Xmx16g -XX:+UseParallelGC -XX:ParallelGCThreads=1" MarkDuplicates -I ./bamfiles/97-51.sorted.bam -O ./bamfiles/97-51.sorted.dedupped.bam -M ./metrics/97-51.sorted.dedupped.metrics --CREATE_INDEX true -VALIDATION_STRINGENCY LENIENT  2>>./logs/97-51.sorted.log 1>>./logs/97-51.sorted.log 
rm -f ./bamfiles/97-51.sorted.bam ./bamfiles/97-51.sorted.bam.bai  

##4.	Split'N'Trim and reassign mapping qualities (Specific for RNA-seq)
##https://gatk.broadinstitute.org/hc/en-us/articles/360035531192?id=3891
java -Xms14g -Xmx18g -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R hg19.fa -I SAMPLE.rmdup.sorted.bam -o SAMPLE.rmdup.splitNTrim.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

##5.	Base quality score recalibration (BQSR)
gatk-4.1.0.0/gatk --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=1" BaseRecalibrator -R hg38.fa -I ./bamfiles/97-51.sorted.dedupped.bam --known-sites dbsnp_138.hg38.vcf.gz --known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz  -O ./otherfiles/97-51.sorted.dedupped.recal.grp --default-base-qualities 1 2>>./logs/97-51.sorted.dedupped.log 1>>./logs/97-51.sorted.dedupped.log 
gatk-4.1.0.0/gatk --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=1" ApplyBQSR -R hg38.fa -I ./bamfiles/97-51.sorted.dedupped.bam -bqsr ./otherfiles/97-51.sorted.dedupped.recal.grp -O ./bamfiles/97-51.sorted.dedupped.recal.bam 2>>./logs/97-51.sorted.dedupped.log 1>>./logs/97-51.sorted.dedupped.log  
rm -f ./bamfiles/97-51.sorted.dedupped.bam  ./bamfiles/97-51.sorted.dedupped.bai  

##6.	Call genotypes (GATK to GVCF)
gatk-4.1.0.0/gatk --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=1" HaplotypeCaller -R hg38.fa -I ./bamfiles/97-51.sorted.dedupped.recal.bam --dbsnp dbsnp_138.hg38.vcf.gz -O ./gvcf/97-51.g.vcf.gz --emit-ref-confidence GVCF  2>>./logs/97-51.log 1>>./logs/97-51.log

##7.	GVCF to VCF
gatk --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=1" GenotypeGVCFs -R hg38.fa -V gvcf/combine.21188.g.vcf.gz -O ./vcf/chroms/95729.1.vcf.gz --dbsnp dbsnp_138.hg38.vcf.gz -L chr1  1>>./logs/GenotypeGVCFs.log  2>>./logs/GenotypeGVCFs.log

##8.	VQSR
gatk-4.1.0.0/gatk VariantFiltration -V ./vcf/allchroms.95729.vcf.gz -O ./vcf/vqsr.filter.95729.vcf.gz --filter-expression "DP < 570" --filter-name "lowDP"  1>>./logs/vqsr.log 2>>./logs/vqsr.log  
gatk-4.1.0.0/gatk --java-options "-Xmx16g -XX:+UseParallelGC -XX:ParallelGCThreads=1" VariantRecalibrator -R hg38.fa -mode SNP -V ./vcf/vqsr.filter.95729.vcf.gz -an QD -an FS -an MQRankSum -an ReadPosRankSum -an MQ -an SOR --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz --resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_138.hg38.vcf.gz  -tranches-file ./vcf/vqsr/95729-gatk_snp.tranches -rscript-file ./vcf/vqsr/95729-gatk_snp.plots.R --tranche 90.0 --tranche 93.0 --tranche 95.0 --tranche 97.0 --tranche 99.0 --tranche 99.9 --tranche 100.0 -O ./vcf/vqsr/95729-gatk_snp.recal  -an DP  1>>./logs/vqsr.log 2>>./logs/vqsr.log  
gatk-4.1.0.0/gatk --java-options "-Xmx16g -XX:+UseParallelGC -XX:ParallelGCThreads=1" ApplyVQSR -R hg38.fa -mode SNP -V ./vcf/vqsr.filter.95729.vcf.gz --recal-file ./vcf/vqsr/95729-gatk_snp.recal --tranches-file ./vcf/vqsr/95729-gatk_snp.tranches -O ./vcf/vqsr-t95-SNPrecal.95729.vcf.gz --truth-sensitivity-filter-level 95  1>>./logs/vqsr.log 2>>./logs/vqsr.log  
gatk-4.1.0.0/gatk --java-options "-Xmx16g -XX:+UseParallelGC -XX:ParallelGCThreads=1" VariantRecalibrator -R hg38.fa -mode INDEL --max-gaussians 4 -V ./vcf/vqsr-t95-SNPrecal.95729.vcf.gz -an QD -an FS -an MQRankSum -an ReadPosRankSum -an SOR --resource:mills,known=true,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -O ./vcf/vqsr/95729-gatk_indel.recal -tranches-file ./vcf/vqsr/95729-gatk_indel.tranches -rscript-file ./vcf/vqsr/95729-gatk_indel.plots.R --tranche 90.0 --tranche 93.0 --tranche 95.0 --tranche 97.0 --tranche 99.0 --tranche 99.9 --tranche 100.0  -an DP   1>>./logs/vqsr.log 2>>./logs/vqsr.log  
gatk --java-options "-Xmx16g -XX:+UseParallelGC -XX:ParallelGCThreads=1" ApplyVQSR -R hg38.fa -mode INDEL -V ./vcf/vqsr-t95-SNPrecal.95729.vcf.gz --recal-file ./vcf/vqsr/95729-gatk_indel.recal -tranches-file ./vcf/vqsr/95729-gatk_indel.tranches -O ./vcf/vqsr-t95-bothrecal.95729.vcf.gz --truth-sensitivity-filter-level 95.0  1>>./logs/vqsr.log 2>>./logs/vqsr.log

##9.	Genotype Imputation
##https://imputationserver.sph.umich.edu/index.html

##10.	Variant annotation (ANNOVAR)
##http://www.openbioinformatics.org/annovar
/opt/tools/seq-analysis/annovar201602/table_annovar.pl INFILE /zs32/data-analysis/reflib/annovar_humandb -outfile OUTFILE -buildver hg19 -protocol refGene,cytoBand,avsnp147,exac03,exac03nonpsych,exac03nontcga,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,cg69,hrcr1,kaviar_20150923,clinvar_20161128,genomicSuperDups,gerp++elem,gerp++gt2,caddgt20,cadd13gt20,revel,mcap,ljb26_all -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,f,f,f,f,f,f --remove --otherinfo --verbose --vcfinput -nastring .
