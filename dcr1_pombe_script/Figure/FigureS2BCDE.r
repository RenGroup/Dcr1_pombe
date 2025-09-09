###################################Figure S2B ssDRIP correlation
computeMatrix reference-point --referencePoint TES -S $ssDrip_WT_dir $ssDrip_Dcr1_dir -R $conding_dir -b 500 -a 500 --skipZeros -bs 50 -p 15 -o ./ssDrip_coding_TES500_matrix.mat.gz --maxThreshold 1000 --missingDataAsZero
#R
library(ggplot2)
library(ggpubr)
enrich_matrix<-read.table("./ssDrip_coding_TES500_matrix.mat.gz",skip=1,header=F,sep="\t")
# "group_labels":["genes"],"group_boundaries":[0,4096],"sample_labels":["DG20_21IP_input","DG690_691IP_input"],"sample_boundaries":[0,20,40]}

term_gene=read.table("/work/home/path/ref/yeast/Dcr1_matched_terminated_other/dcr1_terminated_genes.bed",header=F,sep="\t")
term_gene$V2<-term_gene$V2-1
term_gene_loc=paste(term_gene[,1],term_gene[,2],term_gene[,3],sep=":")
enrich_matrix[is.na(enrich_matrix)]<-0
control_mean<-rowMeans(enrich_matrix[,7:26])
treat_mean<-rowMeans(enrich_matrix[,27:46])
plot_matrix<-data.frame(loc=paste(enrich_matrix[,1],enrich_matrix[,2],enrich_matrix[,3],sep=":"),treat_mean_CPM=as.numeric(treat_mean),control_mean_CPM=as.numeric(control_mean))
plot_matrix$type="others"
plot_matrix$type[match(term_gene_loc,plot_matrix$loc)]<-"Dcr1-terminated genes"
plot_sort<-rbind(plot_matrix[-match(term_gene_loc,plot_matrix$loc),],plot_matrix[match(term_gene_loc,plot_matrix$loc),])

#2d plot
pdf("./ssDrip_coding_TES500_2d.pdf",width=4,height=4)
ggplot(plot_matrix,aes(x=control_mean_CPM,y=treat_mean_CPM)) +
  geom_bin2d(bins=70) +
  scale_fill_continuous(type = "viridis") +
  geom_abline(slope=1)+
  labs(x = "ssDRIP WT RPM log2(IP/Input)", y = "ssDRIP dcr1Δ RPM log2(IP/Input)")+theme_bw()+
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
dev.off() 

###################################Figure S2C GC Boxplot
computeMatrix reference-point --referencePoint TES -S $GC $AT -R ${dir}dcr1_terminated_genes.bed ${dir}dcr1_expression_matched_genes.bed ${dir}pombe_coding_minus_dcr1_matched_terminated.bed -b 500 -a 500 --skipZeros -bs 50 -p 15 -o ./computeMatrix/TES500_GC_term_matrix.mat.gz
 
plotHeatmap -m ./computeMatrix/TES500_GC_term_matrix.mat.gz -out ./computeMatrix/TES500_GC_term_heatmap.pdf --colorMap RdBu_r --dpi 300 --boxAroundHeatmaps no --heatmapHeight 9 --heatmapWidth 5 --missingDataColor 1 --samplesLabel GC AT --regionsLabel term match others --perGroup
 
#boxplot
enrich_matrix<-read.table("./computeMatrix/TES500_GC_term_matrix.mat.gz",skip=1,header=F,sep="\t")
#"group_labels":["dcr1_terminated_genes.bed","dcr1_expression_matched_genes.bed","pombe_coding_minus_dcr1_matched_terminated.bed"],"group_boundaries":[0,207,585,4096],"sample_labels":["GC.sort","AT.sort"],"sample_boundaries":[0,20,40]}

cotrol_enrich_matrix<-enrich_matrix[,1:26]
cotrol_enrich_matrix[is.na(cotrol_enrich_matrix)]<-0
cotrol_enrich_mean<-data.frame(gene=cotrol_enrich_matrix[,4],CPM_sum=rowMeans(cotrol_enrich_matrix[,7:26]))
cotrol_enrich_mean$Features<-"Dcr1-terminated"
cotrol_enrich_mean$Features[208:585]<-"Expression-matched"
cotrol_enrich_mean$Features[586:4096]<-"others"

my_comparisons <- list( c("Dcr1-terminated", "Expression-matched"), c("Dcr1-terminated", "others"), c("Expression-matched", "others") )
p<-ggboxplot(cotrol_enrich_mean, x = "Features", y = "CPM_sum",
          color = "Features", palette = "jco",outlier.shape = NA,bxp.errorbar=T)+ 
  stat_compare_means(comparisons = my_comparisons)

  ggpar(p,
        legend = "right",
        xlab = "Gene type",
        ylab =  "GC content")+rotate_x_text(11)
dev.off()
write.table(cotrol_enrich_mean,"./sup_table_TESGC_Content_boxplot.tsv",quote=F,sep="\t",row.names=F)
###################################Figure S2D WT RNAPII Boxplot
computeMatrix reference-point --referencePoint TES -S ${RNAPII_dir} -R ${term_dir}dcr1_terminated_genes.bed ${term_dir}dcr1_expression_matched_genes_chr.bed ${term_dir}pombe_coding_minus_dcr1_matched_terminated.bed.bed -b 200 -a 200 --skipZeros -bs 50 -p 15 -o ./computeMatrix/RNAPII_WT2023_TES200_select_matrix.mat.gz --maxThreshold 500
 
#boxplot
enrich_matrix<-read.table("./computeMatrix/RNAPII_WT2023_TES200_select_matrix.mat.gz",skip=1,header=F,sep="\t")
# "group_labels":["dcr1_terminated_genes.bed","dcr1_expression_matched_genes_chr.bed","pombe_coding_minus_dcr1_matched_terminated.bed.bed"],"group_boundaries":[0,207,585,4096],"sample_labels":["WT_pS2_rep1IP_input"],"sample_boundaries":[0,8]}

cotrol_enrich_matrix<-enrich_matrix[,1:14]
cotrol_enrich_matrix[is.na(cotrol_enrich_matrix)]<-0
cotrol_enrich_mean<-data.frame(gene=cotrol_enrich_matrix[,4],CPM_sum=rowSums(cotrol_enrich_matrix[,7:14]))
cotrol_enrich_mean$Features<-"Dcr1-terminated"
cotrol_enrich_mean$Features[208:585]<-"Expression-matched"
cotrol_enrich_mean$Features[586:4096]<-"others"

my_comparisons <- list( c("Dcr1-terminated", "Expression-matched"), c("Dcr1-terminated", "others"), c("Expression-matched", "others") )
p<-ggboxplot(cotrol_enrich_mean, x = "Features", y = "CPM_sum",
          color = "Features", palette = "jco",outlier.shape = NA,bxp.errorbar=T)+ 
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
  
  pdf("./boxplot/RNAPII_WT_TES200_filterGene_boxplot2023_CPM.pdf",width=6,height=5)
# Visualize: Specify the comparisons you want
  ggpar(p,
        legend = "right",
        xlab = "Gene type",
        ylab = "WT RNAPII RPM log2(IP/Input)")+rotate_x_text(11)
  # Add global p-value
dev.off()
write.table(cotrol_enrich_mean,"/work/home/path/Yeast/dcr1_Figure/TES200_filterGene/RNAPII/boxplot/sup_tableRNAPII_WT_TES200_filterGene_boxplot2023.tsv",quote = F,row.names = F,sep="\t")


###################################Figure S2E RNA-Seq boxplot

library(ggplot2)
library(ggpubr)
library(dplyr)
RPKM_WT_RNA<-read.table("RPKM_WT_RNAseq.txt",header=F,sep="\t")
diff_gene<-read.table("dcr1_terminated_genes.bed",header=F,sep="\t")
match_gene<-read.table("dcr1_expression_matched_genes.bed",header=F,sep="\t")
others_gene<-read.table("pombe_coding_minus_dcr1_matched_terminated.bed",header=F,sep="\t")

RPKM_WT_RNA$feature<-"others"
RPKM_WT_RNA$feature[match(diff_gene[,4],RPKM_WT_RNA[,1])]<-"diff"
RPKM_WT_RNA$feature[match(match_gene[,4],RPKM_WT_RNA[,1])]<-"match"
RPKM_WT_RNAM<-RPKM_WT_RNA[-which(RPKM_WT_RNA$V5>400),]
RPKM_WT_RNAM<-RPKM_WT_RNAM[order(RPKM_WT_RNAM$feature),]

my_comparisons <- list( c("diff", "match"), c("diff", "others"), c("match", "others") )
p<-ggboxplot(RPKM_WT_RNAM, x = "feature", y = "V5",
          color = "feature", palette = "Dark2",outlier.shape = NA,bxp.errorbar=T)+ 
  stat_compare_means(comparisons = my_comparisons)

  ggpar(p,
        legend = "right",
        xlab = "Gene type",
        ylab = "WT RNA-seq RPKM") +rotate_x_text(11)

write.table(RPKM_WT_RNAM,"/work/home/path/Yeast/Dcr1_RNA_seq/raw_data_2014ren/bowtie2/boxplot/supplemental_table_wtRNA_RPKM_box.tsv",quote = F,row.names = F,sep="\t")