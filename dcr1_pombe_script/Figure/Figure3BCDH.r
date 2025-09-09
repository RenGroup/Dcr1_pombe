###################################figure3B
computeMatrix reference-point --referencePoint TES -S ./DNAPd_DNAPe_compare.bw -R ${term_dir}dcr1_terminated_genes.bed ${term_dir}dcr1_expression_matched_genes.bed ${term_dir}pombe_coding_minus_dcr1_matched_terminated.bed -b 500 -a 500 --skipZeros -bs 10 -p 10 -o ./computeMatrix/DNAPD_DNAPE_enrich_compare_TES.mat.gz --maxThreshold 500 --missingDataAsZero

library(matrixStats)
enrich_matrix<-read.table("./computeMatrix/DNAPD_DNAPE_enrich_compare_TES_matrix_250204.mat.gz",skip=1,header=F,sep="\t")
# "group_labels":["dcr1_terminated_genes.bed","dcr1_expression_matched_genes.bed","pombe_coding_minus_dcr1_matched_terminated.bed"],"group_boundaries":[0,207,585,4096],"sample_labels":["DNAPd_DNAPe_compare"],"sample_boundaries":[0,100]}

cotrol_enrich_matrix<-enrich_matrix[,c(1:106)]
cotrol_enrich_matrix[is.na(cotrol_enrich_matrix)]<-0
 
cotrol_enrich_mean<-data.frame(gene=cotrol_enrich_matrix[,4],CPM_sum=do.call(pmax, cotrol_enrich_matrix[,7:106]))
 
cotrol_enrich_mean$Features<-"Dcr1-terminated"
cotrol_enrich_mean$Features[208:585]<-"Expression-matched"
cotrol_enrich_mean$Features[586:4096]<-"others"

cols <- c("#F76D5E", "#FFFFBF", "#72D8FF")
#  density
p<-ggplot(cotrol_enrich_mean, aes(x =CPM_sum , fill = Features)) +
  geom_density(alpha = 0.3) + 
  scale_fill_manual(values = cols)+
  guides(fill = guide_legend(title = "Class"))+xlab("DNAP usage prefference compare")+coord_cartesian(xlim =  c(-0.5, 3))

ggpar(p,
        legend = "right")

# 单独统计other和term差异
p.value<-wilcox.test(subset(cotrol_enrich_mean,Features=="Dcr1-terminated")[,2],subset(cotrol_enrich_mean,Features=="others")[,2])$p.value
p.adjust(p.value) 


###################################figure3C
computeMatrix reference-point --referencePoint TES  -S ${DNAPd_DNAPe_compare} -R dcr1_terminated_CD.bed dcr1_terminated_HO.bed -b 1000 -a 1000 --skipZeros -bs 50 -p 15 -o ./computeMatrix/TES_DNAP_compare_matrix.mat.gz --maxThreshold 500 --missingDataAsZero

plotProfile -m ./computeMatrix/TES_DNAP_compare_matrix.mat.gz -out ./computeMatrix/TES_DNAP_compare_plot.pdf --plotFileFormat pdf --outFileNameData ./computeMatrix/TES_DNAP_compare_myProfile.tab --regionsLabel term_CD term_HO --samplesLabel DNAPd_DNAPe_compare

##smooth
library(ggplot2)
library(tidyr)
library(cowplot)
library(dplyr)
library(reshape2)
library(ggpubr)

tab.data<-read.table("./computeMatrix/TES_DNAP_compare_myProfile.tab",header = F,sep="\t")
tab.data<-t(tab.data)
dim(tab.data)
# [1] 42  4
tab.rowname<-tab.data[3:42,2]
tab.data.term<-as.data.frame(tab.data[3:42,c(2,3)])
tab.data.match<-as.data.frame(tab.data[3:42,c(2,4)])

tab.data.term$type<-c("Dcr1-terminated Co-direction genes")
tab.data.match$type<-c("Dcr1-terminated Head-on genes")

rownames(tab.data.term)<-tab.rowname
rownames(tab.data.match)<-tab.rowname

colnames(tab.data.term)<-c("loc","Value","type")
colnames(tab.data.match)<-c("loc","Value","type")

all_data<-rbind(tab.data.term,tab.data.match)

pdf("./computeMatrix/smooth_DNAPcompare_TES.pdf",width=6,height=3)
ggplot(all_data,aes(as.numeric(loc),as.numeric(Value),colour=type))+ 
  geom_smooth(span = 0.3, se = F)+
  scale_x_continuous(breaks=c(1,20,40), labels =c("-1000bp","TES","1000bp"))+
  labs(x = "Position", y = "DNAPD/DNAPE preference")+theme_bw()+
  theme(panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank(),  
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"))
  dev.off() #()


###################################figure3D
computeMatrix reference-point --referencePoint TES -S DNAPD_DNAPE_pref_90.bw -R ${dir}dcr1_terminated_genes.bed ${dir}dcr1_expression_matched_genes.bed ${dir}pombe_coding_minus_dcr1_matched_terminated.bed  -b 500 -a 500 --skipZeros -bs 50 -p 15 -o ./computeMatrix/TES500_DNAPD_DNAPE_90MIN_gene_matrix.mat.gz --maxThreshold 500 --missingDataAsZero 

##R
library(matrixStats)
enrich_matrix<-read.table("./computeMatrix/TES500_DNAPD_DNAPE_90MIN_gene_matrix.mat.gz",skip=1,header=F,sep="\t")
# "group_labels":["dcr1_terminated_genes.bed","dcr1_expression_matched_genes.bed","pombe_coding_minus_dcr1_matched_terminated.bed"],"group_boundaries":[0,207,585,4096],"sample_labels":["DNAPD_DNAPE_pref_90"],"sample_boundaries":[0,20]}

cotrol_enrich_matrix<-enrich_matrix[,c(1:26)]
cotrol_enrich_matrix[is.na(cotrol_enrich_matrix)]<-0
cotrol_enrich_mean<-data.frame(gene=cotrol_enrich_matrix[,4],CPM_sum=do.call(pmax, cotrol_enrich_matrix[,7:26]))

cotrol_enrich_mean$Features<-"Dcr1-terminated"
cotrol_enrich_mean$Features[205:585]<-"Expression-matched"
cotrol_enrich_mean$Features[586:4096]<-"others"

cols <- c("#F76D5E", "#FFFFBF", "#72D8FF")
p<-ggplot(cotrol_enrich_mean, aes(x =CPM_sum , fill = Features)) +
  geom_density(alpha = 0.3) + 
  scale_fill_manual(values = cols)+
  guides(fill = guide_legend(title = "Class"))+xlab("DNAP usage prefference compare")+coord_cartesian(xlim =  c(-0.5, 3))
ggpar(p,
        legend = "right")

p.value<-wilcox.test(subset(cotrol_enrich_mean,Features=="Dcr1-terminated")[,2],subset(cotrol_enrich_mean,Features=="others")[,2])$p.value
p.adjust(p.value) 


#################################figure3H
computeMatrix reference-point --referencePoint TES -S Rad52_S_exp_dcr1_WT.bw -R ${dir}dcr1_terminated_genes.bed ${dir}dcr1_expression_matched_genes.bed ${dir}pombe_coding_minus_dcr1_matched_terminated.bed -b 500 -a 500 --skipZeros -bs 50 -p 15 -o ./computeMatrix/TES_Rad52_S_exp_dcr1_WT_compare_term_gene_matrix.mat.gz  --maxThreshold 500
plotHeatmap  -m ./computeMatrix/TES_Rad52_S_exp_dcr1_WT_compare_term_gene_matrix.mat.gz -out ./computeMatrix/TES_Rad52_S_exp_dcr1_WT_compare_term_heatmap.pdf --colorMap RdBu_r --dpi 300 --boxAroundHeatmaps no --heatmapHeight 9 --heatmapWidth 5 --missingDataColor 1 --regionsLabel term match other --samplesLabel Rad52_S_exp_dcr1_wt

####boxplot
enrich_matrix<-read.table("./computeMatrix/TES_Rad52_S_exp_dcr1_WT_compare_term_gene_matrix.mat.gz",skip=1,header=F,sep="\t")
# "group_labels":["dcr1_terminated_genes.bed","dcr1_expression_matched_genes.bed","pombe_coding_minus_dcr1_matched_terminated.bed"],"group_boundaries":[0,207,585,4096],"sample_labels":["Rad52_S_exp_dcr1_WT"],"sample_boundaries":[0,20]}
cotrol_enrich_matrix<-enrich_matrix[,1:26]
cotrol_enrich_matrix[is.na(cotrol_enrich_matrix)]<-0
cotrol_enrich_mean<-data.frame(gene=cotrol_enrich_matrix[,4],CPM_sum=rowSums( cotrol_enrich_matrix[,7:26]))
cotrol_enrich_mean$Features<-"Dcr1-terminated"
cotrol_enrich_mean$Features[205:585]<-"Expression-matched"
cotrol_enrich_mean$Features[586:4096]<-"others"

my_comparisons <- list( c("Dcr1-terminated", "Expression-matched"), c("Dcr1-terminated", "others"), c("Expression-matched", "others") )
p<-ggboxplot(cotrol_enrich_mean, x = "Features", y = "CPM_sum",
          color = "Features", palette = "jco",outlier.shape = NA,bxp.errorbar=T)+ 
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
  
  pdf("./computeMatrix/boxplot/TES_Rad52_S_exp_dcr1_WT_compare_boxplot_CPM.pdf",width=4,height=4)
# Visualize: Specify the comparisons you want
  ggpar(p,
        legend = "right",
        xlab = "Gene type",
        ylab = "Rad52 RPM S_log2(dcr1Δ/WT)-exp_log2(dcr1Δ/WT)")+rotate_x_text(11)
  # Add global p-value
dev.off()

write.table(cotrol_enrich_mean,"./computeMatrix/boxplot/sup_table_TES_Rad52_S_exp_dcr1_WT_compare_box.tsv",quote=F,row.names=F,sep="\t")