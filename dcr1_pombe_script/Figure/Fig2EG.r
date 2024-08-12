####Fig2E Dhp1 WT enrichment
computeMatrix scale-regions -S $WT_dir -R ${dir}dcr1_terminated_genes.bed ${dir}dcr1_expression_matched_genes.bed ${dir}pombe_coding_minus_dcr1_matched_terminated.bed -b 1000 -a 1000 --regionBodyLength 1000 --skipZeros -bs 50 -p 15 -o ./computeMatrix/dhp1_WTm_enrichment_matrix.mat.gz --maxThreshold 1000

plotProfile -m ./computeMatrix/dhp1_WTm_enrichment_matrix.mat.gz -out ./computeMatrix/dhp1_WTm_enrichment_plot.pdf --plotFileFormat pdf --outFileNameData ./computeMatrix/dhp1_WTm_enrichment_myProfile.tab

#smooth
library(ggplot2)
library(tidyr)
library(cowplot)
library(dplyr)
library(reshape2)
library(ggpubr)
tab.data<-read.table("dhp1_WTm_enrichment_myProfile.tab",header = F,sep="\t")
tab.data<-t(tab.data)
dim(tab.data)
# [1] 62  5
tab.rowname<-tab.data[3:62,2]
tab.data.term<-as.data.frame(tab.data[3:62,c(2,3)])
tab.data.match<-as.data.frame(tab.data[3:62,c(2,4)])
tab.data.other<-as.data.frame(tab.data[3:62,c(2,5)])
tab.data.term$type<-c("Dcr1-terminated genes")
tab.data.match$type<-c("Expression-matched genes")
tab.data.other$type<-c("All other genes")
rownames(tab.data.term)<-tab.rowname
rownames(tab.data.match)<-tab.rowname
rownames(tab.data.other)<-tab.rowname
colnames(tab.data.term)<-c("loc","Value","type")
colnames(tab.data.match)<-c("loc","Value","type")
colnames(tab.data.other)<-c("loc","Value","type")
all_data<-rbind(tab.data.term,tab.data.match,tab.data.other)

pdf("/path/computeMatrix/dhp1_WTm_enrichment_smooth.pdf")
ggplot(all_data,aes(as.numeric(loc),as.numeric(Value),colour=type))+ 
  geom_smooth(span = 0.3, se = F)+
  scale_x_continuous(breaks=c(1,20,40,60), labels =c("-1000bp","TSS","TES","1000bp"))+
  labs(x = "Position", y = "dhp1 WT ChIP-seq log2(IP/Input)")+theme_bw()+
  theme(panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank(),  
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black"))
  dev.off() #()

####Fig2G dhp1 WT/dcr1Δ compare boxplot

computeMatrix reference-point --referencePoint TES -S ${compare_dir} -R ${dir}dcr1_terminated_genes.bed ${dir}dcr1_expression_matched_genes.bed ${dir}pombe_coding_minus_dcr1_matched_terminated.bed -b 500 -a 500 --skipZeros -bs 50 -p 15 -o ./computeMatrix/dhp1_TES_WTm_dcr1m_compare_matrix.mat.gz --maxThreshold 1000

#boxplot
enrich_matrix<-read.table("./computeMatrix/dhp1_TES_WTm_dcr1m_compare_matrix.mat.gz",skip=1,header=F,sep="\t")
#"group_labels":["dcr1_terminated_genes.bed","dcr1_expression_matched_genes.bed","pombe_coding_minus_dcr1_matched_terminated.bed"],"group_boundaries":[0,235,698,5205],"sample_labels":["dhp1_WTm_dcr1\u0394m_compare"],"sample_boundaries":[0,20]}


cotrol_enrich_matrix<-enrich_matrix[,1:26]
cotrol_enrich_matrix[is.na(cotrol_enrich_matrix)]<-0
cotrol_enrich_mean<-data.frame(gene=cotrol_enrich_matrix[,4],CPM_mean=rowSums(cotrol_enrich_matrix[,7:26]))
cotrol_enrich_mean$Features<-"Dcr1-terminated"
cotrol_enrich_mean$Features[236:698]<-"Expression-matched"
cotrol_enrich_mean$Features[699:5205]<-"others"

my_comparisons <- list( c("Dcr1-terminated", "Expression-matched"), c("Dcr1-terminated", "others"), c("Expression-matched", "others") )
p<-ggboxplot(cotrol_enrich_mean, x = "Features", y = "CPM_mean",
          color = "Features", palette = "jco",outlier.shape = NA,bxp.errorbar=T)+ 
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
  
  pdf("./computeMatrix/boxplot/dhp1_wt_dcr1d_compare_TESsum_boxplot_240606_CPM.pdf",width=4,height=4)
  ggpar(p,
        legend = "right",
        xlab = "Gene type",
        ylab = "dhp1 RPM log2(WT/dcr1Δ)")+rotate_x_text(11)
  # Add global p-value
dev.off()
write.table(cotrol_enrich_mean,"./computeMatrix/boxplot/supplemental_table_dhp1_wt_dcr1_compareBoxplot.tsv",quote = F,row.names = F,sep="\t")