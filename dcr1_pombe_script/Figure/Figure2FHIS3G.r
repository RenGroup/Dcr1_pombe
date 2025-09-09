###################################figure2F
#### shell 
 omputeMatrix scale-regions -S $WT_dir -R ${dir}dcr1_terminated_genes.bed ${dir}dcr1_expression_matched_genes.bed ${dir}pombe_coding_minus_dcr1_matched_terminated.bed -b 1000 -a 1000 --regionBodyLength 1000 --skipZeros -bs 50 -p 15 -o ./computeMatrix/dhp1_WTm_enrichment_matrix.mat.gz --maxThreshold 1000

plotProfile -m ./computeMatrix/dhp1_WTm_enrichment_matrix.mat.gz -out ./computeMatrix/dhp1_WTm_enrichment_plot.pdf --plotFileFormat pdf --outFileNameData ./computeMatrix/dhp1_WTm_enrichment_myProfile.tab

###R
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

###################################figure2H 2I S3G
#### shell 
computeMatrix reference-point --referencePoint TES -S ${dhp1_dir}dcr1Δ_dhp1_mean_IP_input.bw -R ${term_dir}dcr1_terminated_genes.bed ${term_dir}dcr1_expression_matched_genes_chr.bed ${term_dir}pombe_coding_minus_dcr1_matched_terminated.bed.bed -b 200 -a 200 --skipZeros -bs 50 -p 10 -o ./computeMatrix/dcr1m_Dhp1_enrich_TES200_matrix.mat.gz --maxThreshold 500 --missingDataAsZero

computeMatrix reference-point --referencePoint TES -S ${dhp1_dir}WT_dhp1_mean_IP_input.bw -R ${term_dir}dcr1_terminated_genes.bed ${term_dir}dcr1_expression_matched_genes_chr.bed ${term_dir}pombe_coding_minus_dcr1_matched_terminated.bed.bed -b 200 -a 200 --skipZeros -bs 50 -p 10 -o ./computeMatrix/WTm_Dhp1_enrich_TES200_matrix.mat.gz --maxThreshold 500 --missingDataAsZero

computeMatrix reference-point --referencePoint TES -S ${RNAPII_dir}dcr1_pS2_IP_input.bw -R ${term_dir}dcr1_terminated_genes.bed ${term_dir}dcr1_expression_matched_genes_chr.bed ${term_dir}pombe_coding_minus_dcr1_matched_terminated.bed.bed -b 200 -a 200 --skipZeros -bs 50 -p 10 -o ./computeMatrix/dcr1_RNAPII_enrich_TES200_matrix.mat.gz --maxThreshold 500 --missingDataAsZero

computeMatrix reference-point --referencePoint TES -S ${RNAPII_dir}WT_pS2_IP_input.bw -R ${term_dir}dcr1_terminated_genes.bed ${term_dir}dcr1_expression_matched_genes_chr.bed ${term_dir}pombe_coding_minus_dcr1_matched_terminated.bed.bed -b 200 -a 200 --skipZeros -bs 50 -p 10 -o ./computeMatrix/WT_RNAPII_enrich_TES200_matrix.mat.gz --maxThreshold 500 --missingDataAsZero

#R##

WT_Dhp1_matrix<-read.table("./computeMatrix/WTm_Dhp1_enrich_TES200_matrix.mat.gz",skip=1,header=F,sep="\t")
# "group_boundaries":[0,207,585,4096],"sample_labels":["WT_dhp1_mean_IP_input"],"sample_boundaries":[0,8]}
dcr1_Dhp1_matrix<-read.table("./computeMatrix/dcr1m_Dhp1_enrich_TES200_matrix.mat.gz",skip=1,header=F,sep="\t")
# "group_boundaries":[0,207,585,4096],"sample_labels":["dcr1_dhp1_mean_IP_input"],"sample_boundaries":[0,8]}
dcr1_RNAPII_matrix<-read.table("./computeMatrix/dcr1_RNAPII_enrich_TES200_matrix.mat.gz",skip=1,header=F,sep="\t")
# "group_boundaries":[0,207,585,4096],"sample_labels":["dcr1_pS2_IP_input"],"sample_boundaries":[0,8]}
WT_RNAPII_matrix<-read.table("./computeMatrix/WT_RNAPII_enrich_TES200_matrix.mat.gz",skip=1,header=F,sep="\t")
# "group_boundaries":[0,207,585,4096],"sample_labels":["WT_pS2_IP_input"],"sample_boundaries":[0,8]}
WT_Dhp1_enrich_matrix<-WT_Dhp1_matrix[,c(1:14)]
WT_Dhp1_enrich_matrix[is.na(WT_Dhp1_enrich_matrix)]<-0
WT_Dhp1_enrich_matrix<-data.frame(gene=WT_Dhp1_enrich_matrix[,4],CPM_sum=rowSums( WT_Dhp1_enrich_matrix[,7:14]))
WT_Dhp1_enrich_matrix$Features<-"Dcr1-terminated"
WT_Dhp1_enrich_matrix$Features[208:585]<-"Expression-matched"
WT_Dhp1_enrich_matrix$Features[586:4096]<-"others"

dcr1_Dhp1_enrich_matrix<-dcr1_Dhp1_matrix[,c(1:14)]
dcr1_Dhp1_enrich_matrix[is.na(dcr1_Dhp1_enrich_matrix)]<-0
dcr1_Dhp1_enrich_matrix<-data.frame(gene=dcr1_Dhp1_enrich_matrix[,4],CPM_sum=rowSums(dcr1_Dhp1_enrich_matrix[,7:14]))
dcr1_Dhp1_enrich_matrix$Features<-"Dcr1-terminated"
dcr1_Dhp1_enrich_matrix$Features[208:585]<-"Expression-matched"
dcr1_Dhp1_enrich_matrix$Features[586:4096]<-"others"

dcr1_RNAPII_enrich_matrix<-dcr1_RNAPII_matrix[,c(1:14)]
dcr1_RNAPII_enrich_matrix[is.na(dcr1_RNAPII_enrich_matrix)]<-0
dcr1_RNAPII_enrich_matrix<-data.frame(gene=dcr1_RNAPII_enrich_matrix[,4],CPM_sum=rowSums(dcr1_RNAPII_enrich_matrix[,7:14]))
dcr1_RNAPII_enrich_matrix$Features<-"Dcr1-terminated"
dcr1_RNAPII_enrich_matrix$Features[208:585]<-"Expression-matched"
dcr1_RNAPII_enrich_matrix$Features[586:4096]<-"others"

WT_RNAPII_enrich_matrix<-WT_RNAPII_matrix[,c(1:14)]
WT_RNAPII_enrich_matrix[is.na(WT_RNAPII_enrich_matrix)]<-0
WT_RNAPII_enrich_matrix<-data.frame(gene=WT_RNAPII_enrich_matrix[,4],CPM_sum=rowSums(WT_RNAPII_enrich_matrix[,7:14]))
WT_RNAPII_enrich_matrix$Features<-"Dcr1-terminated"
WT_RNAPII_enrich_matrix$Features[208:585]<-"Expression-matched"
WT_RNAPII_enrich_matrix$Features[586:4096]<-"others"

WT_Dhp1_RNAPII<-WT_Dhp1_enrich_matrix
WT_Dhp1_RNAPII$CPM_sum<-(WT_Dhp1_enrich_matrix$CPM_sum)-(WT_RNAPII_enrich_matrix$CPM_sum)

dcr1d_Dhp1_RNAPII<-dcr1_Dhp1_enrich_matrix
dcr1d_Dhp1_RNAPII$CPM_sum<-(dcr1_Dhp1_enrich_matrix$CPM_sum)-(dcr1_RNAPII_enrich_matrix$CPM_sum)

WT_dcr1d_Dhp1_RNAPII<-dcr1_Dhp1_enrich_matrix
WT_dcr1d_Dhp1_RNAPII$CPM_sum<-(WT_Dhp1_RNAPII$CPM_sum)-(dcr1d_Dhp1_RNAPII$CPM_sum)

my_comparisons <- list( c("Dcr1-terminated", "Expression-matched"), c("Dcr1-terminated", "others"), c("Expression-matched", "others") )
# Replace different samples.
p<-ggboxplot(dcr1d_Dhp1_RNAPII, x = "Features", y = "CPM_sum",
          color = "Features", palette = "jco",outlier.shape = NA,bxp.errorbar=T)+ 
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
 
ggpar(p,
        legend = "right",
        xlab = "Gene type",
        ylab = "log2(dcr1d_Dhp1/dcr1d_RNAPII)")+rotate_x_text(11)
 

write.table(WT_Dhp1_RNAPII,"./sup_table_WT_Dhp1_normRNAPII_TES200_box.tsv",quote = F,row.names = F,sep="\t")
write.table(dcr1d_Dhp1_RNAPII,"./sup_table_dcr1d_Dhp1_normRNAPII_TES200_box.tsv",quote = F,row.names = F,sep="\t")
write.table(WT_dcr1d_Dhp1_RNAPII,"./sup_table_WT_dcr1d_Dhp1_normRNAPII_TES200_box.tsv",quote = F,row.names = F,sep="\t")
