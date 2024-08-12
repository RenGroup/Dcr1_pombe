####Figure1A dcr1 RIP-seq termination index
#shell

computeMatrix scale-regions -S $RIP_dir -R ${dir}dcr1_terminated_genes.bed ${dir}dcr1_expression_matched_genes.bed ${dir}pombe_coding_minus_dcr1_matched_terminated.bed -b 500 -a 500 --regionBodyLength 1000 --skipZeros -bs 100 -p 15 -o ./term_index/RIP_upd500_body1000_bs100_matrix.mat.gz --maxThreshold 500

#R
library(ggplot2)
library(ggpubr)
enrich_matrix<-read.table("./RIP_upd500_body1000_bs100_matrix.mat.gz",skip=1,header=F,sep="\t")
# "group_labels":["dcr1_terminated_genes.bed","dcr1_expression_matched_genes.bed","pombe_coding_minus_dcr1_matched_terminated.bed"],"group_boundaries":[0,235,696,5159],"sample_labels":["Dcr1_RIP_IP_input"],"sample_boundaries":[0,20]}

cotrol_enrich_matrix<-enrich_matrix[,c(1:26)]
cotrol_enrich_matrix[is.na(cotrol_enrich_matrix)]<-0

#v12开始是TSS 22开始是TES
control_matrix<-cotrol_enrich_matrix
control_matrix[which(control_matrix<0,arr.ind = T)]=0

# Max(TES~TES+500)/mean(TSS+300~TES)
control_matrix$body<-rowMeans(control_matrix[,15:21])
control_matrix$end<-apply(control_matrix[,22:26],1,max)
#body<0 -> 0
control_matrix$body[which(control_matrix$body<0)]=0

control_matrix$type<-"Dcr1-terminated genes"
control_matrix$type[236:696]<-"Expression-matched genes"
control_matrix$type[697:5159]<-"All other genes"

control_index_matrix<-subset(control_matrix,select = c("V1","V2","V3","V4","body","end","type"))
control_index_matrix$term_index<-control_index_matrix$end-control_index_matrix$body
control_index_matrix$term_index[which(control_index_matrix$term_index<0)]<-0

my_comparisons <- list( c("Dcr1-terminated genes", "Expression-matched genes"), c("Dcr1-terminated genes", "All other genes"), c("Expression-matched genes", "All other genes") )
p<-ggboxplot(subset(control_index_matrix), x = "type", y = "term_index",
             color = "type", palette = "jco",outlier.shape = NA,bxp.errorbar=T)+ 
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")

ggpar(p,
      legend = "right",
      xlab = "Gene type",
      ylab = "Dcr1 RIP-seq termination index")+rotate_x_text(11)
  # Add global p-value
write.table(subset(control_index_matrix,select = c("V1","V2","V3","V4","type","term_index")),"supplemental_table_RIP_termindex_240516.tsv",quote = F,row.names = F,sep="\t")


####Figure1B
####RNAPII compare boxplot

computeMatrix reference-point --referencePoint TES -S ${RNAPII_dir} -R ${dir}dcr1_terminated_genes.bed ${dir}dcr1_expression_matched_genes.bed ${dir}pombe_coding_minus_dcr1_matched_terminated.bed -b 500 -a 500 --skipZeros -bs 50 -p 15 -o ./computeMatrix/RNAPII_compare2023_TES_select_matrix.mat.gz --maxThreshold 500
plotProfile -m ./computeMatrix/RNAPII_compare2023_TES_select_matrix.mat.gz -out ./computeMatrix/RNAPII_compare2023_TES_plot.pdf --plotFileFormat pdf --outFileNameData ./computeMatrix/RNAPII_compare2023_TES_myProfile.tab --regionsLabel term match other --samplesLabel RNAPII_compare

#boxplot
enrich_matrix<-read.table("./computeMatrix/RNAPII_compare2023_TES_select_matrix.mat.gz",skip=1,header=F,sep="\t")
#"group_labels":["dcr1_terminated_genes.bed","dcr1_expression_matched_genes.bed","pombe_coding_minus_dcr1_matched_terminated.bed"],"group_boundaries":[0,235,698,5205],"sample_labels":["dcr1_WT_compare"],"sample_boundaries":[0,20]}

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
  
  pdf("./RNAPII_TESsum_boxplot2023_240514_CPM.pdf",width=4,height=4)
# Visualize: Specify the comparisons you want
  ggpar(p,
        legend = "right",
        xlab = "Gene type",
        ylab = "RNAPII RPM log2(dcr1Δ/WT)")+rotate_x_text(11)
  # Add global p-value
dev.off()
write.table(cotrol_enrich_mean,"./supplemental_table_fig1B_RNAPII_dcr1_wt_compareBoxplot.tsv",quote = F,row.names = F,sep="\t")


####Figure1C
cd /work/home/path/Yeast/dcr1_Figure/ssDrip/boxplot
term_dir=/work/home/path/ref/yeast/Dcr1_matched_terminated_other/
computeMatrix reference-point --referencePoint TES -S $ssDrip_WT_dir $DRIPc_WT_RNSIII_dir -R ${term_dir}dcr1_terminated_genes.bed ${term_dir}dcr1_expression_matched_genes.bed ${term_dir}pombe_coding_minus_dcr1_matched_terminated.bed -b 500 -a 500 --skipZeros -bs 50 -p 15 -o ./WT_ssDRIP_DRIPcIII_TES_matrix.mat.gz --missingDataAsZero --maxThreshold 1000 
# "group_labels":["dcr1_terminated_genes.bed","dcr1_expression_matched_genes.bed","pombe_coding_minus_dcr1_matched_terminated.bed"],"group_boundaries":[0,235,698,5195],"sample_labels":["WT_ssDRIP","WT_Dripc"],"sample_boundaries":[0,20,40]}

#R
library(ggplot2)
library(ggpubr)
library(dplyr)
library(rstatix)
enrich_matrix<-read.table("WT_ssDRIP_DRIPcIII_TES_matrix.mat.gz",skip=1,header=F,sep="\t")
# "group_labels":["dcr1_terminated_genes.bed","dcr1_expression_matched_genes.bed","pombe_coding_minus_dcr1_matched_terminated.bed"],"group_boundaries":[0,235,698,5195],"sample_labels":["DG20_1_2_DG21_1_srs_chr","S06_S6_sort_rm_sort"],"sample_boundaries":[0,20,40]}

cotrol_enrich_mean<-data.frame(gene=enrich_matrix[,4],CPM_mean=rowSums(enrich_matrix[,7:26]))
cotrol_enrich_mean$Features<-"Dcr1-terminated"
cotrol_enrich_mean$Features[236:698]<-"Expression-matched"
cotrol_enrich_mean$Features[699:5195]<-"others"
treat_enrich_mean<-data.frame(gene=enrich_matrix[,4],CPM_mean=rowSums(enrich_matrix[,27:46]))
treat_enrich_mean$Features<-"Dcr1-terminated"
treat_enrich_mean$Features[236:695]<-"Expression-matched"
treat_enrich_mean$Features[696:5195]<-"others"
##same figure
cotrol_enrich_mean$class<-"ssDRIP"
treat_enrich_mean$class<-"DRIPc"
treat_enrich_mean.tmp<-subset(treat_enrich_mean)

cotrol_enrich_mean.tmp<-cotrol_enrich_mean

cotrol_enrich_mean.tmp$CPM_mean<-(cotrol_enrich_mean$CPM_mean)/1.3
all_enrich_mean<-rbind(treat_enrich_mean.tmp,cotrol_enrich_mean.tmp)
my_comparisons <- list( c("Dcr1-terminated", "Expression-matched"), c("Dcr1-terminated", "others"), c("Expression-matched", "others") )
all_enrich_mean$CPM_mean<-as.numeric(all_enrich_mean$CPM_mean)
stat.test <- all_enrich_mean %>%
  group_by(class) %>%
  wilcox_test(CPM_mean ~ Features, comparisons = my_comparisons) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p")
stat.test <- stat.test %>%
  add_xy_position(x = "class", dodge = 1) 
print(stat.test, width=Inf)
p<-ggboxplot(all_enrich_mean, x = "class", y = "CPM_mean",
             color = "Features", palette = "jco",outlier.shape = NA,bxp.errorbar=T)+ 
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0, step.increase = 0)  +coord_cartesian(ylim =  c(-0,250))+
  scale_y_continuous(name = "DRIPc-seq RPM",#y1
                     sec.axis = sec_axis( trans=~((.)*1.3), name="ssDRIP-seq IP RPM")) #y2
pdf("./WT_ssDRIP_DRIPcIII_TESsum_240516_boxplot.pdf",width=6,height=4)
  ggpar(p,
        legend = "right",
        xlab = "Gene type",
        ylab = "DRIPc-seq RPM")+rotate_x_text(11)
dev.off()

all_enrich_mean$CPM_mean[which(all_enrich_mean$class=="DRIPc")]<-all_enrich_mean$CPM_mean[which(all_enrich_mean$class=="DRIPc")]*1.3
write.table(all_enrich_mean,"./supplemental_table_fig1C_ssDRIP_DRIPc_Boxplot.tsv",quote = F,row.names = F,sep="\t")

#####Fig1D metaplot
mkdir computeMatrix
computeMatrix reference-point --referencePoint TES -S ${RIP_merge_dir} ${RNAPII_dir} $ssDrip_compare_dir -R ${dir}dcr1_terminated_genes.bed ${dir}dcr1_expression_matched_genes.bed ${dir}pombe_coding_minus_dcr1_matched_terminated.bed -b 500 -a 500 --skipZeros -bs 50 -p 15 -o ./computeMatrix/WT_meta_TES_240509_select_matrix.mat.gz --maxThreshold 500
plotProfile -m ./computeMatrix/WT_meta_TES_240509_select_matrix.mat.gz -out ./computeMatrix/WT_meta_TES_240509_select_plot.pdf --plotFileFormat pdf --outFileNameData ./computeMatrix/WT_meta_TES_240509_myProfile.tab --perGroup --regionsLabel term match other --samplesLabel RIP_merge12 RNAPII_compare ssDRIP_compare

#smooth
library(reshape2)
tab.data<-read.table("WT_meta_TES_240509_myProfile.tab",header = F,sep="\t",fill=T)
tab.data<-t(tab.data)
dim(tab.data)
#62 11
tab.rowname<-tab.data[3:22,2]
tab.data.term<-as.data.frame(tab.data[3:22,c(3,6,9)])

rownames(tab.data.term)<-tab.rowname
colnames(tab.data.term)<-c("dcr1_RIP","RNAPII_log2(dcr1Δ/WT)","ssDRIP_log2(dcr1Δ/WT)")
tab.data.term$'RNAPII_log2(dcr1Δ/WT)'<-(as.numeric(tab.data.term$'RNAPII_log2(dcr1Δ/WT)')/2)
tab.data.term$loc<-tab.rowname
tab.data.term<-melt(tab.data.term,measure.vars=c("dcr1_RIP","RNAPII_log2(dcr1Δ/WT)","ssDRIP_log2(dcr1Δ/WT)"))
colnames(tab.data.term)<-c("loc","type","value")
pdf("/work/home/path/Yeast/dcr1_Figure/Metadata/computeMatrix/WT_meta_240509_TES_smooth.pdf")
ggplot(tab.data.term,aes(as.numeric(loc),as.numeric(value),colour=type))+ 
  geom_smooth(span = 0.5, se = F)+
  scale_x_continuous(breaks=c(1,10,20), labels =c("-500bp","TES","500bp"))+
  labs(x = "Position",y="RPM")+theme_bw()+
  theme(panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+  
  scale_y_continuous(name = "RPM",#y1
                     sec.axis = sec_axis( trans=~((.)*2), name="RNAPII_log2(dcr1Δ/WT)"))#y2
  dev.off() #()
