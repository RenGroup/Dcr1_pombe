###################################figure1A
#### shell 
##term sense
computeMatrix scale-regions -S $RIP_fwd -R $term_fwd -b 500 -a 500 --regionBodyLength 1000 -bs 100 -p 15 -o ./computeMatrix/Dcr1_RIP_term_fwd+_matrix.mat.gz --missingDataAsZero --maxThreshold 500
computeMatrix scale-regions -S $RIP_rev -R $term_rev -b 500 -a 500 --regionBodyLength 1000 -bs 100 -p 15 -o ./computeMatrix/Dcr1_RIP_term_rev-_matrix.mat.gz --missingDataAsZero --maxThreshold 500

# Prior to merging, it is necessary to manually modify the group label names in the matrix; otherwise, all entries will be consolidated into a single group during the final merge.
computeMatrixOperations rbind -m ./computeMatrix/Dcr1_RIP_term_fwd+_matrix.mat.gz ./computeMatrix/Dcr1_RIP_term_rev-_matrix.mat.gz -o ./computeMatrix/Dcr1_RIP_term_sense.mat.gz

##match sense 
computeMatrix scale-regions -S $RIP_fwd -R $match_fwd -b 500 -a 500 --regionBodyLength 1000 -bs 100 -p 15 -o ./computeMatrix/Dcr1_RIP_match_fwd+_matrix.mat.gz --missingDataAsZero --maxThreshold 500
computeMatrix scale-regions -S $RIP_rev -R $match_rev -b 500 -a 500 --regionBodyLength 1000 -bs 100 -p 15 -o ./computeMatrix/Dcr1_RIP_match_rev-_matrix.mat.gz --missingDataAsZero --maxThreshold 500
computeMatrixOperations rbind -m ./computeMatrix/Dcr1_RIP_match_fwd+_matrix.mat.gz ./computeMatrix/Dcr1_RIP_match_rev-_matrix.mat.gz -o ./computeMatrix/Dcr1_RIP_match_sense.mat.gz

##other sense
computeMatrix scale-regions -S $RIP_fwd -R $other_fwd -b 500 -a 500 --regionBodyLength 1000 -bs 100 -p 15 -o ./computeMatrix/Dcr1_RIP_other_fwd+_matrix.mat.gz --missingDataAsZero --maxThreshold 500
computeMatrix scale-regions -S $RIP_rev -R $other_rev -b 500 -a 500 --regionBodyLength 1000 -bs 100 -p 15 -o ./computeMatrix/Dcr1_RIP_other_rev-_matrix.mat.gz --missingDataAsZero --maxThreshold 500
computeMatrixOperations rbind -m ./computeMatrix/Dcr1_RIP_other_fwd+_matrix.mat.gz ./computeMatrix/Dcr1_RIP_other_rev-_matrix.mat.gz -o ./computeMatrix/Dcr1_RIP_other_sense.mat.gz

computeMatrixOperations rbind -m ./computeMatrix/Dcr1_RIP_term_sense.mat.gz ./computeMatrix/Dcr1_RIP_match_sense.mat.gz ./computeMatrix/Dcr1_RIP_other_sense.mat.gz -o ./computeMatrix/Dcr1_RIP_sense.mat.gz

#### R 
library(ggplot2)
library(ggpubr)
enrich_matrix<-read.table("./Dcr1_RIP_sense.mat.gz",skip=1,header=F,sep="\t")
# "group_labels":["term","match","other"],"group_boundaries":[0,207,585,4096], "sample_boundaries":[0,20]}
cotrol_enrich_matrix<-enrich_matrix[,c(1:26)]
cotrol_enrich_matrix[is.na(cotrol_enrich_matrix)]<-0

#v12 TSS 22 TES
control_matrix<-cotrol_enrich_matrix
control_matrix[which(control_matrix<0,arr.ind = T)]=0

# Max(TES~TES+100)/mean(TSS+300~TES)
control_matrix$body<-rowMeans(control_matrix[,15:21])
control_matrix$body[which(control_matrix$body<0)]=0
control_matrix$end<-apply(control_matrix[,22:23],1,max)

control_matrix$type<-"Dcr1-terminated genes"
control_matrix$type[208:585]<-"Expression-matched genes"
control_matrix$type[586:4096]<-"All other genes"

control_index_matrix<-subset(control_matrix,select = c("V1","V2","V3","V4","body","end","type"))
control_index_matrix$term_index<-control_index_matrix$end-control_index_matrix$body
control_index_matrix$term_index[which(control_index_matrix$term_index<0)]<-0

my_comparisons <- list( c("Dcr1-terminated genes", "Expression-matched genes"), c("Dcr1-terminated genes", "All other genes"), c("Expression-matched genes", "All other genes") )
p<-ggboxplot(subset(control_index_matrix), x = "type", y = "term_index",
             color = "type", palette = "jco",outlier.shape = NA,bxp.errorbar=T)+ 
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")

# Visualize: Specify the comparisons you want
ggpar(p,
      legend = "right",
      xlab = "Gene type",
      ylab = "Dcr1 RIP-seq termination index")+rotate_x_text(11)
# Add global p-value
write.table(subset(control_index_matrix,select = c("V1","V2","V3","V4","type","term_index")),"supplemental_table_RIP_termindex.tsv",quote = F,row.names = F,sep="\t")


###################################figure1B
###############Fig1B RNAPII compare boxplot
computeMatrix reference-point --referencePoint TES -S ${RNAPII_dir} -R ${term_dir}dcr1_terminated_genes.bed ${term_dir}dcr1_expression_matched_genes_chr.bed ${term_dir}pombe_coding_minus_dcr1_matched_terminated.bed.bed -b 200 -a 200 --skipZeros -bs 50 -p 15 -o ./computeMatrix/RNAPII_compare2023_TES200_select_matrix.mat.gz --maxThreshold 500
 
#boxplot
enrich_matrix<-read.table("./computeMatrix/RNAPII_compare2023_TES200_select_matrix.mat.gz",skip=1,header=F,sep="\t")
#"group_labels":["dcr1_terminated_genes.bed","dcr1_expression_matched_genes_chr.bed","pombe_coding_minus_dcr1_matched_terminated.bed.bed"],"group_boundaries":[0,207,585,4096],"sample_labels":["dcr1_2_WT1_compare"],"sample_boundaries":[0,8]}

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
  
  pdf("./RNAPII_TES200_boxplot2023_CPM.pdf",width=4,height=4)
# Visualize: Specify the comparisons you want
  ggpar(p,
        legend = "right",
        xlab = "Gene type",
        ylab = "RNAPII RPM log2(dcr1Δ/WT)")+rotate_x_text(11)
  # Add global p-value
dev.off()
 
write.table(cotrol_enrich_mean,"./boxplot/sup_table_fig1B_RNAPII_TES200_dcr1_wt_compareBoxplot.tsv",quote = F,row.names = F,sep="\t")


###################################figure1C
computeMatrix reference-point --referencePoint TES -S $ssDrip_WT_dir $wt_ssDRIP_RNH -R ${term_dir}dcr1_terminated_genes.bed ${term_dir}dcr1_expression_matched_genes_chr.bed ${term_dir}pombe_coding_minus_dcr1_matched_terminated.bed.bed -b 200 -a 200 --skipZeros -bs 50 -p 15 -o ./WT_ssDRIP_RNH_TES200_matrix.mat.gz --missingDataAsZero --maxThreshold 500 

### R
 ##R 50bp sum
library(ggplot2)
library(ggpubr)
library(dplyr)
library(rstatix)
enrich_matrix<-read.table("WT_ssDRIP_RNH_TES200_matrix.mat.gz",skip=1,header=F,sep="\t")
# ,"group_boundaries":[0,207,585,4096],"sample_labels":["DG20_1_2_DG21_1_srs_chr","DG20H_1_2_DG21H_1_2_srs_chr"],"sample_boundaries":[0,8,16]

cotrol_enrich_mean<-data.frame(gene=enrich_matrix[,4],CPM_sum=rowSums(enrich_matrix[,7:14]))
cotrol_enrich_mean$Features<-"Dcr1-terminated"
cotrol_enrich_mean$Features[208:585]<-"Expression-matched"
cotrol_enrich_mean$Features[586:4096]<-"others"
treat_enrich_mean<-data.frame(gene=enrich_matrix[,4],CPM_sum=rowSums(enrich_matrix[,15:22]))
treat_enrich_mean$Features<-"Dcr1-terminated"
treat_enrich_mean$Features[208:585]<-"Expression-matched"
treat_enrich_mean$Features[586:4096]<-"others"
##
cotrol_enrich_mean$class<-"ssDRIP"
treat_enrich_mean$class<-"ssDRIP_RH"
treat_enrich_mean.tmp<-subset(treat_enrich_mean)

cotrol_enrich_mean.tmp<-cotrol_enrich_mean
all_enrich_mean<-rbind(treat_enrich_mean.tmp,cotrol_enrich_mean.tmp)
my_comparisons <- list( c("Dcr1-terminated", "Expression-matched"), c("Dcr1-terminated", "others"), c("Expression-matched", "others") )
all_enrich_mean$CPM_sum<-as.numeric(all_enrich_mean$CPM_sum)
stat.test <- all_enrich_mean %>%
  group_by(class) %>%
  wilcox_test(CPM_sum ~ Features, comparisons = my_comparisons) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

stat.test <- stat.test %>%
  mutate(class = factor(class, levels = unique(all_enrich_mean$class))) %>%
  arrange(class)

stat.test <- stat.test %>%
    rstatix::add_xy_position(x = "class")
stat.test$xmin<-sort(stat.test$xmin)
stat.test$xmax<-sort(stat.test$xmax)

print(stat.test, width=Inf)
p<-ggboxplot(all_enrich_mean, x = "class", y = "CPM_sum",
             color = "Features", palette = "jco",outlier.shape = NA,bxp.errorbar=T)+ 
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0, step.increase = 0)+ coord_cartesian(ylim =  c(0,200))

pdf("./WT_ssDRIP_RNH_IP_TES200sum_boxplot.pdf")
  ggpar(p,
        legend = "right",
        xlab = "Gene type",
        ylab = "RPM")+rotate_x_text(11)
dev.off()
write.table(all_enrich_mean,"./supplemental_table_fig1C_ssDRIP_wt_RH_TES200_Boxplot.tsv",quote = F,row.names = F,sep="\t")

###################################figure1D
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


###################################figure1H
computeMatrix reference-point --referencePoint TES -S other_1.bw SPB455_repm_IP_input.bw $other_2 $RNAPII_Dcr1_dir -b 500 -a 500 -R ${dir}dcr1_terminated_genes.bed ${dir}dcr1_expression_matched_genes.bed ${dir}pombe_coding_minus_dcr1_matched_terminated.bed --skipZeros -bs 10 -p 10 -o ./computeMatrix/TES_term_gene_matrix.mat.gz

#dot plot  
library(ggplot2)
library(ggpubr)

# setwd("D:/master/task/yeast_dicer/Figure/RNAPII/C103/computeMatrix")
enrich_matrix<-read.table("./TES_term_gene_matrix.mat.gz",skip=1,header=F,sep="\t")
# ,"group_labels":["dcr1_terminated_genes.bed","dcr1_expression_matched_genes.bed","pombe_coding_minus_dcr1_matched_terminated.bed"],"sample_labels":["other_1","SPB455_repm_IP_input","other_2","dcr1_pS2_rep2IP_input"],"sample_boundaries":[0,100,200,300,400]}

term_gene=read.table("/work/home/path/ref/yeast/Dcr1_matched_terminated_other/Dcr1_terminated_genes.bed",header=F,sep="\t")

term_gene_loc=paste(term_gene[,1],term_gene[,2],term_gene[,3],sep=":")

enrich_matrix[is.na(enrich_matrix)]<-0
control_mean<-rowMeans(enrich_matrix[,107:206])
treat_mean<-rowMeans(enrich_matrix[,307:406])
plot_matrix<-data.frame(loc=paste(enrich_matrix[,1],enrich_matrix[,2],enrich_matrix[,3],sep=":"),treat_mean_CPM=as.numeric(treat_mean),control_mean_CPM=as.numeric(control_mean))
plot_matrix$type="others"
plot_matrix$type[match(term_gene_loc,plot_matrix$loc)]<-"Dcr1-terminated genes"
plot_sort<-rbind(plot_matrix[-match(term_gene_loc,plot_matrix$loc),],plot_matrix[match(term_gene_loc,plot_matrix$loc),])

pdf("SPB455_dcr1d_TES_cor.pdf")
ggscatter(plot_sort, x = "control_mean_CPM", y = "treat_mean_CPM",color = "type",
          add = "reg.line",
          alpha = 0.5,
          palette = c("red", "grey"),
          add.params = list(color = "black",fill = "lightgray",linetype="dashed")
)+
  stat_cor(method = "pearson", 
           label.x = 0, label.y = 3.8)+
  labs(x = "SPB455 PNAPII RPM log2(IP/Input)", y = "dcr1Δ RNAPII RPM log2(IP/Input)")
dev.off() #()
