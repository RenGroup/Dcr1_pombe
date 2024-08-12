####FigS1B
#ori RIP $ssDRIP compare dotplot

computeMatrix scale-regions -S $RIP_merge_dir $ssDRIP_compare -R ${ori_dir}replication_origin.bed -b 0 -a 0 --regionBodyLength 1000 --skipZeros -bs 10 -p 10 -o ./computeMatrix/ori_RIP_ssDRIP_compare_noupd_matrix.mat.gz --maxThreshold 500 
# "group_labels":["genes"],"group_boundaries":[0,554],"sample_labels":["Dcr1_RIP_merge12_IP_input","dcr1_WT_fc_ssDrip"],"sample_boundaries":[0,100,200]}
enrich_matrix<-read.table("./computeMatrix/ori_RIP_ssDRIP_compare_noupd_matrix.mat.gz",skip=1,header=F,sep="\t")
enrich_matrix[is.na(enrich_matrix)]<-0
control_mean<-rowMeans(enrich_matrix[,7:106])
treat_mean<-rowMeans(enrich_matrix[,107:206])
plot_matrix<-data.frame(loc=paste(enrich_matrix[,1],enrich_matrix[,2],enrich_matrix[,3],sep=":"),treat_mean_CPM=as.numeric(treat_mean),control_mean_CPM=as.numeric(control_mean))

pdf("./ori_RIP_ssDRIP_compare_noupd_dotplot.pdf",width=4,height=4)
ggplot(plot_matrix,aes(x=control_mean_CPM,y=treat_mean_CPM))+
geom_point(alpha = 0.5)+
geom_hline(yintercept = 0,linetype = "dashed")+
geom_vline(xintercept = 0,linetype = "dashed")+
  labs(x = "Dcr1 RIP-seq log2(IP/Input)", y = "ssDRIP-seq log2(dcr1Δ/WT)")+theme_bw()+
  theme(panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"))  
dev.off()


#####FigS1D #sRNA boxplot

computeMatrix reference-point --referencePoint TES -S ${sRNA_dir} -R ${dir}dcr1_terminated_genes.bed ${dir}dcr1_expression_matched_genes.bed ${dir}pombe_coding_minus_dcr1_matched_terminated.bed -b 500 -a 500 --skipZeros -bs 50 -p 15 -o ./computeMatrix/TES_sRNA_WT_dcr1_minus240510_matrix.mat.gz --minThreshold -100 --maxThreshold 500

plotProfile -m ./computeMatrix/TES_sRNA_WT_dcr1_minus240510_matrix.mat.gz -out ./computeMatrix/TES_sRNA_WT_dcr1_minus240510_plot.pdf --plotFileFormat pdf --outFileNameData ./computeMatrix/TES_sRNA_WT_dcr1_minus240510_matrix_myProfile.tab


#boxplot
cd computeMatrix/
enrich_matrix<-read.table("./TES_sRNA_WT_dcr1_minus240510_matrix.mat.gz",skip=1,header=F,sep="\t")
# "group_labels":["dcr1_terminated_genes.bed","dcr1_expression_matched_genes.bed","pombe_coding_minus_dcr1_matched_terminated.bed"],"group_boundaries":[0,234,696,5115],"sample_labels":["dg21_dg690_minus"],"sample_boundaries":[0,20]}

cotrol_enrich_matrix<-enrich_matrix[,1:26]
cotrol_enrich_matrix[is.na(cotrol_enrich_matrix)]<-0
cotrol_enrich_mean<-data.frame(gene=cotrol_enrich_matrix[,4],CPM_mean=rowSums(cotrol_enrich_matrix[,7:26]))
cotrol_enrich_mean$Features<-"Dcr1-terminated"
cotrol_enrich_mean$Features[236:696]<-"Expression-matched"
cotrol_enrich_mean$Features[695:5115]<-"others"

my_comparisons <- list( c("Dcr1-terminated", "Expression-matched"), c("Dcr1-terminated", "others"), c("Expression-matched", "others") )
p<-ggboxplot(cotrol_enrich_mean, x = "Features", y = "CPM_mean",
          color = "Features", outlier.shape = NA,palette = "jco",bxp.errorbar=T)+ 
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")+coord_cartesian(ylim =  c(-10, 40))

  pdf("./computeMatrix/boxplot/TES_sRNA_WT_dcr1_minus_240513binsum_boxplot_CPM.pdf",width=4,height=4)
# Visualize: Specify the comparisons you want
  ggpar(p,
        legend = "right",
        xlab = "Gene type",
        ylab = "sRNA WT RPM - dcr1Δ RPM)")+rotate_x_text(11)
  # Add global p-value
dev.off()
write.table(cbind(cotrol_enrich_matrix[,1:3],cotrol_enrich_mean),"./supplemental_table_sRNA_boxplot.tsv",quote = F,row.names = F,sep="\t")


#####FigS1G####ssDRIP wt dcr1d 2d dotplot
computeMatrix reference-point --referencePoint TES -S $ssDrip_WT_dir $ssDrip_Dcr1_dir -R $conding_dir -b 500 -a 500 --skipZeros -bs 50 -p 15 -o ./ssDrip_coding_TES500_matrix.mat.gz --maxThreshold 1000 --missingDataAsZero
#R
library(ggplot2)
library(ggpubr)
enrich_matrix<-read.table("./ssDrip_coding_TES500_matrix.mat.gz",skip=1,header=F,sep="\t")
# "group_labels":["genes"],"group_boundaries":[0,5205],"sample_labels":["DG20_21IP_input","DG690_691IP_input"],"sample_boundaries":[0,20,40]}

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

