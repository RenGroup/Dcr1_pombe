###################################Figure S1 CDE
#### shell 
####NET-seq
computeMatrix scale-regions -S $sample -R ${dir}dcr1_terminated_genes.bed ${dir}dcr1_expression_matched_genes.bed ${dir}pombe_coding_minus_dcr1_matched_terminated.bed -b 1000 -a 1000 --regionBodyLength 1000 --skipZeros -bs 50 -p 15 -o ./computeMatrix/termgene_compare_matrix.mat.gz --maxThreshold 1000 --missingDataAsZero
 
# # termgene_compare_matrix.mat.gz
# {"upstream":[1000],"downstream":[1000],"body":[1000],"bin size":[50],"ref point":[null],"verbose":false,"bin avg type":"mean","missing data as zero":true,"min threshold":null,"max threshold":1000.0,"scale":1,"skip zeros":true,"nan after end":false,"proc number":15,"sort regions":"keep","sort using":"mean","unscaled 5 prime":[0],"unscaled 3 prime":[0],"group_labels":["dcr1_terminated_genes.bed","dcr1_expression_matched_genes.bed","pombe_coding_minus_dcr1_matched_terminated.bed"],"group_boundaries":[0,207,585,4096],"sample_labels":["dcr1_WT_0619_0521_compare_mean_skip0"],"sample_boundaries":[0,60]}

#R scale by gene body
matrix<-read.table("termgene_compare_matrix.mat.gz",skip=1,header=F,sep="\t")
matrix[sapply(matrix, is.numeric)] <- lapply(matrix[sapply(matrix, is.numeric)], function(x) {
  x[is.nan(x)] <- 0
  return(x)
})
# TSS 7:26 
# body 27:46 
# TES 47:66 

##scale 
body_Scale<-rowMeans(matrix[,30:43])
matrix[,7:66]<-(matrix[,7:66]-body_Scale)
matrix[] <- lapply(matrix, function(x) {
  x[is.na(x) | is.infinite(x)] <- 0
  return(x)
})

write.table(matrix,"Norm_TSS200_TES200_termgene_dcr1_WT_merge_add_compare_matrix.mat",quote=F,col.names=F,row.names=F,sep="\t")

###TES ＋500 boxplot 【47：56】
cotrol_enrich_matrix<-matrix[,c(1:6,47:56)]
cotrol_enrich_matrix[is.na(cotrol_enrich_matrix)]<-0
cotrol_enrich_mean<-data.frame(gene=cotrol_enrich_matrix[,4],CPM_sum=rowSums(cotrol_enrich_matrix[,7:16]))
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
        ylab = "log2(dcr1s/WTR2) Retention Index")+rotate_x_text(11)
dev.off()
write.table(cotrol_enrich_mean,"./sup_table/supplemental_table_NET_retentionIndex.tsv",quote = F,row.names = F,sep="\t")

###
cotrol_enrich_matrix<-matrix[,c(1:6,37:56)]
write.table(cotrol_enrich_matrix,"TES500_Norm_TSS200_TES200_termgene_dcr1_WT_merge_add_compare_matrix.mat",quote=F,col.names=F,row.names=F,sep="\t")

##shell 
#### You need to create a header for the matrix by yourself.
#for example
cat > TES500_termgene_dcr1_WT_merge_add_compare.header
# @{"upstream":[500],"downstream":[500],"body":[0],"bin size":[50],"ref point":["TES"],"verbose":false,"bin avg type":"mean","missing data as zero":true,"min threshold":null,"max threshold":500.0,"scale":1,"skip zeros":true,"nan after end":false,"proc number":15,"sort regions":"keep","sort using":"mean","unscaled 5 prime":[0],"unscaled 3 prime":[0],"group_labels":["dcr1_terminated_genes.bed","dcr1_expression_matched_genes.bed","pombe_coding_minus_dcr1_matched_terminated.bed"],"group_boundaries":[0,207,585,4096],"sample_labels":["dcr1_WT_merge_add_compare_1_skip0"],"sample_boundaries":[0,20]}

cat TES500_termgene_dcr1_WT_merge_add_compare.header TES500_Norm_TSS200_TES200_termgene_dcr1_WT_merge_add_compare_matrix.mat > TES500_Norm_TSS200_TES200_Termgene_dcr1_WT_merge_add_compare_matrix.mat
gzip TES500_Norm_TSS200_TES200_Termgene_dcr1_WT_merge_add_compare_matrix.mat
rm TES500_Norm_TSS200_TES200_termgene_dcr1_WT_merge_add_compare_matrix.mat 

plotHeatmap -m TES500_Norm_TSS200_TES200_Termgene_dcr1_WT_merge_add_compare_matrix.mat.gz -out ./per_TES500_Norm_TSS200_TES200_Termgene_merge_add_dcr1_WT_compare_heatmap.pdf --colorMap RdBu_r --dpi 300 --boxAroundHeatmaps no --heatmapHeight 9 --heatmapWidth 5 --missingDataColor 1 --plotTitle 'Body normalization of NET-seq' --samplesLabel termindex_dcr1_WT_merge_add --regionsLabel term match other --perGroup

###################################Figure S1 G
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


###################################Figure S1 H
computeMatrix scale-regions -S ${RNAPII_Compare_dir} -R ${dir}dcr1_terminated_genes.bed -b 1000 -a 1000 --regionBodyLength 1000 --skipZeros -bs 50 -p 15 -o ./computeMatrix/${RNAPII_Compare}_term_matrix.mat.gz --maxThreshold 500 --missingDataAsZero
plotHeatmap -m ./computeMatrix/${RNAPII_Compare}_term_matrix.mat.gz  -out ./computeMatrix/${RNAPII_Compare}_termgene_heatmap.pdf --colorMap RdBu_r --dpi 300 --boxAroundHeatmaps no --heatmapHeight 10 --heatmapWidth 5 --missingDataColor 1 --samplesLabel ${RNAPII_Compare}  --regionsLabel term 

###################################Figure S1 I
computeMatrix scale-regions -S ${sRNA_dir} -R ${dir}dcr1_terminated_genes.bed -b 1000 -a 1000 --regionBodyLength 1000 --skipZeros -bs 50 -p 15 -o ./computeMatrix/${sRNA}_term_matrix.mat.gz --maxThreshold 500 --missingDataAsZero
plotHeatmap -m ./computeMatrix/${sRNA}_term_matrix.mat.gz  -out ./computeMatrix/${sRNA}_termgene_heatmap.pdf --colorMap RdBu_r --dpi 300 --boxAroundHeatmaps no --heatmapHeight 10 --heatmapWidth 5 --missingDataColor 1 --samplesLabel ${sRNA}  --regionsLabel term --zMax 2.5 --zMin -2.5

###################################Figure S1 J
computeMatrix reference-point --referencePoint TES -S $sRNA_dir -R ${term_dir}dcr1_terminated_genes.bed ${term_dir}dcr1_expression_matched_genes_chr.bed ${term_dir}pombe_coding_minus_dcr1_matched_terminated.bed.bed -b 200 -a 200 --skipZeros -bs 50 -p 15 -o ./computeMatrix/TES200_sRNA_WT_dcr1_minus_matrix.mat.gz --minThreshold -100 --maxThreshold 100
 
#boxplot
 
enrich_matrix<-read.table("./computeMatrix/TES200_sRNA_WT_dcr1_minus_matrix.mat.gz",skip=1,header=F,sep="\t")
# ,"pombe_coding_minus_dcr1_matched_terminated.bed.bed"],"group_boundaries":[0,205,571,3744],"sample_labels":["dg21_dg690_minus"],"sample_boundaries":[0,8]}

cotrol_enrich_matrix<-enrich_matrix[,1:14]
cotrol_enrich_matrix[is.na(cotrol_enrich_matrix)]<-0
cotrol_enrich_mean<-data.frame(gene=cotrol_enrich_matrix[,4],CPM_sum=rowSums(cotrol_enrich_matrix[,7:14]))
cotrol_enrich_mean$Features<-"Dcr1-terminated"
cotrol_enrich_mean$Features[206:571]<-"Expression-matched"
cotrol_enrich_mean$Features[572:3744]<-"others"

my_comparisons <- list( c("Dcr1-terminated", "Expression-matched"), c("Dcr1-terminated", "others"), c("Expression-matched", "others") )
p<-ggboxplot(cotrol_enrich_mean, x = "Features", y = "CPM_sum",
          color = "Features", outlier.shape = NA,palette = "jco",bxp.errorbar=T)+ 
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")+coord_cartesian(ylim =  c(-2, 20))
 

  ggpar(p,
        legend = "right",
        xlab = "Gene type",
        ylab = "sRNA WT RPM - dcr1Δ RPM)")+rotate_x_text(11)

write.table(cbind(cotrol_enrich_matrix[,1:3],cotrol_enrich_mean),"./computeMatrix/boxplot/sup_table_TES200_sRNA_WT_dcr1_boxplot.tsv",quote = F,row.names = F,sep="\t")
