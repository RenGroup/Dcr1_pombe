###################################Figure S4A DNAP Ori
computeMatrix reference-point --referencePoint center -S $WT_compare -R ${ori} -b 1000 -a 1000 --skipZeros -bs 50 -p 15 -o ./computeMatrix/ori_WT_DNAPD_DNAPE_IP_compare_matrix.mat.gz --maxThreshold 1000 --missingDataAsZero
plotHeatmap -m ./computeMatrix/ori_WT_DNAPD_DNAPE_IP_compare_matrix.mat.gz -out ./computeMatrix/ori_center_WT_DNAPD_DNAPE_IP_compare_heatmap.pdf --colorMap RdBu_r --dpi 300 --boxAroundHeatmaps no --heatmapHeight 9 --heatmapWidth 5 --missingDataColor 1 --regionsLabel replication_originals --samplesLabel WT_DNAPD_DNAPE_IP

###################################Figure S4C Rqh1 boxplot
enrich_matrix<-read.table("./computeMatrix/compare_Gene_matrix.mat.gz",skip=1,header=F,sep="\t")
#"group_labels":["dcr1_terminated_genes.bed","dcr1_expression_matched_genes.bed","pombe_coding_minus_dcr1_matched_terminated.bed"],"group_boundaries":[0,207,585,4096],"sample_labels":["Rqh1_dcr1_WT_mean_compare"],"sample_boundaries":[0,20]}

cotrol_enrich_matrix<-enrich_matrix[,1:26]
cotrol_enrich_matrix[is.na(cotrol_enrich_matrix)]<-0

cotrol_enrich_mean<-data.frame(gene=cotrol_enrich_matrix[,4],CPM_sum=rowSums(cotrol_enrich_matrix[,7:26]))
 

cotrol_enrich_mean$Features<-"Dcr1-terminated"
cotrol_enrich_mean$Features[208:585]<-"Expression-matched"
cotrol_enrich_mean$Features[586:4096]<-"others"

my_comparisons <- list( c("Dcr1-terminated", "Expression-matched"), c("Dcr1-terminated", "others"), c("Expression-matched", "others") )
p<-ggboxplot(subset(cotrol_enrich_mean,CPM_sum>(-1)), x = "Features", y = "CPM_sum",
          color = "Features", palette = "jco",outlier.shape = NA,bxp.errorbar=T)+ 
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
  
  pdf("./boxplot/Rqh1_dcr1_WT_compare_Sum_boxplot.pdf",width=4,height=4)
# Visualize: Specify the comparisons you want
  ggpar(p,
        legend = "right",
        xlab = "Gene type",
        ylab = "Rqh1 RPM log2(dcr1Δ/WT)")+rotate_x_text(11)
  # Add global p-value
dev.off()

write.table(cotrol_enrich_mean,"./boxplot/sup_table_Rqh1_ChIP_compare_gene.txt",quote=F,col.names=T,row.names=F,sep="\t")