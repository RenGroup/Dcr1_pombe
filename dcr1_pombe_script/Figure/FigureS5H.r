###################################Figure S5H Rad51 ChIP-seq boxplot
###shell
computeMatrix reference-point --referencePoint TES -S ${Rad51_wt_dcr1} -R ${dir}dcr1_terminated_genes.bed ${dir}dcr1_expression_matched_genes.bed ${dir}pombe_coding_minus_dcr1_matched_terminated.bed -b 500 -a 500 --skipZeros -bs 50 -p 15 -o ./computeMatrix/Rad51_TES500_termgene_compare_matrix.mat.gz --maxThreshold 500

###R
enrich_matrix<-read.table("./Rad51_TES500_termgene_compare_matrix.mat.gz",skip=1,header=F,sep="\t")
 
# ,"group_labels":["dcr1_terminated_genes.bed","dcr1_expression_matched_genes.bed","pombe_coding_minus_dcr1_matched_terminated.bed"],"group_boundaries":[0,207,585,4096],"sample_labels":["DG21_DG690_mean_compare"],"sample_boundaries":[0,20]}

cotrol_enrich_matrix<-enrich_matrix[,1:26]
cotrol_enrich_matrix[is.na(cotrol_enrich_matrix)]<-0
cotrol_enrich_mean<-data.frame(gene=cotrol_enrich_matrix[,4],RPM_mean=rowSums(cotrol_enrich_matrix[,7:26]))
cotrol_enrich_mean$Features<-"Dcr1-terminated"
cotrol_enrich_mean$Features[208:585]<-"Expression-matched"
cotrol_enrich_mean$Features[586:4096]<-"others"

my_comparisons <- list( c("Dcr1-terminated", "Expression-matched"), c("Dcr1-terminated", "others"), c("Expression-matched", "others") )
p<-ggboxplot(cotrol_enrich_mean, x = "Features", y = "RPM_mean",
          color = "Features", palette = "jco",outlier.shape = NA,bxp.errorbar=T)+ coord_cartesian(ylim =  c(-7, 40))+ 
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")

  ggpar(p,
        legend = "right",
        xlab = "Gene type",
        ylab = "Rad51 RPM log2(WT/dcr1Δ)")+rotate_x_text(11)
