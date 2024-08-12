####Fig 3C #tDNA dotplot RNAPIIcompare ssDRIPWT

computeMatrix scale-regions -S $RNAPII_dir $ssDrip_WT_dir $RIP_merge_dir -R $tRNA_dir -b 0 -a 0 --regionBodyLength 500 --skipZeros -bs 10 -p 10 -o ./tDNA_RNAPIIcomp_ssDRIPwt_RIP_noupd_matrix.mat.gz

enrich_matrix<-read.table("./tDNA_RNAPIIcomp_ssDRIPwt_RIP_noupd_matrix.mat.gz",skip=1,header=F,sep="\t")
#"group_labels":["genes"],"group_boundaries":[0,171],"sample_labels":["polII_compare","DG20_21IP_input","Dcr1_RIP_merge12_IP_input"],"sample_boundaries":[0,50,100,150]}

enrich_matrix[is.na(enrich_matrix)]<-0
control_mean<-rowMeans(enrich_matrix[,7:56])
treat_mean<-rowMeans(enrich_matrix[,57:106])

RIP_mean<-apply(enrich_matrix[,107:156],1,max)
plot_matrix<-data.frame(loc=paste(enrich_matrix[,1],enrich_matrix[,2],enrich_matrix[,3],sep=":"),treat_mean_CPM=as.numeric(treat_mean),control_mean_CPM=as.numeric(control_mean),RIP_mean=as.numeric(RIP_mean))
plot_matrix<-plot_matrix[order(plot_matrix$RIP_mean,decreasing=F),]


pdf("./tDNA_RNAPIIcomp_ssDRIPwt_RIPmax_color_dotplot.pdf",width=5,height=4)
ggplot(plot_matrix,aes(x=control_mean_CPM,y=treat_mean_CPM))+
xlim(-2.5,2.5)+ylim(-2.5,5)+
geom_point(aes(color=RIP_mean))+
geom_hline(yintercept = 0 ,linetype = "dashed") +
geom_vline(xintercept = 0,linetype = "dashed") +
scale_colour_gradient2(low = "gray", high = "red",mid="#fddbc8", midpoint = 0.5)+
  labs(x = "RNAPII log2(dcr1Δ/WT)", y = "WT ssDRIP log2(IP/Input)",colour ="Dcr1 RIP log2(IP/Input)")+theme_bw()+
  theme(panel.grid.major = element_blank(),   
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"))  
dev.off()


####fig3D#########RNAPII compare Cumulative tDNA

computeMatrix scale-regions -S $RNAPII_WT_dir $RNAPII_Dcr1_dir -R $tRNA_dir -b 0 -a 0 --regionBodyLength 500 --skipZeros -bs 10 -p 10 -o ./RNAPII_tDNA_noupd_dcr1d_WTenrich_matrix.mat.gz
#R
enrich_matrix<-read.table("./RNAPII_tDNA_noupd_dcr1d_WTenrich_matrix.mat.gz",skip=1,header=F,sep="\t")
# "group_labels":["genes"],"group_boundaries":[0,171],"sample_labels":["WT_pS2_rep1IP_input","dcr1_pS2_rep2IP_input"],"sample_boundaries":[0,100,200]}
enrich_matrix[is.na(enrich_matrix)]<-0
control_mean<-rowMeans(enrich_matrix[,7:106])
treat_mean<-rowMeans(enrich_matrix[,107:206])

plot_matrix<-data.frame(loc=paste(enrich_matrix[,1],enrich_matrix[,2],enrich_matrix[,3],sep=":"),dcr1d=as.numeric(treat_mean),WT=as.numeric(control_mean))
plot_matrix<-melt(plot_matrix)
colnames(plot_matrix)<-c("loc","feature","value")

p<-ggplot(plot_matrix,
       aes(
         x=value,
         color=feature
         ))+
  stat_ecdf( geom="smooth", se = F )+theme_bw()+
  theme(panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

pdf("./RNAPII_tDNA_Cumulative_240606.pdf",width=6,height=4)
  ggpar(p,
        legend = "right",
        xlab = "RNAPII tDNA log2(dcr1Δ/WT)",
        ylab = "Cumulative distribution")
dev.off() #()

####fig 3E ######DNAPe Cumulative tDNA
computeMatrix scale-regions -S $WT_dir $Dcr1_dir -R $tRNA_dir -b 0 -a 0 --regionBodyLength 500 --skipZeros -bs 10 -p 10 -o ./DNAPe_tDNA_noupd_dcr1d12m_WT12_matrix.mat.gz

enrich_matrix<-read.table("./DNAPe_tDNA_noupd_dcr1d12m_WT12_matrix.mat.gz",skip=1,header=F,sep="\t")
# "group_labels":["genes"],"group_boundaries":[0,171],"sample_labels":["FY14792_sr_1IP_2IP_1INPUT_enrich","JRY224_sr_1IP_2IP_1INPUT_enrich"],"sample_boundaries":[0,100,200]}
enrich_matrix[is.na(enrich_matrix)]<-0
control_mean<-rowMeans(enrich_matrix[,7:106])
treat_mean<-rowMeans(enrich_matrix[,107:206])
plot_matrix<-data.frame(loc=paste(enrich_matrix[,1],enrich_matrix[,2],enrich_matrix[,3],sep=":"),dcr1d=as.numeric(treat_mean),WT=as.numeric(control_mean))
plot_matrix<-melt(plot_matrix)
colnames(plot_matrix)<-c("loc","feature","value")
p<-ggplot(plot_matrix,
       aes(
         x=value,
         color=feature
         ))+
  stat_ecdf()

pdf("./DNAPe_tDNA_Cumulative.pdf",width=6,height=4)
# Visualize: Specify the comparisons you want
  ggpar(p,
        legend = "right",
        xlab = "DNAPe tDNA log2(IP/Input)",
        ylab = "Cumulative distribution")
dev.off()

####fig3F#########DNAPd  Cumulative tDNA
computeMatrix scale-regions -S $DNAPd_WT_dir $DNAPd_dcr1_dir -R $tRNA_dir -b 0 -a 0 --regionBodyLength 500 --skipZeros -bs 10 -p 10 -o ./dotplot/DNAPd_tRNA_noupd_dcr1d_WTenrich_240624_matrix.mat.gz
# "group_labels":["genes"],"group_boundaries":[0,171],"sample_labels":["21-8cdc1_WT_WGmean_IP_input","DNAPD_dcr1d_IP24m_INPUT240620m_IP_input"],"sample_boundaries":[0,50,100]}
# R
enrich_matrix<-read.table("./dotplot/DNAPd_tRNA_noupd_dcr1d_WTenrich_240624_matrix.mat.gz",skip=1,header=F,sep="\t")
enrich_matrix[is.na(enrich_matrix)]<-0
control_mean<-rowMeans(enrich_matrix[,7:56])
treat_mean<-rowMeans(enrich_matrix[,57:106])
plot_matrix<-data.frame(loc=paste(enrich_matrix[,1],enrich_matrix[,2],enrich_matrix[,3],sep=":"),dcr1d=as.numeric(treat_mean),WT=as.numeric(control_mean))
plot_matrix<-melt(plot_matrix)
colnames(plot_matrix)<-c("loc","feature","value")

p<-ggplot(plot_matrix,
       aes(
         x=value,
         color=feature
         ))+
  stat_ecdf( )

pdf("./DNAPd_tDNA_Cumulative_240624.pdf",width=6,height=4)
  ggpar(p,
        legend = "right",
        xlab = "DNAPd log2(IP/Input)",
        ylab = "Cumulative distribution")
dev.off()