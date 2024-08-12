####FigS3A WT DNAPD/DNAPE in ori
computeMatrix reference-point --referencePoint center -S $WT_compare -R ${ori} -b 1000 -a 1000 --skipZeros -bs 50 -p 15 -o ./computeMatrix/ori_WT_DNAPD_DNAPE_IP_compare_matrix.mat.gz --maxThreshold 1000 --missingDataAsZero
plotHeatmap -m ./computeMatrix/ori_WT_DNAPD_DNAPE_IP_compare_matrix.mat.gz -out ./computeMatrix/ori_center_WT_DNAPD_DNAPE_IP_compare_heatmap.pdf --colorMap RdBu_r --dpi 300 --boxAroundHeatmaps no --heatmapHeight 9 --heatmapWidth 5 --missingDataColor 1 --regionsLabel replication_originals --samplesLabel WT_DNAPD_DNAPE_IP



####FIGS3B #DNAP compare

computeMatrix reference-point --referencePoint TES -S ./DNAPD_DNAPE_IP_minus.bw -R ${term_dir}dcr1_terminated_genes.bed ${term_dir}dcr1_expression_matched_genes.bed ${term_dir}pombe_coding_minus_dcr1_matched_terminated.bed -b 500 -a 500 --skipZeros -bs 10 -p 10 -o ./computeMatrix/DNAPD_DNAPE_IP_minus_compare_TES_matrix.mat.gz --maxThreshold 500

plotHeatmap -m ./computeMatrix/DNAPD_DNAPE_IP_minus_compare_TES_matrix.mat.gz -out ./computeMatrix/DNAPD_DNAPE_IP_minus_compare_TES_heatmap.pdf --dpi 300 --boxAroundHeatmaps no --heatmapHeight 9 --heatmapWidth 5 --missingDataColor 1 --regionsLabel term match others --samplesLabel "DNAPd_IP_compare-DNAPe_IP_compare"  --perGroup --colorMap RdBu_r
