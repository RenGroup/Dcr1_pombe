export name=$1
export read_dir=$2
export threads=$3

export trim=/work/home/path/soft/Trimmomatic-0.39
export star_genome="/work/home/path/ref/yeast/pombe_STAR_ASM294v2.23"
thread=15
export ref_dir="/work/home/path/ref/yeast/rRNA/pombe_rRNA.bed"
export genome=/work/home/path/ref/yeast/pombe_ASM294v2.23M.genome
export coding_gene=/work/home/path/ref/yeast/coding/pombe_coding.bed
trim=/work/home/zhangyzh/soft/Trimmomatic-0.39
#mkdir
mkdir -p $read_dir/trim_data
mkdir -p $read_dir/fastuniq
mkdir -p $read_dir/star_intrond
mkdir -p $read_dir/star_intrond/bw

#trim
java -jar ${trim}/trimmomatic-0.39.jar PE -threads $threads -phred33 $read_dir/Rawdata/${name}_1.fq.gz $read_dir/Rawdata/${name}_2.fq.gz -baseout $read_dir/trim_data/${name}_trim.fq.gz LEADING:3 TRAILING:3 HEADCROP:10 SLIDINGWINDOW:4:15 MINLEN:56 ILLUMINACLIP:/work/home/zhangyzh/ref/adapter_sum.fa:2:30:10 > $read_dir/trim_data/$name.trim.log 2>&1 
#fastqctofastuniq
cat > $read_dir/fastuniq/${name}_file.txt
${name}_1.fq
${name}_2.fq
^z
fastuniq -i $read_dir/fastuniq/${name}_file.txt -o $read_dir/fastuniq/${name}_1P_fastuniq.fq -p $read_dir/fastuniq/${name}_2P_fastuniq.fq

#star intron 2000
echo data mapping start!
STAR --runThreadN $thread --genomeDir $star_genome --readFilesIn $read_dir/fastuniq/${name}_1P_fastuniq.fq $read_dir/fastuniq/${name}_2P_fastuniq.fq --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 8 --outMultimapperOrder Random --alignEndsType EndToEnd --outFileNamePrefix $read_dir/star_intrond/${name}.trim.fastuniq.intr --limitBAMsortRAM 16176068412 --alignIntronMax 2000
echo data mapping END!

#derRNA
echo mapping data derbed start!
bedtools intersect -v -abam $read_dir/star_intrond/${name}.trim.fastuniq.intrAligned.sortedByCoord.out.bam -b ${ref_dir} > $read_dir/star_intrond/${name}.trim.fastuniq.intr.derbed.bam
echo mapping data derbed end!
echo mapping data sorting start!
samtools sort -@ 12 -m 4G $read_dir/star_intrond/${name}.trim.fastuniq.intr.derbed.bam -o $read_dir/star_intrond/${name}.trim.fastuniq.intr.derbed.sort.bam
samtools index $read_dir/star_intrond/${name}.trim.fastuniq.intr.derbed.sort.bam
rm $read_dir/star_intrond/${name}.trim.fastuniq.intr.derbed.bam


#pileup
echo mapping data pileup start!
factor=`echo "scale=5;1000000/$(samtools view -@ 8 -c $read_dir/star_intrond/${name}.trim.fastuniq.intr.derbed.sort.bam)" | bc`
bedtools genomecov -bg -scale ${factor} -ibam $read_dir/star_intrond/${name}.trim.fastuniq.intr.derbed.sort.bam | sort -k1,1 -k2,2n > $read_dir/star_intrond/bw/${name}.trim.fastuniq.intr.derbed.sort.bg
/work/home/path/soft/bedgraphtobigwig/bedGraphToBigWig $read_dir/star_intrond/bw/${name}.trim.fastuniq.intr.derbed.sort.bg ${genome} $read_dir/star_intrond/bw/${name}.trim.fastuniq.intr.derbed.sort.bw
echo mapping data pileup end!

