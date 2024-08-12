#ssDrip data preprocess
# $1  file name
# $2 : data big file path(out of Rawdata)
# $3:threads
#threads:15

###exprt
export name=$1
export read_dir=$2
export threads=$3


export PATH=/work/home/path/soft/Anaconda3/anaconda3/envs/software/bin:$PATH
export trim=/work/home/path/soft/Trimmomatic-0.39
export PATH=/work/home/path/soft/bowtie2-2.4.2-linux-x86_64/:$PATH

export bowtie2_index="/work/home/path/ref/yeast/pomb_bowtie2/ASM294v2.23"
export genome=/work/home/path/ref/yeast/pombe_ASM294v2.23.genome
export bed_dir="work/home/path/ref/yeast/coding/pombe_coding.bed"


#mkdir
mkdir -p $read_dir/trim_data
mkdir -p $read_dir/bowtie2
mkdir -p $read_dir/bowtie2/fwd_rev
mkdir -p $read_dir/bowtie2/fwd_rev/bw
mkdir -p $read_dir/bowtie2/bw

###trim
echo Rawdata Trim Start!
java -jar ${trim}/trimmomatic-0.39.jar PE -threads $threads -phred33 $read_dir/Rawdata/${name}_R1.fq.gz $read_dir/Rawdata/${name}_R2.fq.gz -baseout $read_dir/trim_data/$name.trim.fq.gz LEADING:3 TRAILING:3 HEADCROP:15 SLIDINGWINDOW:4:15 MINLEN:66 ILLUMINACLIP:/work/home/path/ref/adapter_sum.fa:2:30:10 > $read_dir/trim_data/$name.trim.log 2>&1 
echo Rawdata Trim END!


# bowtie2 mapping
echo Trim data mapping start!
bowtie2 -x ${bowtie2_index} --phred33 -q --no-unal --local -1 $read_dir/trim_data/$name.trim_1P.fq.gz -2 $read_dir/trim_data/$name.trim_2P.fq.gz -p $threads -S $read_dir/bowtie2/$name.trim.sam > $read_dir/bowtie2/$name.trim_map.log 2>&1
samtools view -Shb $read_dir/bowtie2/$name.trim.sam > $read_dir/bowtie2/$name.trim.bam

#sort&rmdup
samtools sort -@ $threads -m 5G $read_dir/bowtie2/$name.trim.bam -o $read_dir/bowtie2/$name.trim.sort.bam
samtools rmdup $read_dir/bowtie2/$name.trim.sort.bam $read_dir/bowtie2/$name.trim.sort.rmdup.bam
#rmdup flagstat
samtools flagstat $read_dir/bowtie2/$name.trim.sort.rmdup.bam > $read_dir/bowtie2/flagstat/${name}.trim.sort.rmdup_flagstat.txt

#deM

samtools view -@ $threads -b $read_dir/bowtie2/$name.trim.sort.rmdup.q20.bam I II III | samtools sort -@ $threads -O bam -o - > $read_dir/bowtie2/$name.trim.sort.rmdup.deM.bam


echo mapping data index start!
samtools index $read_dir/bowtie2/$name.trim.sort.rmdup.deM.bam
rm $read_dir/bowtie2/$name.trim.sam $read_dir/bowtie2/$name.trim.bam $read_dir/bowtie2/$name.trim.sort.bam $read_dir/bowtie2/$name.trim.sort.rmdup.bam

#fwd_rev
#fwd
samtools view -b -f 128 -F 16 $read_dir/bowtie2/$name.trim.sort.rmdup.deM.bam > $read_dir/bowtie2/fwd_rev/$name.fwd1.bam
samtools index $read_dir/bowtie2/fwd_rev/$name.fwd1.bam
samtools view -b -f 80 $read_dir/bowtie2/$name.trim.sort.rmdup.deM.bam > $read_dir/bowtie2/fwd_rev/$name.fwd2.bam
samtools index $read_dir/bowtie2/fwd_rev/$name.fwd2.bam
samtools merge -f $read_dir/bowtie2/fwd_rev/$name.fwd.bam $read_dir/bowtie2/fwd_rev/$name.fwd1.bam $read_dir/bowtie2/fwd_rev/$name.fwd2.bam
samtools index $read_dir/bowtie2/fwd_rev/$name.fwd.bam
rm $read_dir/bowtie2/fwd_rev/$name.fwd1.bam $read_dir/bowtie2/fwd_rev/$name.fwd2.bam $read_dir/bowtie2/fwd_rev/$name.fwd1.bam.bai $read_dir/bowtie2/fwd_rev/$name.fwd2.bam.bai

#rev
samtools view -b -f 144 $read_dir/bowtie2/$name.trim.sort.rmdup.deM.bam > $read_dir/bowtie2/fwd_rev/$name.rev1.bam
samtools index $read_dir/bowtie2/fwd_rev/$name.rev1.bam
samtools view -b -f 64 -F 16 $read_dir/bowtie2/$name.trim.sort.rmdup.deM.bam > $read_dir/bowtie2/fwd_rev/$name.rev2.bam
samtools index $read_dir/bowtie2/fwd_rev/$name.rev2.bam
samtools merge -f $read_dir/bowtie2/fwd_rev/$name.rev.bam $read_dir/bowtie2/fwd_rev/$name.rev1.bam $read_dir/bowtie2/fwd_rev/$name.rev2.bam
samtools index $read_dir/bowtie2/fwd_rev/$name.rev.bam
rm $read_dir/bowtie2/fwd_rev/$name.rev1.bam $read_dir/bowtie2/fwd_rev/$name.rev2.bam $read_dir/bowtie2/fwd_rev/$name.rev1.bam.bai $read_dir/bowtie2/fwd_rev/$name.rev2.bam.bai

#RPM scale->bw
echo mapping data pileup start!
factor=`echo "scale=5;1000000/$(samtools view -@ $threads -c $read_dir/bowtie2/$name.trim.sort.rmdup.deM.bam)" | bc`

bedtools genomecov -bg -scale ${factor} -ibam $read_dir/bowtie2/$name.trim.sort.rmdup.deM.bam | sort -k1,1 -k2,2n > $read_dir/bowtie2/bw/$name.trim.sort.rmdup.deM.bg
/work/home/path/soft/bedgraphtobigwig/bedGraphToBigWig_Linux $read_dir/bowtie2/bw/$name.trim.sort.rmdup.deM.bg ${genome} $read_dir/bowtie2/bw/$name.trim.sort.rmdup.deM.bw

#fwd_rev pileup
bedtools genomecov -bg -scale ${factor} -ibam $read_dir/bowtie2/fwd_rev/$name.fwd.bam | sort -k1,1 -k2,2n > $read_dir/bowtie2/fwd_rev/bw/$name.fwd.bg
/work/home/path/soft/bedgraphtobigwig/bedGraphToBigWig_Linux $read_dir/bowtie2/fwd_rev/bw/$name.fwd.bg ${genome} $read_dir/bowtie2/fwd_rev/bw/$name.fwd.bw

bedtools genomecov -bg -scale ${factor} -ibam $read_dir/bowtie2/fwd_rev/$name.rev.bam | sort -k1,1 -k2,2n > $read_dir/bowtie2/fwd_rev/bw/$name.rev.bg
/work/home/path/soft/bedgraphtobigwig/bedGraphToBigWig_Linux $read_dir/bowtie2/fwd_rev/bw/$name.rev.bg ${genome} $read_dir/bowtie2/fwd_rev/bw/$name.rev.bw
