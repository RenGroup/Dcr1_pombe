#RNAPII ChIP-seq data preprocess
# $1  file name :filename(_1.fq.gz )
# $2 : data big file path(out of Rawdata)
# $3:threads
#threads:15

export name=$1
export read_dir=$2
export threads=$3
export trim=/work/home/path/soft/Trimmomatic-0.39
export PATH=/work/home/path/soft/bowtie2-2.4.2-linux-x86_64/:$PATH
export bowtie2_index="/work/home/path/ref/yeast/pomb_bowtie2/ASM294v2.23"
export genome=/work/home/path/ref/yeast/pombe_ASM294v2.23.genome

mkdir $read_dir/trim_data
mkdir $read_dir/trim_data/qc
mkdir $read_dir/bowtie2
mkdir $read_dir/bowtie2/bw

###trim
echo Rawdata Trim Start!
java -jar ${trim}/trimmomatic-0.39.jar PE -threads $threads -phred33 $read_dir/Rawdata/${name}_1.fq.gz $read_dir/Rawdata/${name}_2.fq.gz -baseout $read_dir/trim_data/$name.trim.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 ILLUMINACLIP:/work/home/path/ref/adapter_sum.fa:2:30:10 > $read_dir/trim_data/$name.trim.log 2>&1
echo Rawdata Trim END!

echo Trim data qc start!
fastqc -o $read_dir/trim_data/qc/ --noextract -t $threads -f fastq $read_dir/trim_data/$name.trim*.fq.gz

echo Trim data qc end!

#bowtie2 mapping
echo Trim data mapping start!
bowtie2 -x ${bowtie2_index} -1 $read_dir/trim_data/$name.trim_1P.fq.gz -2 $read_dir/trim_data/$name.trim_2P.fq.gz -p $threads --phred33 --local --no-unal -q -S $read_dir/bowtie2/$name.trim.sam > $read_dir/bowtie2/$name.trim_map.log 2>&1
samtools view -Shb $read_dir/bowtie2/$name.trim.sam > $read_dir/bowtie2/$name.trim.bam
echo Trim data mapping END!

echo mapping data sorting start!
samtools sort -@ 12 -m 4G $read_dir/bowtie2/$name.trim.bam -o $read_dir/bowtie2/$name.trim.sort.bam
echo mapping data rmdup start!
samtools rmdup $read_dir/bowtie2/$name.trim.sort.bam $read_dir/bowtie2/$name.trim.sort.rmdup.bam
echo mapping data index start!
samtools index $read_dir/bowtie2/$name.trim.sort.rmdup.bam
rm $read_dir/bowtie2/$name.trim.sam $read_dir/bowtie2/$name.trim.bam $read_dir/bowtie2/$name.trim.sort.bam

#RPM scale->bw
echo mapping data pileup start!
factor=`echo "scale=5;1000000/$(samtools view -@ 8 -c $read_dir/bowtie2/$name.trim.sort.rmdup.bam)" | bc`
bedtools genomecov -bg -scale ${factor} -ibam $read_dir/bowtie2/$name.trim.sort.rmdup.bam | sort -k1,1 -k2,2n > $read_dir/bowtie2/bw/$name.trim.sort.rmdup.bg
/work/home/path/soft/bedgraphtobigwig/bedGraphToBigWig $read_dir/bowtie2/bw/$name.trim.sort.rmdup.bg ${genome} $read_dir/bowtie2/bw/$name.trim.sort.rmdup.bw
echo mapping data pileup end!
