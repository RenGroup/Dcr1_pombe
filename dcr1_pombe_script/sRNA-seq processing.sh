#sRNA data preprocess
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

mkdir -p $read_dir/trim_data
mkdir -p $read_dir/bowtie2
mkdir -p $read_dir/bowtie2/bw

java -jar ${trim}/trimmomatic-0.39.jar SE -threads 15 -phred33 $read_dir/Rawdata/$name.fastq $read_dir/trim_data/${name}.fq.gz SLIDINGWINDOW:4:15 MINLEN:15 > $read_dir/trim_data/$name.trim.log 2>&1

bowtie2 -x ${bowtie2_index} -r $read_dir/trim_data/${name}.fq.gz --no-unal -p $threads -q -S $read_dir/bowtie2/${name}.sam > $read_dir/bowtie2/${name}_log.txt 2>&1
samtools view -Shb $read_dir/bowtie2/${name}.sam > $read_dir/bowtie2/${name}.bam
#sort
samtools sort -m 4G $read_dir/bowtie2/${name}.bam -o $read_dir/bowtie2/${name}.sort.bam

#bam index
samtools index $read_dir/bowtie2/${name}.sort.bam
#rmdup
samtools rmdup $read_dir/bowtie2/${name}.sort.bam $read_dir/bowtie2/${name}.sort.rmdup.bam

#pileup
export total_reads=$(samtools view -c $read_dir/bowtie2/${name}.sort.rmdup.bam)
export norm_factor=$(echo 'scale=5;1000000/'$total_reads | bc)

bedtools genomecov -bg -scale $norm_factor -ibam $read_dir/bowtie2/${name}.sort.rmdup.bam -g $genome -strand + | sort -k 1,1 -k 2,2n > $read_dir/bowtie2/bw/$name.rmdup_plus.bg
bedtools genomecov -bg -scale $norm_factor -strand "-" -ibam $read_dir/bowtie2/${name}.sort.rmdup.bam -g $genome | sort -k 1,1 -k 2,2n > $read_dir/bowtie2/bw/$name.rmdup_minus.bg

/work/home/path/soft/bedgraphtobigwig/bedGraphToBigWig $read_dir/bowtie2/bw/$name.rmdup_plus.bg ${genome} $read_dir/bowtie2/bw/$name.rmdup_plus.bw
/work/home/path/soft/bedgraphtobigwig/bedGraphToBigWig $read_dir/bowtie2/bw/$name.rmdup_minus.bg ${genome} $read_dir/bowtie2/bw/$name.rmdup_minus.bw

