###exprt
export name=$1
export read_dir=$2
export threads=$3

export PATH=/work/home/path/soft/Anaconda3/anaconda3/envs/software/bin:$PATH
export trim=/work/home/path/soft/Trimmomatic-0.39

export hisat2_index="/work/home/path/ref/yeast/ASM294v2.23/hisat2_ASM294v2.23/hisat2_S_pombe.2.23"
export genome=/work/home/path/ref/yeast/pombe_ASM294v2.23.genome

mkdir -p ${read_dir}/Rawdata/qc
mkdir -p ${read_dir}/trim_data
mkdir -p ${read_dir}/trim_data/qc
mkdir -p ${read_dir}/hisat2
mkdir -p ${read_dir}/hisat2/fwd_rev
mkdir -p ${read_dir}/hisat2/flagstat
mkdir -p ${read_dir}/hisat2/bw
mkdir -p ${read_dir}/hisat2/picard/
mkdir -p ${read_dir}/hisat2/fwd_rev/
mkdir -p ${read_dir}/hisat2/fwd_rev/bw

#qc
fastqc --noextract -t $threads -f fastq $read_dir/Rawdata/$name*fq.gz -o $read_dir/Rawdata/qc/

###trim

java -jar ${trim}/trimmomatic-0.39.jar SE -threads $threads -phred33 $read_dir/Rawdata/${name}.fq.gz ${read_dir}/trim_data/$name.trim.fq.gz SLIDINGWINDOW:4:15 MINLEN:20 ILLUMINACLIP:/work/home/path/ref/polyA_adapter.fa:2:30:10 > $read_dir/trim_data/$name.trim.log 2>&1 
cutadapt -a "A{10}" -g "T{10}" -o ${read_dir}/trim_data/$name.trim.cutadapt.fq.gz ${read_dir}/trim_data/$name.trim.fq.gz -q 20 --minimum-length=20 > $read_dir/trim_data/$name.trim.cutadapt.log

fastqc ${read_dir}/trim_data/$name.trim.fq.gz -o ${read_dir}/trim_data/qc/ -t $threads
fastqc ${read_dir}/trim_data/$name.trim.cutadapt.fq.gz -o ${read_dir}/trim_data/qc/ -t $threads

hisat2 -p $threads --dta -x ${hisat2_index} -U ${read_dir}/trim_data/$name.trim.cutadapt.fq.gz -S ${read_dir}/hisat2/$name.trim.sam > ${read_dir}/hisat2/$name.trim_map.log 2>&1
samtools view -Shb ${read_dir}/hisat2/$name.trim.sam | samtools sort -@ 12 -m 4G > ${read_dir}/hisat2/$name.trim.bam

samtools index ${read_dir}/hisat2/$name.trim.bam

##
samtools view -hb ${read_dir}/hisat2/$name.trim.bam I II III | samtools view -hb -F 0x04 | samtools sort -@ $threads -m 2G > ${read_dir}/hisat2/$name.Pombe.trim.bam
samtools flagstat ${read_dir}/hisat2/$name.Pombe.trim.bam > ${read_dir}/hisat2/flagstat/${name}.Pombe.trim.raw_rmdup_flagstat.txt

##pcr
picard MarkDuplicates \
  QUIET=true INPUT=${read_dir}/hisat2/$name.Pombe.trim.bam OUTPUT=${read_dir}/hisat2/$name.Pombe.trim.rmdup.bam METRICS_FILE=$read_dir/hisat2/picard/$name.trim.sort.metrics \
  REMOVE_DUPLICATES=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=/work/home/path/tmp
samtools flagstat ${read_dir}/hisat2/$name.Pombe.trim.rmdup.bam > $read_dir/hisat2/flagstat/${name}.Pombe.trim.rmdup_flagstat.txt

##uniq mapping reads
samtools view -h ${read_dir}/hisat2/$name.Pombe.trim.rmdup.bam |awk 'BEGIN {header=1} /^@/ {print $0} !/^@/ && $0 ~ /NH:i:1/ && $0 !~ /ZS:i/ {print $0}' > ${read_dir}/hisat2/$name.Pombe.trim.rmdup.uniq.sam

samtools view -Shb ${read_dir}/hisat2/$name.Pombe.trim.rmdup.uniq.sam | samtools sort -@ 5 -O bam -o - > ${read_dir}/hisat2/$name.Pombe.trim.rmdup.uniq.bam

samtools flagstat ${read_dir}/hisat2/$name.Pombe.trim.rmdup.uniq.bam > $read_dir/hisat2/flagstat/${name}.Pombe.trim.sort.rmdup.uniq_flagstat.txt
rm ${read_dir}/hisat2/$name.Pombe.trim.rmdup.uniq.sam

#pileup
echo mapping data pileup start!
factor=`echo "scale=5;1000000/$(samtools view -@ 8 -c ${read_dir}/hisat2/$name.Pombe.trim.rmdup.uniq.bam)" | bc`
bedtools genomecov -bg -scale ${factor} -ibam ${read_dir}/hisat2/$name.Pombe.trim.rmdup.uniq.bam | sort -k1,1 -k2,2n > $read_dir/hisat2/bw/${name}.Pombe.trim.rmdup.uniq.bg
/work/home/path/soft/bedgraphtobigwig/bedGraphToBigWig $read_dir/hisat2/bw/${name}.Pombe.trim.rmdup.uniq.bg ${genome} $read_dir/hisat2/bw/${name}.Pombe.trim.rmdup.uniq.bw
echo mapping data pileup end!

###Fwd_Rev
Pair="${name: -2}"
if [[ "$Pair" == "R1" ]]; then
samtools view -b -f 16 ${read_dir}/hisat2/$name.Pombe.trim.rmdup.uniq.bam > ${read_dir}/hisat2/fwd_rev/${name}_rev.trim.rmdup.uniq.sort.bam
samtools view -b -F 16 ${read_dir}/hisat2/$name.Pombe.trim.rmdup.uniq.bam > ${read_dir}/hisat2/fwd_rev/${name}_fwd.trim.rmdup.uniq.sort.bam
fi

if [[ "$Pair" == "R2" ]]; then
samtools view -b -f 16 ${read_dir}/hisat2/$name.Pombe.trim.rmdup.uniq.bam > ${read_dir}/hisat2/fwd_rev/${name}_fwd.trim.rmdup.uniq.sort.bam
samtools view -b -F 16 ${read_dir}/hisat2/$name.Pombe.trim.rmdup.uniq.bam > ${read_dir}/hisat2/fwd_rev/${name}_rev.trim.rmdup.uniq.sort.bam
fi
#bam-bg-bw scale(total bam)
genome=/work/home/path/ref/yeast/pombe_ASM294v2.23M.genome
factor=`echo "scale=5;1000000/$(samtools view -@ 8 -c ${read_dir}/hisat2/$name.Pombe.trim.rmdup.uniq.bam )" | bc`
bedtools genomecov -bg -scale ${factor} -ibam ${read_dir}/hisat2/fwd_rev/${name}_fwd.trim.rmdup.uniq.sort.bam | sort -k1,1 -k2,2n > ${read_dir}/hisat2/fwd_rev/bw/${name}_fwd.trim.rmdup.uniq.sort.bg
/work/home/path/soft/bedgraphtobigwig/bedGraphToBigWig ${read_dir}/hisat2/fwd_rev/bw/${name}_fwd.trim.rmdup.uniq.sort.bg ${genome} ${read_dir}/hisat2/fwd_rev/bw/${name}_fwd.trim.rmdup.uniq.sort.bw
bedtools genomecov -bg -scale ${factor} -ibam ${read_dir}/hisat2/fwd_rev/${name}_rev.trim.rmdup.uniq.sort.bam | sort -k1,1 -k2,2n > ${read_dir}/hisat2/fwd_rev/bw/${name}_rev.trim.rmdup.uniq.sort.bg
/work/home/path/soft/bedgraphtobigwig/bedGraphToBigWig ${read_dir}/hisat2/fwd_rev/bw/${name}_rev.trim.rmdup.uniq.sort.bg ${genome} ${read_dir}/hisat2/fwd_rev/bw/${name}_rev.trim.rmdup.uniq.sort.bw
