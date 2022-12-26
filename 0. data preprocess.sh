# make directorys
mkdir rawdata cleandata bwa bigwig profile bedgraph macs2

# rename files
for i in `ls *_R1.fq.gz`; do x=${i/_R1.fq.gz/}; mv $i $x.R1.fastq.gz; done
for i in `ls *_R2.fq.gz`; do x=${i/_R2.fq.gz/}; mv $i $x.R2.fastq.gz; done

# generate Bam files
for i in `ls ../rawdata/*.R1.fastq.gz`; do x=${i/.R1.fastq.gz/}; x=${x/..\/rawdata\//}; echo "fastp -i $i -I ../rawdata/$x.R2.fastq.gz  -o ../cleandata/$x.r1.fq.gz -O ../cleandata/$x.r2.fq.gz --detect_adapter_for_pe -w 4 --compression 9 -h ../cleandata/$x.html -j ../cleandata/$x.json &&
bwa mem -M -t 4 /public-supool/home/longz/reference/cs/iwgsc1/bwa_index/cs ../cleandata/$x.r1.fq.gz ../cleandata/$x.r2.fq.gz > ../bwa/$x.sam &&
samtools view -bS -F 1804 -f 2 -q 30 ../bwa/$x.sam | samtools sort - | samtools rmdup -s - ../bwa/$x.raw.bam &&
java -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4 -Xmx8000M -jar ~/bin/picard.jar MarkDuplicates I=../bwa/$x.raw.bam O=../bwa/$x.rmdup.bam M=../bwa/$x.txt REMOVE_DUPLICATES=true
"; done > command.sh
sh command.sh

# call peaks

# 1. SEACR
for i in `ls ../bwa/*.rmdup.bam`; do x=${i/..\/bwa\//}; x=${x/.rmdup.bam/}; 
echo "bedtools bamtobed -bedpe -i $i | awk '\$1==\$4 && \$6-\$2 < 1000 {print \$1,\$2,\$6}' OFS =\"\t\"| sort -k1,1 -k2,2n -k3,3n | bedtools genomecov -bg -i --g /public-supool/home/longz/reference/cs/iwgsc1/cs.genome > ../bedgraph/$x.bedgraph &&
SEACR_1.3.sh ../bedgraph/$x.bedgraph 0.05 norm stringent ../bedgraph/$x
"; done > command.sh

# 2. MACS2 for narrow
for i in `ls ../bwa/*.rmdup.bam`; do x=${i/..\/bwa\//}; x=${x/.rmdup.bam/}; 
echo "macs2 callpeak -t $i -p 1e-3 -f BAMPE -g 14600000000 --keep-dup all -n $x --outdir ../macs2";
done > command.sh

# 2. MACS2 for broad
for i in `ls ../bwa/*.rmdup.bam`; do x=${i/..\/bwa\//}; x=${x/.rmdup.bam/}; 
echo "macs2 callpeak -t $i --broad --broad-cutoff 0.05 -f BAMPE -g 2700000000 --keep-dup all -n $x --outdir ../macs2";
done > command.sh

# 3. overlap peaks
for i in `ls ../bedgraph/*bed`; do x=${i/..\/bedgraph\//}; x=${x/.bed/}; 
echo "bedtools intersect -a ../macs2/$x.bed -b $i -wa > $x.overlaped.bed";
done > command.sh

