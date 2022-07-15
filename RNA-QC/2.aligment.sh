for sample in `cat sample.list`
do
hisat2 -x genome/genome  -1 clean_data/$sample.R1.fq.gz -2 clean_data/$sample.R2.fq.gz -p 8 --rna-strandness RF --fr -S  hisat2/$sample.sam  2>hisat2/$sample.hisat2.log && samtools view -uS  hisat2/$sample.sam | samtools sort -@ 5 -o  hisat2/$sample.bam && samtools index hisat2/$sample.bam
done