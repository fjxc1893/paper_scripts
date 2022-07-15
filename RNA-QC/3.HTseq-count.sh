for sample in `cat sample.list`
do
htseq-count -i gene_id  -f bam -s reverse -r name  hisat2/$sample.bam genome/gene.gtf > Quantification/$sample.counts.t
done