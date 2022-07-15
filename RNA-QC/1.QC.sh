for sample in `cat sample.list`
do
java -jar  -Xmx40g  /software/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar PE  -threads 10  raw/$sample/$sample.R1.fq.gz raw/$sample/$sample.R2.fq.gz clean_data/$sample.R1.fq.gz clean_data/$sample.R1_un.fq.gz clean_data/$sample.R2.fq.gz  clean_data/$sample.R2_un.fq.gz  CROP:150 ILLUMINACLIP:/software/Trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10:8:true LEADING:3 TRAILING:3  SLIDINGWINDOW:4:15 MINLEN:50 
done