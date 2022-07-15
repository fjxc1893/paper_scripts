# -*- coding: utf-8 -*- 
#================================================

trim = 0
if mode == "捕获测序" or mode == "肿瘤基因组" :
    trim = config["param"]["QC"]["fastp"]["trim"]
min_lin = config["param"]["QC"]["fastp"]["length"]

##=======   local work flow   =======##
localrules: Sum_md5, fastp_plot , QC_sum, 

##=======   function   =======##



##=======   rule   =======##

##=======   check fastq   =======##
rule sample_format:
    input:
        sample=data_check
    output:
        RD1 = "Raw_data/{sample}/{sample}.raw.R1.fastq.gz",
        RD2 = "Raw_data/{sample}/{sample}.raw.R2.fastq.gz"
    run:
        raw = namedic[wildcards.sample]
        shell("{mp} && {perl} && perl {name_format} -raw {raw} -n {wildcards.sample} ")


##=======   QC fastp   =======##
rule QC_fastp:
    input:
        RD1 = "Raw_data/{sample}/{sample}.raw.R1.fastq.gz",
        RD2 = "Raw_data/{sample}/{sample}.raw.R2.fastq.gz"
    output:
        CD1 = "Clean_data/{sample}/{sample}.R1.fq.gz",
        CD2 = "Clean_data/{sample}/{sample}.R2.fq.gz",
        json = "result/01.QC/Fastp/{sample}.QC.fastp.json",
        html = "result/01.QC/Fastp/{sample}.QC.fastp.html",
        md5_1 = "result/01.QC/MD5/{sample}.R1.fq.gz.md5",
        md5_2 = "result/01.QC/MD5/{sample}.R2.fq.gz.md5"
    params:
        n_base_limit = 5,
        window_size = 4,
        mean_quality = 20,
        qualified_quality = 15
    threads: 4
    log:
        "logs/01.QC/QC_fastp/{sample}.log"
    benchmark: 
        "benchmarks/01.QC/QC_fastp/{sample}.txt" 
    shell:
        """
        {mp}
        {fastp}
        fastp -a auto --adapter_sequence_r2 auto --detect_adapter_for_pe -f {trim} -F {trim} -w {threads} -i {input.RD1} -I {input.RD2} --cut_front --cut_tail --n_base_limit {params.n_base_limit} --cut_window_size {params.window_size} --cut_mean_quality {params.mean_quality} --length_required {min_lin} --qualified_quality_phred {params.qualified_quality} -o {output.CD1} -O {output.CD2} --json {output.json} --html {output.html} >&2 2>{log}
        {md5sum} {output.CD1} > {output.md5_1}
        {md5sum} {output.CD2} > {output.md5_2}
        """

rule fastp_plot:
    input:
        json = "result/01.QC/Fastp/{sample}.QC.fastp.json"
    output:
        stat = "result/01.QC/QC_stat/{sample}/{sample}.QC.stat.xls",
        quality_png = "result/01.QC/QC_stat/{sample}/{sample}.quality.png",
        acgtn_png = "result/01.QC/QC_stat/{sample}/{sample}.acgtn.png",
        QCsummary_png = "result/01.QC/QC_stat/{sample}/{sample}.QCsummary.png"
    params:
        outdir = "result/01.QC/QC_stat/"
    shell:
        """
        {mp}
        {perl}
        perl {qc_stat} -i {input.json} -s {wildcards.sample} -o {params.outdir} 
        """

rule Sum_md5:
    input:
        md5=expand("result/01.QC/MD5/{sample}.R{t}.fq.gz.md5",sample=samples,t=[1,2])
    output:
        expand("Clean_data/{compact_num}.md5.txt",compact_num=compact_num)
    shell:
        """
        cat {input.md5} | sed -r 's#\s.*/(\S+)/#\t\\1/#g' - > {output} 
        """


##=======   QC fastqc   =======##
rule QC_fastqc:
    input:
        CD1 = rules.QC_fastp.output.CD1,
        CD2 = rules.QC_fastp.output.CD2
    output:
        html1 = "result/01.QC/Fastqc/{sample}.R1_fastqc.html",
        html2 = "result/01.QC/Fastqc/{sample}.R2_fastqc.html",
        zip1 = "result/01.QC/Fastqc/{sample}.R1_fastqc.zip",
        zip2 = "result/01.QC/Fastqc/{sample}.R2_fastqc.zip"
    params:
        outdir = "result/01.QC/Fastqc/"
    threads: 4
    log:
        "logs/01.QC/QC_fastqc/{sample}.log"
    benchmark:
        "benchmarks/01.QC/QC_fastqc/{sample}.txt"
    shell:
        """
        {mp}
        {fastqc}
        fastqc  -t {threads} {input.CD1} {input.CD2} -o {params.outdir} >&2 2>{log}
        """


rule fastqc_multiqc:
    input:
        fastqc_html = expand("result/01.QC/Fastqc/{sample}.R{t}_fastqc.html",sample=samples,t=[1,2]),
        fastqc_zip = expand("result/01.QC/Fastqc/{sample}.R{t}_fastqc.zip",sample=samples,t=[1,2])
    output:
        report_html = "result/01.QC/multiQC_Fastqc/multiqc_report.html"
    params:
        indir = "result/01.QC/Fastqc/",
        outdir = "result/01.QC/multiQC_Fastqc/"
    log:
        "logs/01.QC/fastqc_multiqc/fastqc_multiqc.log"
    benchmark:
        "benchmarks/01.QC/fastqc_multiqc/fastqc_multiqc.txt"
    shell:
        """
        set +eu
        {mp}
        {DNA}
        {multiqc} -f -o {params.outdir} {params.indir} -ip >&2 2>{log}
        set -eu
        """


##=======   QC summary   =======##
rule QC_sum:
    input:
        stat = expand("result/01.QC/QC_stat/{sample}/{sample}.QC.stat.xls",sample=samples)
    output:
        summary = "result/01.QC/QC_stat/QC.summary.xls"
    shell:
        """
        cat {input.stat} | awk '!a[$0]++' > {output.summary}
        """

