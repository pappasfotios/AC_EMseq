SAMPLES = glob_wildcards("raw_reads/{sample}_R1_001.fastq.gz").sample
METH_OUTPUT = "meth_extracted"
print("The sample names to be processed are:", SAMPLES)
#PAIRED = ["R1","R2"]

rule all:
    input: 
        expand("meth_extracted/{sample}_bismark_bt2_pe.deduplicated.bismark.cov.gz", sample=SAMPLES)


rule fastp:
    output:
        read1 = "trimmed_{sample}_R1_001.fastq.gz",
        read2 = "trimmed_{sample}_R2_001.fastq.gz"
    input:
        read1 = "raw_reads/{sample}_R1_001.fastq.gz",
        read2 = "raw_reads/{sample}_R2_001.fastq.gz"
    log:
        "logs_fastp/{sample}.log"    
    shell:
        "fastp -i {input.read1} -I {input.read2} -o {output.read1} -O {output.read2} --html {wildcards.sample}.html 2> {log}"  
        
rule bismark_align:
    output:
        bam = temporary("{sample}_bismark_bt2_pe.bam"),
        report = "{sample}_bismark_bt2_PE_report.txt"
    input:
        read1 = "trimmed_{sample}_R1_001.fastq.gz",
        read2 = "trimmed_{sample}_R2_001.fastq.gz"
    log:
        "bismark_align/{sample}.log"
    threads:
        4    
    shell:
        r"""
        bismark --parallel {threads} --genome /home/csas0002/Ref_genomes/Salvelinus -1 {input.read1} -2 {input.read2} 2> {log}
        mv trimmed_{wildcards.sample}_R1_001_bismark_bt2_pe.bam {wildcards.sample}_bismark_bt2_pe.bam
        mv trimmed_{wildcards.sample}_R1_001_bismark_bt2_PE_report.txt {wildcards.sample}_bismark_bt2_PE_report.txt
        """

rule bismark_deduplicate:
    output:
        bam_dedup = "{sample}_bismark_bt2_pe.deduplicated.bam",
        report = "{sample}_bismark_bt2_pe.deduplication_report.txt"
    input:
        bam = "{sample}_bismark_bt2_pe.bam"
    log:
        "bismark_deduplicate/{sample}.log"    
    shell:
        "deduplicate_bismark --paired --bam {input.bam} 2> {log}"

rule meth_extractor:
    output:
        meth = "meth_extracted/{sample}_bismark_bt2_pe.deduplicated.bismark.cov.gz" 
        #cytosine_report = "meth_extracted/{sample}_bismark_bt2_pe.deduplicated.bismark.cov.gz"
    input:
        bam= "{sample}_bismark_bt2_pe.deduplicated.bam"
    log:
        "meth_extractor/{sample}.log"
    shell:
        r"""
        bismark_methylation_extractor --buffer_size 20G --scaffolds --paired-end --parallel 4 --gzip --comprehensive --cytosine_report --bedgraph {input.bam} --output {METH_OUTPUT} 2> {log}
        """
        
