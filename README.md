## Analysis of EM-seq data derived from sperm DNA
Run filtering, alignment and methylation extraction pipeline (uses Bismark): **[Raw reads to methylation Snakemake pipeline](https://github.com/pappasfotios/AC_EMseq/blob/main/Reads_to_meth_pipeline/Snakefile)**

Double-mask alignments and call genetic variants: **[Alignment double masking and SNP calling Snakemake pipeline](https://github.com/pappasfotios/AC_EMseq/blob/main/SNP_calling/Snakefile)**

Detect CpG islands and their respective shores in the reference genome: **[CpG island and shore detection](https://github.com/pappasfotios/AC_EMseq/blob/main/Utilities/CpGislandsAndShores.py)**

Explorative data analysis of CpG methylation landscape : **[Exploratory analysis](https://github.com/pappasfotios/AC_EMseq/tree/main/Exploratory_analysis)**
![CpG_summary_new](https://github.com/pappasfotios/AC_EMseq/assets/49454378/9b7f29aa-39ce-427d-9bbf-a9d55b59616d)

Comparison of genetic relationships and epigenetic similarities : **[Methylation by genetic background](https://github.com/pappasfotios/AC_EMseq/tree/main/MethylationByGenotype)**

Comethylation networks analysis considering promoters, CpG islands and first introns: **[Comethylation analysis](https://github.com/pappasfotios/AC_EMseq/tree/main/ComethylationNetworks_associations)**

Case-control association tests for sperm density and velocity considering clusters of adjacent CpGs: **[Annotation-free association](https://github.com/pappasfotios/AC_EMseq/tree/main/AnnotationFree_association)**

Bash script orchestrating annotation of signals : **[annotate.sh](https://github.com/pappasfotios/AC_EMseq/blob/main/Utilities/annotate.sh)** (must previously creat Diamond DB with human (hs) and zebrafish (zf) proteome)



```{shell}
bash annotate.sh <input_bed_file> <slop_value in bp> <database (zf or hs for zebrafish or human protein database respectively)>
```

Conda environments used for different analyses: **[Environments](https://github.com/pappasfotios/AC_EMseq/tree/main/Environments)**
