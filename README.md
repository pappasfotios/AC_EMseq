## Analysis of EM-seq data
Run filtering, alignment and methylation extraction pipeline (uses Bismark): **Snakefile**

```{shell}
snakemake --cores 4 --keep-going
```

Detect CpG islands and their respective shores in the reference genome: **CpGislandsAndShores.py**

Explorative data analysis of CpG methylation landscape : **methExplore.R**
![CpG_summary_new](https://github.com/pappasfotios/AC_EMseq/assets/49454378/9b7f29aa-39ce-427d-9bbf-a9d55b59616d)

Co-methylation of CpGs by physical distance : **comethylation_decay.R**

Comparison of genetic and epigenetic similarity matrices : **methVSgen.R**

Example of co-methylation module analysis considering promoters (1kb upstream regions): **comethyl_promoters1kb.Rmd**

Case-control association tests for sperm density and velocity: **Assoc.R**

Script to annotate signals : **annotate.sh** (have previously created Diamond DB with human (hs) and zebrafish (zf) proteome)



```{shell}
bash annotate.sh <input_bed_file> <slop_value in bp> <database (zf or hs for zebrafish or human protein database respectively)>
```
