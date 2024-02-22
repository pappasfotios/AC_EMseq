## Analysis of EM-seq data

Detect CpG islands and their respective shores in the reference genome: CpGislandsAndShores.py

Explorative data analysis of CpG methylation landscape : methExplore.R

Co-methylation of CpGs by physical distance : comethylation_decay.R

Comparison of genetic and epigenetic similarity matrices : methVSgen.R

Example of co-methylation module analysis considering promoters (!kb upstream regions): comethyl_promoters1kb.Rmd

Script to annotate signals : annotate.sh (have previously created Diamond DB with human and zebrafish proteome)



```{shell}
bash annotate.sh <input_bed_file> <slop_value> <database (zf or hs for zebrafish or human protein database respectively)>
```
