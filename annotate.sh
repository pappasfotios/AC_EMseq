#!/bin/bash

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <input_bed_file> <slop_value> <database (zf or hs for zebrafish or human protein database respectively)>"
    exit 1
fi

input_bed="$1"
slop_value="$2"
database="$3"

# Activate conda environment
source /home/fotis/miniconda3/etc/profile.d/conda.sh
conda activate biotools

output_dir="./go${database}"
mkdir -p "$output_dir"

# Set database path
if [ "$database" == "zf" ]; then
    diamond_db="/home/fotis/analysis/EMseq/enrichment/zebrafishProt2"
elif [ "$database" == "hs" ]; then
    diamond_db="/home/fotis/analysis/EMseq/enrichment/human"
else
    echo "Invalid database option. Use 'zf' for zebrafish or 'hs' for human."
    exit 1
fi

# Run bedtools slop
bedtools slop -i "$input_bed" -b "$slop_value" -g /home/fotis/analysis/Genome/Salvelinus_hybrid/ncbi_dataset/data/GCF_002910315.2/GCF_002910315.2_ASM291031v2_genomic.fna.fai > "${output_dir}/regions_slopped.bed"

# Intersect
bedtools intersect -a "${output_dir}/regions_slopped.bed" -b /home/fotis/analysis/EMseq/enrichment/gene_bodies/Salpinus_genes.gff -wb > "${output_dir}/regions_intersected.bed"

cat "${output_dir}/regions_intersected.bed" | cut -f 1,7,8 > "${output_dir}/filter_intersected.bed"

awk -F';' '/Name=/ {for (i=1; i<=NF; i++) if ($i ~ /Name=/) {split($i, name, "="); print name[2]; break}}' "${output_dir}/regions_intersected.bed" > "${output_dir}/gene_names.txt"

# Get fasta
bedtools getfasta -fi /home/fotis/analysis/Genome/Salvelinus_hybrid/ncbi_dataset/data/GCF_002910315.2/GCF_002910315.2_ASM291031v2_genomic.fna -bed "${output_dir}/filter_intersected.bed" -fo "${output_dir}/gen_regions.fa"

# BlastX
diamond blastx -d "${diamond_db}" -q "${output_dir}/gen_regions.fa" -e 10e-10 -p 8 -o "${output_dir}/${database}_hits"

sort -k1,1 -k11,11g "${output_dir}/${database}_hits" | sort --merge -u -k1,1 | awk '{print $2}' | cut -d "." -f 1 > "${output_dir}/${database}_proteins"
