import pandas as pd
import subprocess

# Activate the biotools environment
activate_biotools = "source /home/fotis/miniconda3/etc/profile.d/conda.sh && conda activate biotools"
subprocess.run(activate_biotools, shell=True, executable="/bin/bash")
bedtools_path = "/home/fotis/miniconda3/envs/biotools/bin/bedtools"
newcpgreport_path = "/home/fotis/miniconda3/envs/biotools/bin/newcpgreport"

# CpG islands detection with strict options
command1 = [newcpgreport_path, "-sequence", "/home/fotis/analysis/Genome/Salvelinus_hybrid/ncbi_dataset/data/GCF_002910315.2/GCF_002910315.2_ASM291031v2_genomic.fna", "-window", "100", "-shift", "1", "-minlen", "200", "-minoe", "0.65", "-minpc", "55.", "-outfile", "sa_cpgi.report"]
subprocess.run(command1)

file_path = "./sa_cpgi.report"
data = []
current_chromosome = None

# Iterate over report lines
with open(file_path, 'r') as file:
    for line in file:
        if line.startswith("ID"):
            current_chromosome = line.split()[1]
        elif line.startswith("FT") and "CpG island" in line:
            entry = line.split()
            feature_info_index = entry.index("CpG") + 2  # Index where CpG island info starts
            start, end = map(int, entry[feature_info_index].split(".."))
            data.append([current_chromosome, start, end])

columns = ["Chromosome", "Start", "End"]
newcpg_data = pd.DataFrame(data, columns=columns)

newcpg_data.to_csv(path_or_buf="./cpgi2.bed", sep="\t", header=False, index=False)

# CpG shores
command1 = [bedtools_path, "flank", 
            "-i", "cpgi2.bed", 
            "-g", "/home/fotis/analysis/Genome/Salvelinus_hybrid/ncbi_dataset/data/GCF_002910315.2/GCF_002910315.2_ASM291031v2_genomic.fna.fai", 
            "-b", "2000"]

output_file1 = "cpgi_shores2.bed"
with open(output_file1, "w") as output:
    # Run the bedtools command and capture the output
    subprocess.run(command1, stdout=output)

# Islands and their shores
command2 = [bedtools_path, "slop", 
            "-i", "cpgi2.bed", 
            "-g", "/home/fotis/analysis/Genome/Salvelinus_hybrid/ncbi_dataset/data/GCF_002910315.2/GCF_002910315.2_ASM291031v2_genomic.fna.fai", 
            "-b", "2000"]

output_file2 = "IslandsPlusShores2.bed"
with open(output_file2, "w") as output:
    # Run the bedtools command and capture the output
    subprocess.run(command2, stdout=output)
