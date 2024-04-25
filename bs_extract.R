library(Biobase)
library(comethyl)
library(GenomicRanges)
library(openxlsx)
library(dplyr)
library(readr)

setwd("C:\\Users\\fopa0001\\Downloads\\OneDrive_1_04-12-2023")

options(stringsAsFactors = FALSE)
Sys.setenv(R_THREADS = 1)
WGCNA::enableWGCNAThreads(nThreads = 4)


colData <- read.xlsx("pheno_comethyl.xlsx", rowNames = TRUE)
colData$Density[colData$Density>10000] <- 10000
colData <- colData %>% rename(Concentration = Density)

chr = read.table(file="chrom_list.txt",
                 header = F,
                 stringsAsFactors = F)

chr_list <- chr$V1

bs <- getCpGs(colData,
              path = "./",
              chroms = chr_list,
              file = "Unfiltered_BSseq.rds")

bs <- filterCpGs(bs, cov = 10, perSample = 0.85, file = "Filtered_BSseq.rds")

# Use transitions information from vcf to filter out SNPs
transitions <- read.table("transitions052al.bed", header=F)
names(transitions) <- c("chr", "start", "end")
gr_transitions <- makeGRangesFromDataFrame(transitions)

overlaps <- findOverlaps(gr_transitions, granges(bs), ignore.strand=T)

bs_indices <- subjectHits(overlaps)

bs <- bs[-bs_indices]

write_rds(bs, file = "noSNPs_BSseq.rds")
