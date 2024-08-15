setwd("C:\\Users\\fopa0001\\Downloads\\OneDrive_1_04-12-2023")

library(data.table)
library(GenomicRanges)
library(dplyr)
library(ggplot2)

source("C:\\Users\\fopa0001\\Documents\\AC_EMseq\\comethylation_decay_function.R")

# Import data
bs <- readRDS("noSNPs_BSseq.rds")
genome <- read.table("//wsl.localhost/Ubuntu/home/fotis/analysis/Genome/Salvelinus_hybrid/ncbi_dataset/data/GCF_002910315.2/GCF_002910315.2_ASM291031v2_genomic.fna.fai",
                     header = F)
genome$start <- rep(1,nrow(genome))
names(genome)[1:2] <- c("chr", "end")
genome <- genome[grep("NC", genome$chr),]
genome <- genome[genome$chr != "NC_000861.1",]

genome <- makeGRangesFromDataFrame(genome)

m <- bsseq::getMeth(bs, type = "raw", what = "perBase")
r <- as.data.frame(granges(bs))
r$seqnames <- as.character(r$seqnames)

cov <- bsseq::getCoverage(bs, type = "Cov")

SDmask <- matrixStats::rowSds(m) > 0.05
COVmask <- matrixStats::rowMins(cov) >= 5

r <- r[complete.cases(m) & SDmask & COVmask,]
r <- as.data.frame(r)
m <- m[complete.cases(m) & SDmask & COVmask,]

short_range <- comethylation_decay(bin_widths = seq(100, 1000, by=100), nregions = 2000, length_sd = 0, genome = genome,
                                 m = m, r = r, seed = 2118)


long_range <- comethylation_decay(bin_widths = seq(10000, 100000, by=10000), nregions = 1000, length_sd = 0, genome = genome,
                                   m = m, r = r, seed = 1837)


reduced_short <- setDT(short_range)[, .(value= mean(abs_cor)),by=width]

png(filename = "AbsCometh_short.png", width = 1800, height = 1800, res=300)
ggplot(data=reduced_short, aes(x=width, y=value)) + theme_light() + 
  geom_line(col="darkred", lwd=1) +
  geom_smooth(col="grey20", lwd=0.7) + 
  ylab("Absolute correlation") + 
  xlab("Bin width (Kb)")
dev.off()

png(filename = "ComethBars_short.png", width = 1800, height = 1800, res=300)
ggplot(data=short_range, aes(x=width, y=cor, group=width)) + #ylim(c(-0.005,0.1)) +
  geom_boxplot(fill="#E30B5D", outlier.shape = NA) + xlab("Bin width (Kb)") + ylab("Correlation") +
  theme_minimal()
dev.off()

png(filename = "Abs_ComethBars_short.png", width = 1800, height = 1800, res=300)
ggplot(data=short_range, aes(x=width, y=abs_cor, group=width)) + #ylim(c(0.11,0.17)) +
  geom_boxplot(fill="#E30B5D", outlier.shape = NA) + xlab("Bin width (Kb)") + ylab("Absolute correlation") +
  theme_bw()
dev.off()

reduced_long <- setDT(long_range)[, .(value= mean(abs_cor)),by=width]

png(filename = "AbsCometh_long.png", width = 1800, height = 1800, res=300)
ggplot(data=reduced_long, aes(x=width, y=value)) + theme_light() + 
  geom_line(col="darkred", lwd=1) +
  geom_smooth(col="grey20", lwd=0.7) + 
  ylab("Absolute correlation") + 
  xlab("Bin width (Kb)")
dev.off()

png(filename = "ComethBars_long.png", width = 1800, height = 1800, res=300)
ggplot(data=long_range, aes(x=width, y=cor, group=width)) + ylim(c(-0.01,0.1)) +
  geom_boxplot(fill="#E30B5D", outlier.shape = NA) + xlab("Bin width (Kb)") + ylab("Correlation") +
  theme_minimal()
dev.off()

png(filename = "Abs_ComethBars_long.png", width = 1800, height = 1800, res=300)
ggplot(data=long_range, aes(x=width, y=abs_cor, group=width)) + ylim(c(0.11,0.17)) +
  geom_boxplot(fill="#E30B5D", outlier.shape = NA) + xlab("Bin width (Kb)") + ylab("Absolute correlation") +
  theme_bw()
dev.off()
