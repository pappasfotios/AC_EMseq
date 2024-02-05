setwd("C:\\Users\\fopa0001\\Downloads\\OneDrive_1_04-12-2023")

library(ggplot2)
library(data.table)
library(GenomicRanges)
library(dplyr)

# Import data
bs <- readRDS("Filtered_BSseq.rds")

CpGislands <- read.table("cpgi.bed")
colnames(CpGislands) <- c("chr","start","end")

CpGshores <- read.table("cpgi_shores.bed")
colnames(CpGshores) <- c("chr","start","end")

gene_bodies <- read.table("Salpinus_genes.gff", stringsAsFactors = F, sep="\t",)
gene_bodies <- gene_bodies[,c(1,4,5,7)]
colnames(gene_bodies) <- c("chr","start","end","strand")
gene_bodies$uID <- paste0(gene_bodies$chr,":",gene_bodies$start)
gene_bodies$uID[gene_bodies$strand=="-"] <- paste0(gene_bodies$chr[gene_bodies$strand=="-"],":",gene_bodies$end[gene_bodies$strand=="-"])
gene_bodies <- gene_bodies[!duplicated(gene_bodies),]

promoters <- read.table("Salpinus_promoters_1kb.gff", stringsAsFactors = F, sep="\t",)
promoters$GeneID <- sub(".*GeneID:(\\d+).*", "\\1", promoters$V9)
promoters <- promoters[,c(1,4,5,7,10)]
colnames(promoters) <- c("chr","start","end","strand","GeneID")
promoters <- promoters[!duplicated(promoters),]

# Exons
gene_features <- read.table("genomic.gff", stringsAsFactors = F, sep="\t",)
gene_features$GeneID <- sub(".*GeneID:(\\d+).*", "\\1", gene_features$V9)
gene_features <- gene_features[gene_features[,3]=="exon",c(1,4,5,7,10)]
colnames(gene_features) <- c("chr","start","end","strand","GeneID")
gene_features <- gene_features[!duplicated(gene_features),]
gene_features <- gene_features[!gene_features$chr=="NC_000861.1",]
gene_features$GeneID <- as.integer(gene_features$GeneID)

first_exons <- setDT(gene_features)[, .(chr = chr, 
                                  start = ifelse(strand == "+", min(start), max(start)),
                                  end = ifelse(strand == "+", min(end), max(end))),
                              by = GeneID][, .(chr = first(chr), start = first(start), end = first(end)), by = GeneID]

# Introns
gene_features <- gene_features[order(gene_features$GeneID,gene_features$start, gene_features$end),]
geneIDs_with_multiple_entries_pos <- setDT(gene_features)[strand == "+", .(N = .N), by = GeneID][N > 1, GeneID]
geneIDs_with_multiple_entries_neg <- setDT(gene_features)[strand == "-", .(N = .N), by = GeneID][N > 1, GeneID]

pos_introns <- setDT(gene_features)[GeneID %in% geneIDs_with_multiple_entries_pos,
                                    .(chr = chr[1],
                           start = ifelse(end < dplyr::lead(start), end + 1, ifelse(end < dplyr::lead(start,n=2L), dplyr::lead(end +1), dplyr::lead(end +1, n=2L))),
                           end = ifelse(end < dplyr::lead(start),dplyr::lead(start - 1), ifelse(end < dplyr::lead(start,n=2L), dplyr::lead(end +1, n=2L), dplyr::lead(end +1, n=2L))),
                           strand = strand[1]),
                         by = GeneID]
pos_introns <- pos_introns[complete.cases(pos_introns),]

neg_introns <- setDT(gene_features)[GeneID %in% geneIDs_with_multiple_entries_neg, 
                                    .(chr = chr[1],
                                      start = ifelse(end < dplyr::lead(start), end + 1, ifelse(end < dplyr::lead(start,n=2L), dplyr::lead(end +1), dplyr::lead(end +1, n=2L))),
                                      end = ifelse(end < dplyr::lead(start),dplyr::lead(start - 1), ifelse(end < dplyr::lead(start,n=2L), dplyr::lead(end +1, n=2L), dplyr::lead(end +1, n=2L))),
                                      strand = strand[1]),
                                    by = GeneID]
neg_introns <- neg_introns[complete.cases(neg_introns),]

introns <- rbind(pos_introns,neg_introns)
introns <- introns[order(chr, start)]
introns$width <- introns$end - introns$start + 1

### Check validity
introns <- introns[introns$width>=20,]

## First introns
first_introns <- setDT(introns)[, .(chr = chr, 
                                    start = ifelse(strand == "+", min(start), max(start)),
                             end = ifelse(strand == "+", min(end), max(end))),
                         by = GeneID][, .(chr = first(chr), start = first(start), end = first(end)), by = GeneID]

# Intergenic regions
gene_bodies_gr <- makeGRangesFromDataFrame(gene_bodies,
                                           keep.extra.columns = F,
                                           seqnames.field = "chr")
## reduce?
intergenic <- as.data.frame(gaps(reduce(gene_bodies_gr)))

# Function to extract region means and SDs
get_meth_stats <- function(ant){
  
  m <- bsseq::getMeth(bs, ant ,type = 'raw', what = 'perRegion')
  
  means <- matrixStats::rowMeans2(m)
  sds <- matrixStats::rowSds(m)
  
  return(data.frame(means=means, sds=sds))
}


# Extract region means and SDs
islands_ms <- get_meth_stats(CpGislands)

shores_ms <- get_meth_stats(CpGshores)

genes_ms <- get_meth_stats(gene_bodies[,1:3])

promoters_ms <- get_meth_stats(promoters)
promoters_ms$GeneID <- as.character(promoters$GeneID)

exons_ms <- get_meth_stats(setDF(gene_features))

introns_ms <- get_meth_stats(setDF(introns[,2:4]))

intergenic_ms <- get_meth_stats(intergenic)

first_introns_ms <- get_meth_stats(setDF(first_introns[,2:4]))
first_introns_ms$GeneID <- as.character(first_introns$GeneID)

combo <- dplyr::inner_join(first_introns_ms, promoters_ms, by="GeneID")  ## further explore promoter - first intron - nearest CpG island relationship

# Combine in long format
data <- rbind(data.frame(feature="CpG_islands", means=islands_ms$means, sds=islands_ms$sds), 
              data.frame(feature="CpG_shores", means=shores_ms$means, sds=shores_ms$sds), 
              data.frame(feature="Gene_bodies", means=genes_ms$means, sds=genes_ms$sds),
              data.frame(feature="Proximal_upstream_regions", means=promoters_ms$means, sds=promoters_ms$sds),
              data.frame(feature="Exons", means=exons_ms$means, sds=exons_ms$sds),
              data.frame(feature="Intergenic", means=intergenic_ms$means, sds=intergenic_ms$sds),
              data.frame(feature="Introns", means=introns_ms$means, sds=introns_ms$sds),
              data.frame(feature="First_introns", means=first_introns_ms$means, sds=first_introns_ms$sds)
              )

# Violin plots
means_violin <- ggplot(data, aes(x=feature, y=means)) + 
  geom_violin(aes(color=feature), linewidth=1) +
  #geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.3, binwidth = 0.0004, aes(fill=feature), alpha=0.3) + 
  ylab("Averaged mean methylation per region") + 
  xlab("Genomic feature") + 
  theme_bw() + 
  theme(legend.position = "none") 

sds_violin <- ggplot(data, aes(x=feature, y=sds)) + 
  geom_violin(aes(color=feature), linewidth=1) + 
  ylab("Sd of mean methylation per region") +
  xlab("Genomic feature") + 
  theme_bw() + 
  theme(legend.position = "none") 

cowplot::plot_grid(means_violin,sds_violin, ncol = 1, align = "vh", labels = LETTERS[1:2])

# Over- and hypo- methylated features
under_promoters <- cbind(promoters[promoters_ms$means<0.05 & promoters_ms$sds<0.05,1:3], promoters_ms[promoters_ms$means<0.05 & promoters_ms$sds<0.05,])
under_promoters <- under_promoters[complete.cases(under_promoters),]
write.table(under_promoters[,1:3], file = "hypomethylated_promoters.bed", sep = "\t", quote = F, row.names = F, col.names = F)

over_promoters <- cbind(promoters[promoters_ms$means>0.95 & promoters_ms$sds<0.05,1:3], promoters_ms[promoters_ms$means>0.95 & promoters_ms$sds<0.05,])
over_promoters <- over_promoters[complete.cases(over_promoters),]
write.table(over_promoters[,1:3], file = "hypermethylated_promoters.bed", sep = "\t", quote = F, row.names = F, col.names = F)
