setwd("C:\\Users\\fopa0001\\Downloads\\OneDrive_1_04-12-2023")

library(Biobase)
library(comethyl)
library(GenomicRanges)
library(tidyverse)

# Read file
bs <- readRDS("Filtered_BSseq.rds")

# Set binary traits
bs@colData$statusDen <- rep("high", nrow(bs@colData))
bs@colData$statusDen[order(bs@colData$Density)[21:27]] <- NA
bs@colData$statusDen[order(bs@colData$Density)[1:20]] <- "low"

bs@colData$statusVCL <- rep("high", nrow(bs@colData))
bs@colData$statusVCL[order(bs@colData$VCL)[21:27]] <- NA
bs@colData$statusVCL[order(bs@colData$VCL)[1:20]] <- "low"


#r <- promoters[,1:3]
#ac_annot <- makeGRangesFromDataFrame(r, 
#                                     keep.extra.columns = F,
#                                     ignore.strand = TRUE)

# Get and filter CpG regions
regions <- getRegions(bs, n=5) #, custom=ac_annot)

regions<- filterRegions(regions, covMin = 10, methSD = 0.1, file = "Filtered_Regions.txt")

# Get methylation per region
meth <- getRegionMeth(regions, 
                      bs = bs, 
                      file = "Region_Methylation.rds")

mod <- model.matrix(~1, data = pData(bs))

methAdj <- adjustRegionMeth(meth, 
                            PCs = mod,
                            file = "Adjusted_Region_Methylation.rds")

#r$seqnames <- as.character(r$seqnames)
#c <- bsseq::getCoverage(bs, regions = promoters[,1:3], type = "Cov", what = "perRegionAverage")
#r <- r[complete.cases(m) & matrixStats::rowSds(m) > 0.1 & matrixStats::rowMins(c) >= 10 ,]
#r <- as.data.frame(r)
#m <- m[complete.cases(m) & matrixStats::rowSds(m) > 0.1 & matrixStats::rowMins(c) >= 10,]

# Association tests for sperm density
den_assoc <- data.frame(pValue=as.numeric(), diff=as.numeric())

for (i in seq(1,ncol(methAdj))) {
  
  ttest <- t.test(methAdj[,i]~bs$statusDen)
  den_assoc <- rbind(den_assoc, data.frame(pValue=ttest$p.value, diff=ttest$estimate[1] - ttest$estimate[2]))
  
}

den_assoc <- cbind(den_assoc, regions)
den_assoc$adjP <- p.adjust(den_assoc$pValue, method = "BH")

# Association tests for sperm velocity
vcl_assoc <- data.frame(pValue=as.numeric(), diff=as.numeric())

for (i in seq(1,ncol(methAdj))) {
  
  ttest <- t.test(methAdj[,i]~bs$statusVCL)
  vcl_assoc <- rbind(vcl_assoc, data.frame(pValue=ttest$p.value, diff=ttest$estimate[1] - ttest$estimate[2]))
  
}

vcl_assoc <- cbind(vcl_assoc, regions)
vcl_assoc$adjP <- p.adjust(vcl_assoc$pValue, method = "BH")


# Plots
vlc2 <- ggplot(vcl_assoc, aes(x=diff,y=-log10(pValue))) + 
  geom_point(aes(col=adjP > 0.05)) + 
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05/nrow(regions)), col="limegreen", linetype = "dashed") + 
  scale_color_manual(values = c("red", "grey"))  + 
  ggtitle(label = "high vs low Velocity") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
  ylim(-0.2, 8)
vlc2


vlc1 <- ggplot(den_assoc, aes(x=diff,y=-log10(pValue))) + 
  geom_point(aes(col=adjP > 0.05)) +
  theme_minimal() + 
  geom_hline(yintercept = -log10(0.05/nrow(regions)), col="limegreen", linetype = "dashed") + 
  scale_color_manual(values = c("red", "grey"))  + 
  ggtitle(label = "high vs low Density") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
  ylim(-0.2, 8)
vlc1


cowplot::plot_grid(vlc1,NULL,vlc2, ncol=3, align = "vh", rel_widths = c(1,0.1,1))


den_assoc$chrom <- den_assoc$chr
den_assoc$chrom[grep("NW", den_assoc$chr)] <- "UK"
den_assoc$chrom <- match(den_assoc$chrom, unique(den_assoc$chrom))

qqman::manhattan(den_assoc, chr = "chrom", bp="start", p="pValue", snp = "RegionID", 
                 chrlabs = c(unique(den_assoc$chr)[1:39],"UA"), logp = T, genomewideline = -log10(0.05/nrow(den_assoc)), 
                 suggestiveline = -log10(mean(max(den_assoc$pValue[den_assoc$adjP<0.05]),min(den_assoc$pValue[den_assoc$adjP>0.05]))))
qqman::qq(den_assoc$pValue)



vcl_assoc$chrom <- vcl_assoc$chr
vcl_assoc$chrom[grep("NW", vcl_assoc$chr)] <- "UK"
vcl_assoc$chrom <- match(vcl_assoc$chrom, unique(vcl_assoc$chrom))

qqman::manhattan(vcl_assoc, chr = "chrom", bp="start", p="pValue", snp = "RegionID", 
                 chrlabs = c(unique(vcl_assoc$chr)[1:39],"UA"), logp = T, genomewideline = -log10(0.05/nrow(vcl_assoc)), 
                 suggestiveline = -log10(mean(max(vcl_assoc$pValue[vcl_assoc$adjP<0.05]),min(vcl_assoc$pValue[vcl_assoc$adjP>0.05]))))
qqman::qq(vcl_assoc$pValue)





calculate_EST <- function(methylationMatrix, groupingVariable) {
  # Remove samples with NA in the grouping variable
  validIndices <- !is.na(groupingVariable)
  
  # Split the methylation data into two groups based on the grouping variable
  group1 <- groupingVariable == "high"
  group2 <- groupingVariable == "low"
  
  # across variance
  numerator <- matrixStats::rowSds(methylationMatrix[,validIndices])
  
  # within variance
  overallMean <- (matrixStats::rowSds(methylationMatrix[,group1])^2 + matrixStats::rowSds(methylationMatrix[,group2])^2) / 2
  
  # Calculate E_ST for each CpG site
  E_ST <- numerator^2 / (numerator^2 + overallMean)
  
  # Ensure the output is a numeric vector
  return(as.numeric(E_ST))
}


calculate_EST <- function(methylationMatrix, groupingVariable) {
  # Ensure 'matrixStats' library is available
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop("The 'matrixStats' package is required but not installed.")
  }
  
  # Filter out samples with NA in groupingVariable
  validSamples <- !is.na(groupingVariable)
  methylationMatrix <- methylationMatrix[, validSamples]
  groupingVariable <- groupingVariable[validSamples]
  
  # Convert grouping variable to factor if it's not already
  groupingVariable <- factor(groupingVariable)
  
  # Check for exactly two levels in groupingVariable
  if(length(levels(groupingVariable)) != 2) {
    stop("Grouping variable must have exactly two levels.")
  }
  
  # Calculate group means
  groupMeans <- tapply(1:ncol(methylationMatrix), groupingVariable, function(i) rowMeans(methylationMatrix[, i], na.rm = TRUE))
  
  # Calculate between-group variance (squared difference of means)
  betweenGroupVariance <- (groupMeans[[1]] - groupMeans[[2]])^2
  
  # Calculate within-group variance for each group and then average them
  withinGroupVariance <- sapply(groupMeans, function(mean, group) {
    diffSq <- (methylationMatrix[, groupingVariable == group] - mean)^2
    rowMeans(diffSq, na.rm = TRUE)
  }, group=levels(groupingVariable))
  
  averageWithinGroupVariance <- rowMeans(withinGroupVariance, na.rm = TRUE)
  
  # Calculate E_ST for each CpG site
  E_ST <- betweenGroupVariance / (betweenGroupVariance + averageWithinGroupVariance)
  
  # Return the E_ST values
  return(E_ST)
}

# Usage example
# Assuming 'methAdj' is your matrix with samples as columns and CpG sites as rows,
# and 'bs$statusDen' is the binary factor indicating sperm density status for each sample:
E_ST_values <- calculate_EST(t(methAdj), bs$statusDen)
head(E_ST_values)



E_ST_values <- calculate_EST(t(methAdj), bs$statusDen)
