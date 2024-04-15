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
  ggtitle(label = "High vs Low Curvilinear Velocity (VCL)") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
  labs(x="Difference of means", y=expression(-log[10](p))) +
  ylim(-0.2, 8)
vlc2


vlc1 <- ggplot(den_assoc, aes(x=diff,y=-log10(pValue))) + 
  geom_point(aes(col=adjP > 0.05)) +
  theme_minimal() + 
  geom_hline(yintercept = -log10(0.05/nrow(regions)), col="limegreen", linetype = "dashed") + 
  scale_color_manual(values = c("red", "grey"))  + 
  ggtitle(label = "High vs Low Concentration") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
  labs(x="Difference of means", y=expression(-log[10](p))) +
  ylim(-0.2, 8)
vlc1


cowplot::plot_grid(vlc1,NULL,vlc2, ncol=3, align = "vh", rel_widths = c(1,0.1,1))


den_assoc$chrom <- den_assoc$chr
den_assoc$chrom[grep("NW", den_assoc$chr)] <- "UK"
den_assoc$chrom <- match(den_assoc$chrom, unique(den_assoc$chrom))

qqman::manhattan(den_assoc, chr = "chrom", bp="start", p="pValue", snp = "RegionID", col = c("skyblue3","grey40"),
                 chrlabs = c(paste0("LG",c(1,2,3,"4p","4q.1:29","4q.2",5,"6.1","6.2",seq(7,28),seq(30,37))),"UA"), logp = T, genomewideline = -log10(0.05/nrow(den_assoc)), 
                 suggestiveline = 0.5*(-log10(max(den_assoc$pValue[den_assoc$adjP<0.05])) + -log10(min(den_assoc$pValue[den_assoc$adjP>0.05]))), 
                 las=2, xlab=" ", main="Concentration")
qqman::qq(den_assoc$pValue)



vcl_assoc$chrom <- vcl_assoc$chr
vcl_assoc$chrom[grep("NW", vcl_assoc$chr)] <- "UK"
vcl_assoc$chrom <- match(vcl_assoc$chrom, unique(vcl_assoc$chrom))

qqman::manhattan(vcl_assoc, chr = "chrom", bp="start", p="pValue", snp = "RegionID", col = c("skyblue3","grey40"),
                 chrlabs = c(paste0("LG",c(1,2,3,"4p","4q.1:29","4q.2",5,"6.1","6.2",seq(7,28),seq(30,37))),"UA"), logp = T, genomewideline = -log10(0.05/nrow(vcl_assoc)), 
                 suggestiveline = 0.5*(-log10(max(vcl_assoc$pValue[vcl_assoc$adjP<0.05])) + -log10(min(vcl_assoc$pValue[vcl_assoc$adjP>0.05]))), 
                 las=2, xlab=" ", main="Curvilinear Velocity (VCL)")
qqman::qq(vcl_assoc$pValue)



sig_regions <- cbind(methAdj[,colnames(methAdj) %in% den_assoc$RegionID[den_assoc$adjP < 0.05]], methAdj[,colnames(methAdj) %in% vcl_assoc$RegionID[vcl_assoc$adjP < 0.05]])
write.table(sig_regions,file = "sigAssocRegions.txt")
