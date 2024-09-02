setwd("C:\\Users\\fopa0001\\Downloads\\OneDrive_1_04-12-2023")

library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(RColorBrewer)
library(pedigree)
library(ggplot2)
library(GGally)

bs <- readRDS("noSNPs_BSseq.rds")

# ERM
m <- bsseq::getMeth(bs, type = "raw", what = "perBase")
cov <- bsseq::getCoverage(bs, type = "Cov")

SDmask <- matrixStats::rowSds(m) > 0.1
COVmask <- matrixStats::rowMins(cov) >= 10

m <- m[complete.cases(m) & SDmask & COVmask,]

MRM <- scale(t(m)) %*% t(scale(t(m)))
MRM <- as.matrix(MRM/nrow(m))

colnames(MRM) <- seq(1,47)
rownames(MRM) <- seq(1,47)

corrplot::corrplot(MRM, method = "square", diag = T, is.corr = F)

ped <- read.table("//wsl.localhost/Ubuntu/home/fotis/analysis/BLUP/ac_ped_updated_2023.txt", header = T, stringsAsFactors = T)
ped <- ped[!duplicated(ped$Id),]

tagIDs <- openxlsx::read.xlsx("C:/Users/fopa0001/Downloads/Arctic_charr_EMseq_pheno.xlsx")
tagIDs$Id_tag <- as.factor(tagIDs$Id_tag)
tagIDs <- tagIDs[-grep(398, tagIDs$Id_seq),]

# A matrix
sub <- ped$Id %in% tagIDs$Id_tag
A <- makeA(ped = ped[,1:3], which = sub)

A <- read.table(file = "./A.txt")
A_IDs <- unique(A$V1)                                    
A$Id1 <- ped$Id[A$V1]
A$Id2 <- ped$Id[A$V2]

rA <- reshape(A[,3:5], idvar = "Id1", timevar = "Id2", direction = "wide")
row.names(rA) <- rA$Id1
rA <- rA[,2:48]
names(rA) <- row.names(rA)

bs@colData$Id_seq <- rep(NA, nrow(bs@colData))
subset_rowsA <- grep("UC-", row.names(bs@colData))
subset_rowsA <- subset_rowsA[!is.na(subset_rowsA)]

bs@colData$Id_seq[subset_rowsA] <- substr(rownames(bs@colData)[subset_rowsA], 9, 18)
bs@colData$Id_seq[-subset_rowsA] <- substr(rownames(bs@colData)[-subset_rowsA], 9, 21)
colData <- dplyr::inner_join(as.data.frame(bs@colData), tagIDs[,1:2], by="Id_seq")
colData$key <- seq(1,47)

row.names(rA) <- colData$key[colData$Id_tag == names(rA)]
names(rA) <- row.names(rA)

A <- reshape::melt(rA,)

corrplot::corrplot(as.matrix(rA), method = "square", diag = F, is.corr = F, type = "lower")

# Genomic relationships
GRM1 <- read.table(file = "plink.rel", header = F)
grm1IDs <- read.table(file = "plink.rel.id", header = F)

subset_rowsG <- grep("UC-", grm1IDs[,1])
subset_rowsG <- subset_rowsG[!is.na(subset_rowsG)]

colnames(GRM1)[subset_rowsG] <- substr(grm1IDs[subset_rowsG,1], 9, 18)
colnames(GRM1)[-subset_rowsG] <- substr(grm1IDs[-subset_rowsG,1], 9, 21)

rownames(GRM1) <- colnames(GRM1)

GRM1 <- GRM1[bs$Id_seq,bs$Id_seq]  # Reorder columns and rows
GRM1 <- as.matrix(GRM1)

colnames(GRM1) <- seq(1,47)
rownames(GRM1) <- colnames(GRM1)

# Kinships from pedigree
ped <- ped %>%
  mutate(FID = as.factor(interaction(Sire, Dam, drop = TRUE))) %>%
  mutate(FID = as.numeric(FID)) # Convert FID to numeric for PLINK compatibility

ped$Sex <- 0
ped$Phenotype <- -9

ped_for_plink <- ped %>%
  select(FID, Id, Sire, Dam, Sex, Phenotype) %>%
  rename(IID = Id, PID = Sire, MID = Dam)

tagIDs <- tagIDs %>% rename(IID = Id_tag)

ped_for_plink <- ped_for_plink[ped_for_plink$IID %in% tagIDs$IID,]

ped_for_plink <- merge(ped_for_plink, tagIDs[,1:2], by = "IID")

ped_for_plink <- ped_for_plink[match(bs$Id_seq, ped_for_plink$Id_seq),]

score_matrix <- matrix(0, nrow = nrow(ped_for_plink), ncol = nrow(ped_for_plink))

rownames(score_matrix) <- ped_for_plink$IID
colnames(score_matrix) <- ped_for_plink$IID

for (i in seq(1,ncol(score_matrix))) {
  for (j in seq(1,ncol(score_matrix))) {
    if (i == j) {
      # Same IID
      score_matrix[i, j] <- 1
    } else {
      # Different IIDs
      if (ped_for_plink$FID[i] == ped_for_plink$FID[j]) {
        # Same FID
        score_matrix[i, j] <- 0.5
      } else if (ped_for_plink$PID[i] == ped_for_plink$PID[j] | ped_for_plink$MID[i] == ped_for_plink$MID[j]) {
        # Different FID but common PID or MID
        score_matrix[i, j] <- 0.25
      } # Otherwise, score remains 0 (different IIDs, FIDs, and no common PID or MID)
    }
  }
}


# Plots
meth_gen <- data.frame(MRM = as.numeric(MRM[lower.tri(MRM, diag = F)]),
                       A_mat = as.numeric(rA[lower.tri(rA, diag = F)]), 
                       GRM = as.numeric(GRM1[lower.tri(GRM1, diag = F)]),
                       kinship = as.numeric(score_matrix[lower.tri(score_matrix, diag = F)])
                       )

meth_gen$kinship[meth_gen$kinship==0.5] <- "full-sibs"
meth_gen$kinship[meth_gen$kinship==0.25] <- "half-sibs"
meth_gen$kinship[meth_gen$kinship==0] <- "non-sibs"

png(filename = "matrix_pairs12.png", width = 1250, height = 1200, res = 300)
ggpairs(meth_gen, columns = c(1,2), aes(color=kinship, alpha=0.9)) + 
  scale_fill_manual(values = c("turquoise4", "salmon3", "purple3")) + 
  scale_color_manual(values = c("turquoise4", "salmon3", "purple3"))
dev.off()

png(filename = "matrix_pairs13.png", width = 1250, height = 1200, res = 300)
ggpairs(meth_gen, columns = c(1,3), aes(color=kinship, alpha=0.9)) + 
  scale_fill_manual(values = c("turquoise4", "salmon3", "purple3")) + 
  scale_color_manual(values = c("turquoise4", "salmon3", "purple3"))
dev.off()

## Complex Heatmap
col1 = colorRamp2(c(min(meth_gen$A_mat), median(meth_gen$A_mat), max(meth_gen$A_mat)), c("white", "skyblue","darkblue"))
col2 = colorRamp2(c(min(meth_gen$MRM), median(meth_gen$MRM), max(meth_gen$MRM)), c("white", "lightpink","darkred"))
col3 = colorRamp2(c(min(meth_gen$GRM), median(meth_gen$GRM), max(meth_gen$GRM)), c("white", "thistle","purple4"))


ht1 = Heatmap(as.matrix(rA), rect_gp = gpar(type="none"), col=col1,
              cluster_rows = FALSE, cluster_columns = FALSE,
              cell_fun = function(j, i, x, y, w, h, fill) {
                if(i > j) {
                  grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "darkgrey", lwd=2))
                }
              }, heatmap_legend_param = list(legend_height=unit(5,"cm"), title = expression(bold('A'['MAT']))))

ht2 = Heatmap(MRM,name = "MRM" ,rect_gp = gpar(type = "none"), col = col2, column_labels = rep(" ", ncol(MRM)),
              cluster_rows = FALSE, cluster_columns = FALSE,
              cell_fun = function(j, i, x, y, w, h, fill) {
                if(i < j) {
                  grid.rect(x, y, w, h, gp = gpar(fill = fill,col = "darkgrey", lwd=2))
                }
              }, heatmap_legend_param = list(legend_height=unit(5,"cm")))

ht3 = Heatmap(GRM1,name = "GRM", rect_gp = gpar(type="none"), col=col3,
              cluster_rows = FALSE, cluster_columns = FALSE,
              cell_fun = function(j, i, x, y, w, h, fill) {
                if(i > j) {
                  grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "darkgrey", lwd=2))
                }
              }, heatmap_legend_param = list(legend_height=unit(5,"cm")))

draw(ht1 + ht2, ht_gap = unit(-300, "mm"))
draw(ht3 + ht2, ht_gap = unit(-300, "mm"))

## Tanglegram
d1 <- as.matrix(1 - rA) %>% dist() %>% hclust(method = "ward.D") %>% as.dendrogram()

d2 <- as.matrix(1 - MRM) %>% dist() %>% hclust(method = "ward.D") %>% as.dendrogram()

dl <- dendlist(d1 %>%
                 set("branches_lty", 1) %>%
                 set("branches_k_color", value = brewer.pal(12,"Paired")[c(1:10,12)], k=11),
               d2 %>%
                 set("branches_lty", 1) %>%
                 set("branches_k_color", value = brewer.pal(12,"Paired")[c(1:10,12)], k=11)
)

set.seed(3958)
x <- dl %>% untangle(method = "random", R = 10)
x %>% plot(main="Hierarchical clustering: A-matrix vs Methylation correlation matrix")
