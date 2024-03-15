setwd("C:\\Users\\fopa0001\\Downloads\\OneDrive_1_04-12-2023")

library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(RColorBrewer)
library(pedigree)
library(ggplot2)

bs <- readRDS("Filtered_BSseq.rds")

# ERM
m <- bsseq::getMeth(bs, type = "raw", what = "perBase")
cov <- bsseq::getCoverage(bs, type = "Cov")

SDmask <- matrixStats::rowSds(m) > 0.1
COVmask <- matrixStats::rowMins(cov) >= 10

m <- m[complete.cases(m) & SDmask & COVmask,]

erm1 <- cor(m)
erm2 <- scale(t(m)) %*% t(scale(t(m)))
erm2 <- erm2/mean(diag(erm2))

colnames(erm1) <- seq(1,47)
rownames(erm1) <- seq(1,47)

corrplot::corrplot(erm2, method = "square", diag = T, is.corr = F)

ped <- read.table(file.choose(), header = T, stringsAsFactors = T)
ped <- ped[!duplicated(ped$Id),]

tagIDs <- openxlsx::read.xlsx(file.choose())
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
subset_rows <- grep("UC-", row.names(bs@colData))
subset_rows <- subset_rows[!is.na(subset_rows)]

bs@colData$Id_seq[subset_rows] <- substr(rownames(bs@colData)[subset_rows], 9, 18)
bs@colData$Id_seq[-subset_rows] <- substr(rownames(bs@colData)[-subset_rows], 9, 21)
colData <- dplyr::inner_join(as.data.frame(bs@colData), tagIDs[,1:2], by="Id_seq")
colData$key <- seq(1,47)

row.names(rA) <- colData$key[colData$Id_tag == names(rA)]
names(rA) <- row.names(rA)

A <- reshape::melt(rA,)

corrplot::corrplot(as.matrix(rA), method = "square", diag = F, is.corr = F, type = "lower")

meth_gen <- data.frame(ERM = as.numeric(erm2[lower.tri(erm2, diag = F)]), A_mat = as.numeric(rA[lower.tri(rA, diag = F)]))
ggplot(data = meth_gen, aes(x=A_mat, y=ERM)) + geom_point(color="magenta3", pch=19) + geom_smooth(method=lm , color="black", fill="skyblue3", se=TRUE) + theme_light()

## Complex Heatmap
col1 = colorRamp2(c(min(meth_gen$A_mat), median(meth_gen$A_mat), max(meth_gen$A_mat)), c("white", "skyblue","darkblue"))
col2 = colorRamp2(c(min(meth_gen$ERM), median(meth_gen$ERM), max(meth_gen$ERM)), c("white", "lightpink","darkred"))

ht1 = Heatmap(as.matrix(rA),name = "A-matrix", rect_gp = gpar(type="none"), col=col1,
              cluster_rows = FALSE, cluster_columns = FALSE,
              cell_fun = function(j, i, x, y, w, h, fill) {
                if(i > j) {
                  grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "darkgrey", lwd=2))
                }
              })

ht2 = Heatmap(erm2,name = "MRM" ,rect_gp = gpar(type = "none"), col = col2, column_labels = rep(" ", ncol(erm1)) ,
              cluster_rows = FALSE, cluster_columns = FALSE,
              cell_fun = function(j, i, x, y, w, h, fill) {
                if(i < j) {
                  grid.rect(x, y, w, h, gp = gpar(fill = fill,col = "darkgrey", lwd=2))
                }
              })

draw(ht1 + ht2, ht_gap = unit(-300, "mm"))

## Tanglegram
d1 <- as.matrix(1 - rA) %>% dist() %>% hclust(method = "ward.D") %>% as.dendrogram()

d2 <- as.matrix(1 - erm1) %>% dist() %>% hclust(method = "ward.D") %>% as.dendrogram()

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
