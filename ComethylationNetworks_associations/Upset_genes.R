setwd("//wsl.localhost/Ubuntu/home/fotis/analysis/EMseq/enrichment/")

library(ComplexHeatmap)
library(circlize)

matrix_build <- function(db="hs") {
  
  suffix <- paste0("go",db,"/",db,"_proteins")
  
  promoters_mot <- read.table(paste0("./promoters/mot/", suffix), header = F)
  promoters_den <- read.table(paste0("./promoters/den/", suffix), header = F)
  Islands_mot_pos <- read.table(paste0("./IslandsPlus/mot/pos/", suffix), header = F)
  Islands_mot_neg <- read.table(paste0("./IslandsPlus/mot/neg/", suffix), header = F)
  Islands_den <- read.table(paste0("./IslandsPlus/den/", suffix), header = F)
  firstIntrons_mot <- read.table(paste0("./first_intron/", suffix), header = F)
  
  list <- list(Promoters_motility = promoters_mot$V1, 
               Promoters_concentration = promoters_den$V1,
               Islands_motility = c(Islands_mot_pos$V1, Islands_mot_neg$V1),
               Islands_concentration = Islands_den$V1,
               FirstIntrons_motility = firstIntrons_mot$V1)
  
  return(make_comb_mat(list, mode = "intersect"))
  
}

m_zf <- matrix_build(db="zf")
m_hs <- matrix_build()

UpSet(m_zf[comb_degree(m_zf) >= 2], comb_order = order(comb_size(m_zf[comb_degree(m_zf) >= 2])))
UpSet(m_hs[comb_degree(m_hs) >= 2], comb_order = order(comb_size(m_hs[comb_degree(m_hs) >= 2])))
