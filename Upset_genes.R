setwd("//wsl.localhost/Ubuntu/home/fotis/analysis/EMseq/enrichment/")

library(ComplexHeatmap)
library(circlize)

matrix_build <- function(db="hs") {
  
  suffix <- paste0("go",db,"/",db,"_proteins")
  
  promoters_mot <- read.table(paste0("./promoters/1kb/vcl/", suffix), header = F)
  promoters_den <- read.table(paste0("./promoters/1kb/density/", suffix), header = F)
  IslandsShores_mot <- read.table(paste0("./IslandsPlus/", suffix), header = F)
  firstIntrons_mot <- read.table(paste0("./first_intron/", suffix), header = F)
  
  list <- list(Promoters_velocity = promoters_mot$V1, 
               Promoters_density = promoters_den$V1,
               IslandsPlusShores_velocity = IslandsShores_mot$V1,
               FirstIntrons_velocity = firstIntrons_mot$V1)
  
  return(make_comb_mat(list, mode = "intersect"))
  
}

m_zf <- matrix_build(db="zf")
m_hs <- matrix_build()

UpSet(m_zf[comb_degree(m_zf) >= 2], comb_order = order(comb_size(m_zf[comb_degree(m_zf) >= 2])))
UpSet(m_hs[comb_degree(m_hs) >= 2], comb_order = order(comb_size(m_hs[comb_degree(m_hs) >= 2])))
