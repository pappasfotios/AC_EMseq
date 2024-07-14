library(data.table)
library(GenomicRanges)
library(dplyr)
library(ggplot2)

# Import data
bs <- readRDS("noSNPs_BSseq.rds")
genome <- read.table(file.choose(), header = F)
genome$start <- rep(1,nrow(genome))
names(genome)[1:2] <- c("chr", "end")
genome <- genome[grep("NC", genome$chr),]
genome <- genome[genome$chr != "NC_000861.1",]

genome <- makeGRangesFromDataFrame(genome)

m <- bsseq::getMeth(bs, type = "raw", what = "perBase")
r <- as.data.frame(granges(bs))
r$seqnames <- as.character(r$seqnames)

cov <- bsseq::getCoverage(bs, type = "Cov")

SDmask <- matrixStats::rowSds(m) > 0.1
COVmask <- matrixStats::rowMins(cov) >= 10

r <- r[complete.cases(m) & SDmask & COVmask,]
r <- as.data.frame(r)
m <- m[complete.cases(m) & SDmask & COVmask,]
  
bin_widths <- seq(10000, 1000000, by=10000)

cors <- data.frame(cor = numeric(), abs_cor = numeric(), width = integer(), n = as.integer())

for (j in bin_widths) {
  
  bins <- regioneR::createRandomRegions(nregions = 1000,
                                        length.mean = j,
                                        length.sd = 0,
                                        genome = genome)
  bins <- as.data.frame(bins)
  bins$seqnames <- as.character(bins$seqnames)
  
  for (i in seq(1, nrow(bins))) {
    meth <- m[r$seqnames == bins$seqnames[i] & r$start >= bins$start[i] & r$end <= bins$end[i], ]
    w <- bins$width[i]
    
    if (!is.null(meth) && !is.null(dim(meth)) && !is.na(w) && w > 0) {
      if (!is.na(dim(meth)[1]) && dim(meth)[1] > 5 && dim(meth)[1] > w / 5000) {
        cometh <- cor(t(meth))
        c <- mean(cometh[lower.tri(cometh, diag = FALSE)])
        abs_c <- mean(abs(cometh[lower.tri(cometh, diag = FALSE)]))
        n <- ncol(cometh)
      } else {
        c <- NA
        abs_c <- NA
        n <- NA
      }
    } else {
      c <- NA
      abs_c <- NA
      n <- NA
    }
    
    cors <- rbind(cors, data.frame(cor = c, abs_cor = abs_c, width = w, n = n))
  }
}


cors <- cors[complete.cases(cors),]
cors$width <- cors$width / 1000

reduced <- setDT(cors)[, .(value= mean(abs_cor)),by=width]

png(filename = "AbsCometh.png", width = 1200, height = 1000, res=300)
ggplot(data=reduced, aes(x=width, y=value)) + theme_light() +
  geom_smooth(col="darkred", lwd=1) + 
  ylab("Absolute correlation") + 
  xlab("Bin width")
dev.off()

png(filename = "AbsCometh.png", width = 1600, height = 1000, res=300)
ggplot(data=cors, aes(x=width, y=cor, group=width)) + 
  geom_boxplot(fill="#E30B5D") + xlab("Bin width (kb)") + ylab("Correlation") +
  theme_minimal()
dev.off()


