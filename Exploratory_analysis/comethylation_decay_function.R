comethylation_decay <- function(bin_widths = seq(10000, 500000, by=10000),
                                      nregions = 1000,
                                      length_sd = 0,
                                      genome,
                                      m,
                                      r,
                                      seed = 3789) {
  
  # Initialize a data frame to store results
  cors <- data.frame(cor = numeric(), abs_cor = numeric(), width = integer(), n = as.integer())
  
  # Set seed for reproducibility
  set.seed(seed)
    
  # Loop over each bin width
  for (j in bin_widths) {
      
    # Generate random regions
    if (j <= 300) {
      
      bins <- regioneR::createRandomRegions(nregions = 3*nregions,
                                            length.mean = j,
                                            length.sd = length_sd,
                                            genome = genome)
      
    } else {
      
      bins <- regioneR::createRandomRegions(nregions = nregions,
                                            length.mean = j,
                                            length.sd = length_sd,
                                            genome = genome)
      
    }
    
    bins <- as.data.frame(bins)
    bins$seqnames <- as.character(bins$seqnames)
      
    # Loop over each region
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
        
      # Append the results to the data frame
      cors <- rbind(cors, data.frame(cor = c, abs_cor = abs_c, width = w, n = n))
    }
    }
  
  # Filter out incomplete cases and adjust width scale to Kb
  cors <- cors[complete.cases(cors),]
  cors$width <- cors$width / 1000
  
  return(cors)
}