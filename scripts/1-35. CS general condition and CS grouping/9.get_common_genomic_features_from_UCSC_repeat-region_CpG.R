# 根据UCSC中数据库对rmsk区域进行处理


# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  suppressPackageStartupMessages(library(qs))
  suppressPackageStartupMessages(library(GenomicFeatures))
}

# load the genome windows bin
if (T) {
  file <- "data/saved_data/1.windowsNoBlack.withid.qs"
  bin <- qread(file, nthreads = 6)
}

# define function to calculate overlaped bin ID
if (T) {
  overID <- function(target_region, bin) {
    # test data
    if (F) {
      target_region <- tss
    }
    
    # find overlap
    if (T) {
      overlap <- findOverlaps(target_region, bin, type="any", ignore.strand=T)
    }
    
    # format the results
    if (T) {
      bin_tss <- bin$ID[unique(subjectHits(overlap))]
      bin_tss <- sort(bin_tss)
    }
    
    return(bin_tss)
  }
}

# the repeat element
if (T) {
  # format dat
  if (T) {
    file <- "data/saved_data/9.ucsc_common_repeats_format_dat.qs"
    
    if (! file.exists(file)) {
      input <- "data/raw_data/UCSC/rmsk_hg38.tsv.gz"
      
      # read the data
      if (T) {
        dat <- data.table::fread(input, header = T, sep = "\t", data.table = F)
        head(dat)  
      }
      
      # only stay chr1-22, XY
      if (T) {
        unique(dat$genoName)
        dat <- dat[dat$genoName %in% paste0("chr", c(1:22, "X", "Y")), ]
      }
      
      # data filter based on alignment swScore
      if (T) {
        # Smith Waterman alignment score
        summary(dat$swScore)
        boxplot(dat$swScore, outline = F)
        plot(density(dat$swScore[dat$swScore<100]))
        
        cutoff <- quantile(dat$swScore, 0.2)
        dat <- dat[dat$swScore > cutoff, ]
      }
      
      # select column
      if (T) {
        length(unique(dat$repName))# 1741
        length(unique(dat$repClass))# 20
        length(unique(dat$repFamily))# 60
        dat <- dat[, c(6:8,10:13)]
        head(dat)  
      }
      
      # filter out low confident region
      if (T) {
        unique(dat$repClass)
        dat <- dat[!grepl("?", dat$repClass, fixed = T), ]
        dat <- dat[!grepl("?", dat$repFamily, fixed = T), ]
        unique(dat$repClass)
      }
      
      # save the data
      qsave(dat, file, nthreads = 6)
    }
    if (file.exists(file)) {
      dat <- qread(file, nthreads = 6)  
    }
  }
  
  file <- "data/saved_data/9.common_ucsc_structure_repeats_consistent.bin.qs"
  if (! file.exists(file)) {
    # construct GRanges object
    if (T) {
      head(dat)
      sapply(dat[, 5:7], unique)
      
      dat <- GRanges(seqnames = dat$genoName,
                     ranges = IRanges(start = dat$genoStart+1,
                                      end = dat$genoEnd),
                     type = dat$repClass)
    }
    
    # overlap with bin
    if (T) {
      dat <- lapply(split(dat, dat$type), function(x) {
        overID(x, bin)
      })
    }
    
    # convert into dataframe
    if (T) {
      n_num <- lengths(dat)
      name <- names(dat)
      dat <- data.frame(ID = unlist(dat), 
                        Type = rep(name, n_num))
      rownames(dat) <- 1:nrow(dat)
    }
    
    qsave(dat, file, nthreads = 6)
  }
}

# the CpG island
if (T) {
  file <- "data/saved_data/9.common_ucsc_structure_cpgislandext_quality_consistent.bin.qs"
  
  if (!file.exists(file)) {
    input <- "data/raw_data/UCSC/cpgIslandExt_hg38.txt.gz"
    dat <- data.table::fread(input, header = T, sep = "\t", data.table = F)
    head(dat)
    summary(dat$obsExp)
    cutoff <- quantile(dat$obsExp, 0.2)
    dat <- dat[dat$obsExp>cutoff, ]
    
    dat <- GRanges(seqnames = dat$chrom,
                   ranges = IRanges(start = dat$chromStart+1,
                                    end = dat$chromEnd))
    dat <- overID(dat, bin)
    qsave(dat, file, nthreads = 6)
  }
}
