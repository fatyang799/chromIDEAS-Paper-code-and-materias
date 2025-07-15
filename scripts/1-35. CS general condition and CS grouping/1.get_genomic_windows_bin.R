# 以IDEAS的output为标准
# 得到所有bin的坐标
# 由于IDEAS的bin坐标是bed格式，故坐标是0-base的
# 在GRange对象中是1-based，故将start+1转变

# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(qs)
  suppressPackageStartupMessages(library(GenomicFeatures))
}

# construct the genome windows bin
if (T) {
  file <- "data/saved_data/1.windowsNoBlack.withid.qs"
  
  if (! file.exists(file)) {
    dat <- data.table::fread(file = "data/raw_data/2.states/chromIDEAS.state", header = T, sep = " ", fill = T, skip = "", data.table=F)
    head(dat)
    
    bin <- GRanges(seqnames = dat[,2],
                   ranges = IRanges(start = dat[,3]+1,
                                    end = dat[,4]),
                   ID = dat[,1])
    qsave(bin, file = file, nthreads = 6)
  }
}
