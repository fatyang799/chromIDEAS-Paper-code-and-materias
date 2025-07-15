# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(stringr)
  library(ggplot2)
  library(reshape2)
  suppressMessages(library(GenomicFeatures))
  suppressMessages(library(ChIPseeker))
}

# load the txdb file
if (T) {
  file <- "data/raw_data/gencode.v40.annotation.sqlite"
  txdb <- loadDb(file)
  txdb
}

# load the peaks regions
if (T) {
  files <- list.files(path = "data/raw_data/8.paper/PMC8175737/histonePK", pattern = "Peak", full.names = T)
  peaks <- lapply(files, function(x) {
    # x <- files[1]
    pk <- read.table(x, header = F, sep = "\t")
    pk <- pk[pk$V8>(-log10(0.001)), ]
    
    pk <- pk[pk$V1 %in% paste0("chr",c(1:22,"X","Y")),]
    pk <- GRanges(seqnames = pk$V1,
                  ranges = IRanges(start = pk$V2+1, end = pk$V3),
                  # strand = pk$V6,
                  name = pk$V4,
                  score = pk$V5,
                  signalValue = pk$V7,
                  pValue = pk$V8,
                  qValue = pk$V9,
                  id = str_split(basename(x), "[.]", simplify = T)[, 1])
    
    return(pk)
  })
  
  names(peaks) <- str_split(basename(files), "[.]", simplify = T)[, 1]
  lengths(peaks)
}

# load the blacklist
if (T) {
  file <- "data/raw_data/2.states/hg38-blacklist.v2.bed"
  
  bl <- read.table(file, header = F, sep = "\t", fill = T, comment.char = "")
  bl <- GRanges(seqnames = bl$V1,
                ranges = IRanges(start = bl$V2+1, end = bl$V3))
  bl
}

# target region setting [TSS,TES]
if (T) {
  # 61544
  gene <- genes(txdb)
  
  # findOverlaps with bl
  if (T) {
    hit <- findOverlaps(query=bl, subject=gene, type="any")
    
    overlap_gene <- subjectHits(hit)
    bl_genes <- gene$gene_id[unique(overlap_gene)]
  }
  
  # remove gene overlap with bl (59413)
  gene <- gene[! gene$gene_id %in% bl_genes]
  
  rm(bl, overlap_gene, bl_genes, hit)
}

# findOverlaps
if (T) {
  peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb, 
                         tssRegion=c(-1000, 1000), verbose=FALSE)
  
  plotAnnoBar(peakAnnoList)
}