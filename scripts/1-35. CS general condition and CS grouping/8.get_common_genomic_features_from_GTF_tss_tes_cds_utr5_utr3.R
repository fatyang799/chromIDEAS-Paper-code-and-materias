# 根据GENCODE的GTF文件对GenomicFeatures进行处理

###### definition for regions ######
# bin_tss: TSS±200bp
# bin_tes: TES±200bp
# bin_cds: CDS
# bin_utr5: UTR5
# bin_utr3: UTR3
# bin_exon: EXON
# bin_intron: INTRON
# bin_intergenic: exclude region from [tss-400, tes+400]


# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(stringr)
  library(reshape2)
  suppressPackageStartupMessages(library(GenomicFeatures))
  library(qs)
}

# load the genome windows bin
if (T) {
  file <- "data/saved_data/1.windowsNoBlack.withid.qs"
  bin <- qread(file, nthreads = 6)
}

# construct the gtf annotation txdb object
if (T) {
  file <- "data/raw_data/gencode.v40.annotation.sqlite"
  txdb <- loadDb(file)
  
  # only focus on the chr1-22.XY
  seqlevels(txdb) <- c(paste0("chr", 1:22), "chrX", "chrY")
  
  # get all tx
  tx <- transcripts(txdb, columns=c("TXNAME", "GENEID"))
}

# get location site
if (T) {
  chr <- seqnames(tx)
  start <- start(tx)
  end <- end(tx)
  strand <- as.vector(strand(tx))
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

# construct the gtf related genomic feature bin: TSS±200bp
if (T) {
  file <- "data/saved_data/8.common_gtf_structure_tss_quality_consistent.bin.qs"
  
  if (file.exists(file)) {
    mess <- "The result of TSS in GTF has existed."
    print(mess)
  }
  if (!file.exists(file)) {
    tss <- GRanges(seqnames = chr, 
                   ranges = IRanges(start = ifelse(strand == "+", start-200, end-200),
                                    end = ifelse(strand == "+", start+200, end+200)))
    bin_tss <- overID(tss, bin)
    
    qsave(bin_tss, file = file, nthreads = 6)
    
    mess <- "The result of TSS in GTF has completed."
    print(mess)
  }
}

# construct the gtf related genomic feature bin: TES±200bp
if (T) {
  file <- "data/saved_data/8.common_gtf_structure_tes_quality_consistent.bin.qs"
  
  if (file.exists(file)) {
    mess <- "The result of TES in GTF has existed."
    print(mess)
  }
  if (!file.exists(file)) {
    tes <- GRanges(seqnames = chr, 
                   ranges = IRanges(start = ifelse(strand == "+", end-200, start-200),
                                    end = ifelse(strand == "+", end+200, start+200)))
    bin_tes <- overID(tes, bin)
    
    qsave(bin_tes, file = file, nthreads = 6)
    
    mess <- "The result of TES in GTF has completed."
    print(mess)
  }
}

# construct the gtf related genomic feature bin: CDS
if (T) {
  file <- "data/saved_data/8.common_gtf_structure_cds_quality_consistent.bin.qs"
  
  if (file.exists(file)) {
    mess <- "The result of CDS in GTF has existed."
    print(mess)
  }
  if (!file.exists(file)) {
    cds <- cds(txdb, columns=c("TXNAME", "GENEID"))
    bin_cds <- overID(cds, bin)
    
    qsave(bin_cds, file = file, nthreads = 6)
    
    mess <- "The result of CDS in GTF has completed."
    print(mess)
  }
}

# construct the gtf related genomic feature bin: utr5
if (T) {
  file <- "data/saved_data/8.common_gtf_structure_utr5_quality_consistent.bin.qs"
  
  if (file.exists(file)) {
    mess <- "The result of UTR5 in GTF has existed."
    print(mess)
  }
  if (!file.exists(file)) {
    utr5 <- fiveUTRsByTranscript(txdb)
    utr5 <- unlist(utr5)
    utr5 <- granges(utr5)
    utr5
    
    bin_utr5 <- overID(utr5, bin) 
    
    qsave(bin_utr5, file = file, nthreads = 6)
    
    mess <- "The result of UTR5 in GTF has completed."
    print(mess)
  }
}

# construct the gtf related genomic feature bin: utr3
if (T) {
  file <- "data/saved_data/8.common_gtf_structure_utr3_quality_consistent.bin.qs"
  
  if (file.exists(file)) {
    mess <- "The result of UTR3 in GTF has existed."
    print(mess)
  }
  if (!file.exists(file)) {
    utr3 <- threeUTRsByTranscript(txdb)
    utr3 <- unlist(utr3)
    utr3 <- granges(utr3)
    utr3
    
    bin_utr3 <- overID(utr3, bin)
    
    qsave(bin_utr3, file = file, nthreads = 6)
    
    mess <- "The result of UTR3 in GTF has completed."
    print(mess)
  }
}

# construct the gtf related genomic feature bin: exon
if (T) {
  file <- "data/saved_data/8.common_gtf_structure_exon_quality_consistent.bin.qs"
  
  if (file.exists(file)) {
    mess <- "The result of Exon in GTF has existed."
    print(mess)
  }
  if (!file.exists(file)) {
    exon <- exons(txdb)
    exon
    
    bin_exon <- overID(exon, bin)
    
    qsave(bin_exon, file = file, nthreads = 6)
    
    mess <- "The result of Exon in GTF has completed."
    print(mess)
  }
}

# construct the gtf related genomic feature bin: intron
if (T) {
  file <- "data/saved_data/8.common_gtf_structure_intron_quality_consistent.bin.qs"
  
  if (file.exists(file)) {
    mess <- "The result of Intron in GTF has existed."
    print(mess)
  }
  if (!file.exists(file)) {
    intron <- intronsByTranscript(txdb, use.names=T)
    intron <- unlist(intron)
    intron <- granges(intron)
    intron
    
    bin_intron <- overID(intron, bin)
    
    qsave(bin_intron, file = file, nthreads = 6)
    
    mess <- "The result of Intron in GTF has completed."
    print(mess)
  }
}

# construct the gtf related genomic feature bin: intergenic
if (T) {
  file <- "data/saved_data/8.common_gtf_structure_intergenic_quality_consistent.bin.qs"
  
  if (file.exists(file)) {
    mess <- "The result of Intergenic in GTF has existed."
    print(mess)
  }
  if (!file.exists(file)) {
    gene <- genes(txdb)
    gene_ext <- gene+400
    
    bin_gene_ext <- overID(gene_ext, bin)
    bin_intergenic <- setdiff(bin$ID, bin_gene_ext)
    
    qsave(bin_gene_ext, file = file, nthreads = 6)
    
    mess <- "The result of Intergenic in GTF has completed."
    print(mess)
  }
}
