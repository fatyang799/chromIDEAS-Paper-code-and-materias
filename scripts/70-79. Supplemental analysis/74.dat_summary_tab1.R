# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(stringr)
  library(reshape2)
}

# read tab dat
if (T) {
  tab <- read.table("data/raw_data/10.data_summary/Table S1.data summary.csv", header = T, sep = ",", fill = T, comment.char = "")
  tab <- tab[tab$Database == "ENCODE", ]
}

# for fq file
if (T) {
  # read metadata
  if (T) {
    meta <- read.table("data/raw_data/10.data_summary/file_correspondence.txt", header = T, sep = "\t", fill = T, comment.char = "")
    meta <- data.frame(sapply(meta, basename))
  }
  
  # merge dat
  if (T) {
    colnames(tab)
    colnames(meta)
    
    tab$exp_meta <- meta$Source[match(tab$Exp_SourceFqFile, meta$Target)]
    tab$ct_meta <- meta$Source[match(tab$CT_SourceFqFile, meta$Target)]
  }
}

# for fq file from bam file
if (T) {
  # read metadata
  if (T) {
    meta <- read.table("data/raw_data/10.data_summary/9.change_name_detail_info.txt", header = T, sep = "\t", fill = T, comment.char = "")
    meta <- data.frame(sapply(meta, basename))
  }
  
  # fq file from bam format
  if (T) {
    tab$bam_fq_exp <- ifelse(grepl("no_fq_ID", tab$Exp_SourceFqFile), basename(tab$Exp_SourceFqFile), NA)
    tab$bam_fq_ct <- ifelse(grepl("no_fq_ID", tab$CT_SourceFqFile), basename(tab$CT_SourceFqFile), NA)
  }
  
  # merge dat
  if (T) {
    colnames(tab)
    colnames(meta)
    
    tab$exp_meta_bam <- meta$Source[match(tab$bam_fq_exp, meta$Target)]
    tab$ct_meta_bam <- meta$Source[match(tab$bam_fq_ct, meta$Target)]
  }
}

# merge the info
if (T) {
  merged <- tab[, 1:3]
  
  merged$exp <- ifelse(is.na(tab$exp_meta), tab$exp_meta_bam, tab$exp_meta)
  merged$ct <- ifelse(is.na(tab$ct_meta), tab$ct_meta_bam, tab$ct_meta)
}

# format the data 
if (T) {
  head(merged)
  
  merged$exp_ID <- str_extract(merged$Exp_SourceFqFile, "ENCSR.{6}")
  merged$ct_ID <- str_extract(merged$CT_SourceFqFile, "ENCSR.{6}")
  
  merged$exp_rep <- str_extract(merged$exp, "rep[0-9]{1,2}_[0-9]{1,2}")
  merged$ct_rep <- str_extract(merged$ct, "rep[0-9]{1,2}_[0-9]{1,2}")
  
  table(merged$exp_ID, merged$exp_rep)
}

write.table(merged, file = "results/1.tab/74.dat_summary_ENCODE.csv", 
            quote = F, sep = ",", col.names = T, row.names = F)
