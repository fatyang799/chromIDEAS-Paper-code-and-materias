# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(qs)
}

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
up_bin_num <- 5
down_bin_num <- 5
body_bin_num <- 10

# get all gene body bin id matrix
if (T) {
  data_type <- "tx"
  
  gene_body_mat <- qread(paste0("data/saved_data/11.", data_type, "_level_profile_Body_", up_bin_num, "UP_", 
                                down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".bin.qs"), nthreads = 6)
}

# get non_overlap-txs ID
if (T) {
  file <- "data/saved_data/20.txs_with_nonoverlap_region.qs"
  
  non_overlap_txs <- qread(file, nthreads = 6)
}

# get non_overlap-txs body profile ID
if (T) {
  file <- paste0("data/saved_data/47.", data_type, "s_with_nonoverlap_region_matrix_bs", bin_size, ".qs")
  
  if (! file.exists(file)) {
    over <- intersect(non_overlap_txs$tx_id, gene_body_mat$gene_id)
    
    gene_body_mat <- gene_body_mat[match(over, gene_body_mat$gene_id), ]
    
    qsave(gene_body_mat, file, nthreads = 6)
  }
}
