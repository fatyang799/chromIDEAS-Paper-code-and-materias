# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(reshape2)
  library(qs)
}

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
body_bin_num <- 10
rna_levels <- 10
length_leveles <- 10
state_file <- "data/raw_data/2.states/chromIDEAS.state"

# get gene body metadata: gene_body_ID
if (T) {
  file <- "data/saved_data/20.geneBody_metadata_tx_level.qs"
  gene_body_ID <- qread(file)
}

# manual
if (F) {
  ## gene_body_mat[, name]:
  ##     gene_body_ID$strand == "+": 
  ##         gene_body_mat$unit>=1: 
  ##             paste0(
  ##                 pos_strand_rounding(gene_body_ID$tss_BinID + (bin_id -1) * gene_body_mat$unit), 
  ##                 "-", 
  ##                 pos_strand_rounding(gene_body_ID$tss_BinID + (bin_id) * gene_body_mat$unit - 1)
  ##             )
  ##         else: 
  ##             paste0(
  ##                 floor(gene_body_ID$tss_BinID + (bin_id -1) * gene_body_mat$unit + 1e-5), 
  ##                 "-", 
  ##                 floor(gene_body_ID$tss_BinID + (bin_id) * gene_body_mat$unit - 1e-5)
  ##             )
  ##     else: 
  ##         gene_body_mat$unit>=1: 
  ##             paste0(
  ##                 neg_strand_rounding(gene_body_ID$tes_BinID + (body_bin_num+1-(bin_id) -1) * gene_body_mat$unit + 1), 
  ##                 "-", 
  ##                 neg_strand_rounding(gene_body_ID$tes_BinID + (body_bin_num+1-(bin_id)) * gene_body_mat$unit)
  ##             )
  ##         else: 
  ##             paste0(
  ##                 floor(gene_body_ID$tes_BinID + (body_bin_num+1-(bin_id) -1) * gene_body_mat$unit + 1e-5), 
  ##                 "-", 
  ##                 floor(gene_body_ID$tes_BinID + (body_bin_num+1-(bin_id)) * gene_body_mat$unit - 1e-5)
  ##             )
}

# Tx up down matrix
if (T) {
  head(gene_body_ID)
  gene_body_mat <- data.frame(gene_id = gene_body_ID$tx_id)
  gene_body_mat$strand <- gene_body_ID$strand
  gene_body_mat$unit <- ifelse(gene_body_ID$strand == "+", 
                               (gene_body_ID$tes_BinID - gene_body_ID$tss_BinID +1) / body_bin_num, 
                               (gene_body_ID$tss_BinID - gene_body_ID$tes_BinID +1) / body_bin_num)
  
  # define rounding function
  if (T) {
    pos_strand_rounding <- function(x) {
      x_nextL <- x*10
      remaining <- x_nextL %% 10
      
      res <- ifelse(remaining>=5, ceiling(x), floor(x))
      
      return(res)
    }
    neg_strand_rounding <- function(x) {
      x_nextL <- x*10
      remaining <- x_nextL %% 10
      
      res <- ifelse(remaining<=5, floor(x)-1, ceiling(x)-1)
      
      return(res)
    }
  }
  
  for (bin_id in 1:body_bin_num) {
    # bin_id <- 1
    name <- paste0("G", bin_id)
    print(name)
    
    gene_body_mat[, name] <- ifelse(gene_body_ID$strand == "+", 
                                    ifelse(gene_body_mat$unit>=1, 
                                           paste0(
                                             pos_strand_rounding(gene_body_ID$tss_BinID + (bin_id -1) * gene_body_mat$unit), 
                                             "-", 
                                             pos_strand_rounding(gene_body_ID$tss_BinID + (bin_id) * gene_body_mat$unit - 1)
                                           ), 
                                           paste0(
                                             floor(gene_body_ID$tss_BinID + (bin_id -1) * gene_body_mat$unit + 1e-5), 
                                             "-", 
                                             floor(gene_body_ID$tss_BinID + (bin_id) * gene_body_mat$unit - 1e-5)
                                           )), 
                                    ifelse(gene_body_mat$unit>=1, 
                                           paste0(
                                             neg_strand_rounding(gene_body_ID$tes_BinID + (body_bin_num+1-(bin_id) -1) * gene_body_mat$unit + 1), 
                                             "-", 
                                             neg_strand_rounding(gene_body_ID$tes_BinID + (body_bin_num+1-(bin_id)) * gene_body_mat$unit)
                                           ), 
                                           paste0(
                                             floor(gene_body_ID$tes_BinID + (body_bin_num+1-(bin_id) -1) * gene_body_mat$unit + 1e-5), 
                                             "-", 
                                             floor(gene_body_ID$tes_BinID + (body_bin_num+1-(bin_id)) * gene_body_mat$unit - 1e-5)
                                           )))
  }
  
  rm(bin_id, name, pos_strand_rounding, neg_strand_rounding)
}

# read the states data
if (T) {
  state <- data.table::fread(state_file, sep = " ", header = T, data.table = F)
  state <- state[, c(1, 5, 6)]
  colnames(state)[1] <- "ID"
  head(state)
}

# get state percentage matrix
if (T) {
  state_profile <- function(gene_body_mat, state, col) {
    # col <- 2
    # hello info
    if (T) {
      start_mess <- paste0("|", paste(rep("-", 100), collapse = ""), "|\n")
      cat(start_mess)
    }
    
    # prepare for calculation
    if (T) {
      gene_body_mat$num <- 1:nrow(gene_body_mat)
      state_order <- sort(as.numeric(unique(state[, col])))
      
      breaks <- round(seq(min(gene_body_mat$num), max(gene_body_mat$num), length.out=100))
    }
    
    # get info for each region: percentage of each state
    if (T) {
      dat <- apply(gene_body_mat, 1, function(x) {
        # x <- unlist(gene_body_mat[1, ])
        if (as.numeric(x[length(x)]) %in% breaks) {
          if (which(as.numeric(x[length(x)]) == breaks) == 1) {
            cat("|*")
          }
          if (which(as.numeric(x[length(x)]) == breaks) == 100) {
            cat("*|\n")
          }
          if (! which(as.numeric(x[length(x)]) == breaks) %in% c(1, 100)) {
            cat("*")
          }
        }
        
        # genebody
        if (T) {
          genebody_bin <- x[grepl("^G[0-9]{1,2}$", names(x))]
          genebody <- data.frame(t(
            sapply(genebody_bin, function(gb) {
              # gb <- genebody_bin[1]
              start <- as.numeric(strsplit(gb, "-")[[1]][1])
              end <- as.numeric(strsplit(gb, "-")[[1]][2])
              s_dat <- state[start:end, col]
              
              s_p <- sapply(state_order, function(s) {
                # norm the number of bins
                sum(s_dat==s) / length(s_dat)
              })
              names(s_p) <- paste0("S", state_order)
              return(s_p)
            })
          ))
          rownames(genebody) <- paste0(x[1], "@", rownames(genebody))
        }
        
        # return result
        return(genebody)
      })
      dat <- data.frame(do.call(rbind, dat))
    }
    
    return(dat)
  }
  
  head(state)
  
  for (cell in c(cell1, cell2)) {
    file <- paste0("data/saved_data/21.tx_Body(", body_bin_num, "parts)_percentage_based_on_Body(", length_leveles, 
                   "parts)_RNA(", rna_levels, "parts).bs", bin_size, ".", cell, ".qs")
    if (! file.exists(file)) {
      n <- grep(cell, colnames(state))
      dat <- state_profile(gene_body_mat, state, n)
      qsave(dat, file, nthreads = 6)
    }
  }
}


# line1-23: 读入tss_bin以及tes_bin数值（变量gene_body_ID），用于制作gene_body_mat
# line24-42: 基因切块逻辑说明
# line43-70: 在gene_body_mat中记录了当基因分成10等份后，基因每份中bin id的编号。编号具体计算方法如下：
#                1-9部分：【ceiling(tss + (tes-tss) / body_bin_num * (gene_part-1)) 】 - 【ceiling(tss + (tes-tss) / body_bin_num * (gene_part)) -1】
#                10部分：【ceiling(tss + (tes-tss) / body_bin_num * (gene_part-1)) 】 - 【ceiling(tss + (tes-tss) / body_bin_num * (gene_part))】
#                **注意**：为避免重复计算问题，仅最后一个区间为首尾均闭合，其他区间均为【）模式
# line71-78: 读入cd34和thp1的chromatin states数据
# line79-154: 根据gene_body_mat数据，可以计算细胞在每个基因的每个部分的chromatin states比例
