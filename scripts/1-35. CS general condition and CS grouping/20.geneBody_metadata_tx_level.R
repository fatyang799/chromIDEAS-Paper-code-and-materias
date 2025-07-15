# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(reshape2)
  library(qs)
  suppressPackageStartupMessages(library(GenomicFeatures))
}

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
body_bin_num <- 10
rna_levels <- 10
length_leveles <- 10
txdb_sqlite <- "data/raw_data/gencode.v40.annotation.sqlite"
windowID <- "data/saved_data/1.windowsNoBlack.withid.qs"

# get txs not overlap with other txs
if (T) {
  over_cutoff <- 0.1
  file <- "data/saved_data/20.txs_with_nonoverlap_region.qs"
  if (! file.exists(file)) {
    # construct the gtf annotation txdb object
    if (T) {
      txdb <- loadDb(txdb_sqlite)
      
      # only focus on the chr1-22.XY
      seqlevels(txdb) <- c(paste0("chr", 1:22), "chrX", "chrY")
      
      tx <- transcripts(txdb, columns=c("TXNAME", "GENEID"))
    }
    
    # statistics: the number of tx per gene
    if (T) {
      # get gene id
      if (T) {
        table(sapply(tx$GENEID, length))
        tx$gene_id <- sapply(tx$GENEID, c)
        
        gene_tx_num_stat <- data.frame(geneid=tx$gene_id, 
                                       txid=tx$TXNAME)
      }
      
      # statistics
      if (T) {
        # data statisitcs
        if (T) {
          head(gene_tx_num_stat)
          ggdat <- as.data.frame.table(table(table(gene_tx_num_stat$geneid)))
          colnames(ggdat) <- c("N_txs", "N_genes")
          ggdat$N_txs <- as.numeric(ggdat$N_txs)
          
          # separate into <20 and >=20 txs
          if (T) {
            ggdat1 <- ggdat[ggdat$N_txs<20, ]
            ggdat2 <- sum(ggdat$N_genes[ggdat$N_txs>=20])
            ggdat <- rbind(ggdat1, 
                           data.frame(
                             N_txs = ">=20", 
                             N_genes = ggdat2
                           ))
            ggdat$N_txs <- factor(ggdat$N_txs, levels = c(1:19, ">=20"))
            rm(ggdat1, ggdat2)
          }
          
          # convert into percentage
          if (T) {
            head(ggdat)
            ggdat$N_genes_p <- ggdat$N_genes/sum(ggdat$N_genes)*100
            ggdat$N_genes_inter <- round(ggdat$N_genes/sum(ggdat$N_genes)*100, 1)
          }
        }
        
        # number of txs for each gene plot
        if (T) {
          p <- ggplot(ggdat) +
            geom_bar(aes(x=N_txs, y=N_genes_p, fill=N_txs), stat = "identity") +
            geom_text(aes(x=N_txs, y=N_genes_p+1, label=N_genes_inter)) +
            xlab("Number of transcripts per gene") +
            ylab("Genes Percentage (%)") +
            cowplot::theme_cowplot() +
            theme(axis.title = element_text(size = rel(1.1)),
                  axis.text.x = element_text(size = rel(1.1), angle = 60, hjust = 1, vjust = 1),
                  axis.text.y = element_text(size = rel(1.1)),
                  legend.position = "none",
                  strip.text = element_text(size = rel(1.1)))
        }
        
        ggsave("results/2.pic/20.number_of_txs_for_each_gene.pdf", plot = p, width = 8, height = 6)
      }
      
      rm(txdb, tx, ggdat, p)
    }
    
    # get non-overlap tx for each gene
    if (T) {
      gene_body_ID <- qread(paste0("data/saved_data/11.tx_level_Body_region_length_ge3_bins.qs"), 
                            nthreads = 6)
      gene_body_ID$gene_id <- gene_tx_num_stat[match(gene_body_ID$tx_id, gene_tx_num_stat$txid), "geneid"]
      rm(gene_tx_num_stat)
      
      # 50272 genes
      length(unique(gene_body_ID$gene_id))
      non_overlap_txs <- lapply(split(gene_body_ID, gene_body_ID$gene_id), function(gene_dat) {
        # gene_dat <- gene_body_ID[gene_body_ID$gene_id == gene_body_ID$gene_id[12], ]
        tx_bins <- lapply(1:nrow(gene_dat), function(x) {
          gene_dat$tss_BinID[x] : gene_dat$tes_BinID[x]
        })
        
        gene_bins_ids <- do.call(c, tx_bins)
        bins_num_stat <- table(gene_bins_ids)
        over_bins <- as.numeric(names(bins_num_stat[bins_num_stat>1]))
        
        target_txs <- sapply(tx_bins, function(txid) {
          # txid <- tx_bins[[2]]
          tx_len <- length(txid)
          over_len <- sum(txid %in% over_bins)
          over_p <- over_len/tx_len
          
          if (over_p > over_cutoff) {
            return(F)
          } else {
            return(T)
          }
        })
        
        return(gene_dat$tx_id[target_txs])
      })
      
      # 232092 txs
      length(unique(gene_body_ID$tx_id))
      head(non_overlap_txs)
      
      # 0     1         2     3     6 
      # 22836 27396    36     3     1
      table(lengths(non_overlap_txs))
      # 27483
      sum(lengths(non_overlap_txs))
      non_overlap_txs <- do.call(c, non_overlap_txs)
      
      # 27483
      gene_body_ID <- gene_body_ID[gene_body_ID$tx_id %in% non_overlap_txs, ]
    }
    
    qsave(gene_body_ID, file, nthreads = 6)
  }
  if (file.exists(file)) {
    gene_body_ID <- qread(file, nthreads = 6)
  }
  
  table(table(gene_body_ID$gene_id))
}

# get mRNA level
if (T) {
  # read the data
  if (T) {
    dat <- qread("data/saved_data/7.median_rpkm_value_for_thp1_cd34.qs", nthreads = 6)
    dat$tx_id <- rownames(dat)
    
    # only retain txs not overlap with other txs
    dat <- dat[match(gene_body_ID$tx_id, dat$tx_id), ]
  }
  
  # Classification based on gene expression
  if (T) {
    head(dat)
    dat$rna_type_cell2 <- dat$rna_type_cell1 <- NA
    
    non0dat1 <- dat[, cell1][dat[, cell1]>0]
    non0dat2 <- dat[, cell2][dat[, cell2]>0]
    
    for (q in 1:rna_levels) {
      # q <- 1
      for (c in c(1, 2)) {
        # c <- 1
        if (c == 1) {
          non0dat <- non0dat1
        } else if (c == 2) {
          non0dat <- non0dat2
        }
        
        cutoff1 <- quantile(non0dat, (q-1)/10)
        cutoff2 <- quantile(non0dat, q/10)
        
        if (q == rna_levels) {
          dat[, paste0("rna_type_cell", c)] <- ifelse(dat[, c] == 0, "subdat0", 
                                                      ifelse(dat[, c]>=cutoff1 & dat[, c]<=cutoff2, 
                                                             paste0("Q", q), 
                                                             dat[, paste0("rna_type_cell", c)]))
        } else {
          dat[, paste0("rna_type_cell", c)] <- ifelse(dat[, c] == 0, "subdat0", 
                                                      ifelse(dat[, c]>=cutoff1 & dat[, c]<cutoff2, 
                                                             paste0("Q", q), 
                                                             dat[, paste0("rna_type_cell", c)]))
        }
      }
      
      rm(q, c, non0dat, cutoff1, cutoff2)
    }
    
    rm(non0dat1, non0dat2)
  }
  
  # format the factor
  if (T) {
    dat$rna_type_cell1 <- factor(dat$rna_type_cell1, levels = c("subdat0", paste0("Q", 1:rna_levels)))
    dat$rna_type_cell2 <- factor(dat$rna_type_cell2, levels = c("subdat0", paste0("Q", 1:rna_levels)))
  }
  
  # statistics
  if (T) {
    head(dat)
    table(dat$rna_type_cell1)
    table(dat$rna_type_cell2)
    
    head(dat)
    
    # cell1
    if (T) {
      non0min <- min(dat[, cell1][dat[, cell1]>0])
      head(dat)
      
      p <- ggplot(dat) +
        geom_boxplot(aes(x=rna_type_cell1, y=log(thp1+non0min), group=rna_type_cell1, fill=rna_type_cell1), outlier.shape=NA) +
        ggtitle(cell1) +
        labs(x = "mRNA Quantile Level", y = "log(rpkm)") +
        cowplot::theme_cowplot() +
        theme(axis.title = element_text(size = rel(1.2)),
              axis.text = element_text(size = rel(1.2)),
              legend.position = "none")
      
      file <- paste0("results/2.pic/20.non-overlap_tx_", rna_levels, "expression_level_", cell1, ".jpeg")
      ggsave(file, plot = p, width = 8, height = 6)
    }
    
    # cell2
    if (T) {
      non0min <- min(dat[, cell2][dat[, cell2]>0])
      head(dat)
      
      p <- ggplot(dat) +
        geom_boxplot(aes(x=rna_type_cell2, y=log(cd34+non0min), group=rna_type_cell2, fill=rna_type_cell2), outlier.shape=NA) +
        ggtitle(cell2) +
        labs(x = "mRNA Quantile Level", y = "log(rpkm)") +
        cowplot::theme_cowplot() +
        theme(axis.title = element_text(size = rel(1.2)),
              axis.text = element_text(size = rel(1.2)),
              legend.position = "none")
      
      file <- paste0("results/2.pic/20.non-overlap_tx_", rna_levels, "expression_level_", cell2, ".jpeg")
      ggsave(file, plot = p, width = 8, height = 6)
    }
    
    rm(p, file, non0min)
  }
}

# remove the end point
if (T) {
  # load the genome windows bin
  bin <- qread(windowID, nthreads = 6)
  
  bin$Chr <- seqnames(bin)
  ends <- do.call(rbind, tapply(bin$ID, bin$Chr, range))
  ends_bins <- c(ends)
  
  # the ends locate in the tx
  torf <- apply(gene_body_ID, 1, function(x) {
    # x <- unlist(gene_body_ID[1, ])
    
    genebody_start <- ifelse(x["strand"] == "+", as.numeric(x["tss_BinID"])+1, as.numeric(x["tss_BinID"])-1)
    genebody_end <- ifelse(x["strand"] == "+", as.numeric(x["tes_BinID"])-1, as.numeric(x["tes_BinID"])+1)
    torf <- sum(genebody_start:genebody_end %in% ends_bins)>0
    
    return(torf)
  })
  
  # the number of tx location in chromatin ends (0)
  sum(torf)
  mess <- paste0("There are ", sum(torf), " txs location in chromatin ends, these txs will be removed.\n")
  cat(mess)
  
  gene_body_ID <- gene_body_ID[!torf, ]
  
  rm(bin, ends, ends_bins, torf, mess)
}

# add mRNA info
if (T) {
  head(gene_body_ID)
  head(dat)
  
  gene_body_ID <- merge(gene_body_ID, dat, by="tx_id")
  gene_body_ID <- gene_body_ID[, c("tx_id", "tss_BinID", "tes_BinID", "strand", "Len", "gene_id", cell1, cell2, "rna_type_cell1", "rna_type_cell2")]
}

# change gene name 
if (T) {
  gene_body_ID$tx_id <- gsub("[_.]", "-", gene_body_ID$tx_id)
  gene_body_ID$gene_id <- gsub("[_.]", "-", gene_body_ID$gene_id)
}

file <- "data/saved_data/20.geneBody_metadata_tx_level.qs"
qsave(gene_body_ID, file, nthreads = 6)



# line1-202: 考虑到有些tx的tss或者tes可能会存在于blacklist中，因此先对所有的tx进行过滤，仅从与blacklist无关的tx中挑选
#            由于每个gene有多个tx，一个tx的tss上的state可能会被认为是另外一个tx的geneBody上的state，为了避免这种错误归属的问题
#            挑选tx时仅选择在同一个gene内，不同tx之间不存在overlap的tx进行后续分析
# line203-261: 找到thp1和cd34的以tx为单位的mRNA数据，并按照表达量，将所有基因分为subdat0以及Q1-Q10一共11级
# line262-288: 去除[TSS+1, TES-1]范围内存在染色质端点的gene
# line289-297: 合并gene tss以及tes信息 + 基因表达水平信息
# line298-302: 由于后续seurat中对于基因的识别中不能存在"-"或者"."，故这里直接将"[-.]"替换为"-"
