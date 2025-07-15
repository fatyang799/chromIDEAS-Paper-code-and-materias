# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(qs)
  suppressPackageStartupMessages(library(GenomicFeatures))
}

bin_size <- 200
up_bin_num <- 5
down_bin_num <- 5
body_bin_num <- 10
txdb_sqlite <- "data/raw_data/gencode.v40.annotation.sqlite"
windowID <- "data/saved_data/1.windowsNoBlack.withid.qs"

# load the genome windows bin
if (T) {
  bin <- qread(windowID, nthreads = 6)
}

for (data_type in c("gene", "tx")) {
  # data_type <- "tx"
  cat(paste0("Now process ", data_type, ":\n"))
  
  # construct the gtf annotation txdb object
  if (T) {
    txdb <- loadDb(txdb_sqlite)
    
    # only focus on the chr1-22.XY
    seqlevels(txdb) <- c(paste0("chr", 1:22), "chrX", "chrY")
    
    if (data_type == "tx") {
      tx <- transcripts(txdb, columns=c("TXNAME", "GENEID"))
    }
    if (data_type == "gene") {
      tx <- genes(txdb)
    }
  }
  
  # construct the gtf related genomic feature bin: TSS ± n_bins
  if (T) {
    cat(paste0("|", paste(rep("-", 50), collapse = ""), " TSS ", paste(rep("-", 50), collapse = ""), "|\n"))
    
    file <- paste0("data/saved_data/11.", data_type, "_level_profile_TSS_", up_bin_num, "UP_", down_bin_num, "DW.bs", bin_size, ".bin.qs")
    if (! file.exists(file)) {
      # tss regions
      if (T) {
        chr <- seqnames(tx)
        start <- start(tx)
        end <- end(tx)
        strand <- as.vector(strand(tx))
        
        tss <- GRanges(seqnames = chr, 
                       ranges = IRanges(start = ifelse(strand == "+", start, end),
                                        end = ifelse(strand == "+", start, end)))
        rm(chr, start, end)
      }
      
      # TSS bin ID
      if (T) {
        # find overlap
        overlapTSS <- findOverlaps(tss, bin, type="any")
        
        # summary the results
        if (data_type == "tx") {
          tssID <- data.frame(tx_id = tx$TXNAME[queryHits(overlapTSS)], 
                              Bin_ID = (bin$ID)[subjectHits(overlapTSS)], 
                              Strand = strand[queryHits(overlapTSS)])
        }
        if (data_type == "gene") {
          tssID <- data.frame(tx_id = tx$gene_id[queryHits(overlapTSS)], 
                              Bin_ID = (bin$ID)[subjectHits(overlapTSS)], 
                              Strand = strand[queryHits(overlapTSS)])
        }
        
        rm(strand, tss, overlapTSS)
      }
      
      # save the data
      if (T) {
        basic_file <- paste0("data/saved_data/11.", data_type, "_level_TSS_region_bins.qs")
        if (! file.exists(basic_file)) {
          qsave(tssID, file = basic_file, nthreads = 6)
        }
      }
      
      # TSS up down matrix
      if (T) {
        head(tssID)
        tss_mat <- data.frame(gene_id = tssID$tx_id)
        
        total_bin_num <- up_bin_num+1+down_bin_num
        
        for (bin_id in 1:total_bin_num) {
          name <- ifelse(bin_id<up_bin_num+1, paste0("U", up_bin_num-bin_id+1), 
                         ifelse(bin_id==up_bin_num+1, "TSS", paste0("D", bin_id-down_bin_num-1)))
          print(name)
          
          tss_mat[, name] <- ifelse(rep(bin_id<up_bin_num+1, nrow(tssID)), 
                                    (ifelse(tssID$Strand == "+", tssID$Bin_ID-(up_bin_num-bin_id+1), tssID$Bin_ID+(up_bin_num-bin_id+1))), 
                                    ifelse(rep(bin_id==up_bin_num+1, nrow(tssID)), 
                                           tssID$Bin_ID, 
                                           (ifelse(tssID$Strand == "+", tssID$Bin_ID+(bin_id-down_bin_num-1), tssID$Bin_ID-(bin_id-down_bin_num-1)))))
        }
      }
      
      # remove the end point
      if (T) {
        bin$Chr <- seqnames(bin)
        ends <- do.call(rbind, tapply(bin$ID, bin$Chr, range))
        ends_bins <- c(ends)
        
        # the ends locate in the tx
        torf <- apply(tss_mat, 1, function(x) {
          # x <- tss_mat[1, ]
          up <- as.numeric(x[paste0("U", up_bin_num-1)])
          down <- as.numeric(x[paste0("D", down_bin_num-1)])
          
          torf <- sum(up:down %in% ends_bins)>0
          
          return(torf)
        })
        
        # the number of tx location in chromatin ends (gene:1 tx:6)
        sum(torf)
        mess <- paste0("There are ", sum(torf), " ", data_type, " location in chromatin ends, these ", data_type, "s will be removed.\n")
        cat(mess)
        
        tss_mat <- tss_mat[!torf, ]
      }
      
      qsave(tss_mat, file = file, nthreads = 6)
      
      mess <- paste0("The result of [TSS-", up_bin_num, "bins, TSS+", down_bin_num, "bins] has been done.\n")
      cat(mess)
    }
    
    cat(paste0("|", paste(rep("-", 50), collapse = ""), " Done ", paste(rep("-", 50), collapse = ""), "|\n"))
  }
  
  # construct the gtf related genomic feature bin: Body
  if (T) {
    cat(paste0("|", paste(rep("-", 50), collapse = ""), " Body ", paste(rep("-", 50), collapse = ""), "|\n"))
    
    file <- paste0("data/saved_data/11.", data_type, "_level_profile_Body_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".bin.qs")
    if (! file.exists(file)) {
      # tx region
      if (T) {
        chr <- seqnames(tx)
        start <- start(tx)
        end <- end(tx)
        strand <- as.vector(strand(tx))
      }
      
      # TSS bin ID
      if (T) {
        tss <- GRanges(seqnames = chr, 
                       ranges = IRanges(start = ifelse(strand == "+", start, end),
                                        end = ifelse(strand == "+", start, end)))
        
        # find overlap
        overlapTSS <- findOverlaps(tss, bin, type="any")
        
        # summary the results
        if (data_type == "tx") {
          tssID <- data.frame(TSS_ID = tx$TXNAME[queryHits(overlapTSS)], 
                              Bin_ID = (bin$ID)[subjectHits(overlapTSS)], 
                              Strand = strand[queryHits(overlapTSS)])
        }
        if (data_type == "gene") {
          tssID <- data.frame(TSS_ID = tx$gene_id[queryHits(overlapTSS)], 
                              Bin_ID = (bin$ID)[subjectHits(overlapTSS)], 
                              Strand = strand[queryHits(overlapTSS)])
        }
        
        rm(overlapTSS, tss)
      }
      
      # TES bin ID
      if (T) {
        tes <- GRanges(seqnames = chr, 
                       ranges = IRanges(start = ifelse(strand == "+", end, start),
                                        end = ifelse(strand == "+", end, start)))
        
        # find overlap
        overlapTES <- findOverlaps(tes, bin, type="any")
        
        # summary the results
        if (data_type == "tx") {
          tesID <- data.frame(TES_ID = tx$TXNAME[queryHits(overlapTES)], 
                              Bin_ID = (bin$ID)[subjectHits(overlapTES)], 
                              Strand = strand[queryHits(overlapTES)])
        }
        if (data_type == "gene") {
          tesID <- data.frame(TES_ID = tx$gene_id[queryHits(overlapTES)], 
                              Bin_ID = (bin$ID)[subjectHits(overlapTES)], 
                              Strand = strand[queryHits(overlapTES)])
        }
        
        rm(overlapTES, tes, chr, start, end, strand)
      }
      
      # get tx with tes and tss at same time
      if (T) {
        over <- intersect(tesID$TES_ID, tssID$TSS_ID)
        tssID <- tssID[match(over, tssID$TSS_ID), ]
        tesID <- tesID[match(over, tesID$TES_ID), ]
        if (identical(tssID$Strand, tesID$Strand)) {
          gene_body_ID <- data.frame(tx_id = over, 
                                     tss_BinID = tssID$Bin_ID, 
                                     tes_BinID = tesID$Bin_ID, 
                                     strand = tssID$Strand)
          head(gene_body_ID)
        }
        
        gene_body_ID$Len <- ifelse(gene_body_ID$strand == "+", 
                                   gene_body_ID$tes_BinID - gene_body_ID$tss_BinID + 1, 
                                   gene_body_ID$tss_BinID - gene_body_ID$tes_BinID + 1)
        rm(over, tssID, tesID)
      }
      
      # filter txs: length >= 3
      if (T) {
        torf <- gene_body_ID$Len >= 3
        # gene:50186 tx:232092
        sum(torf)
        
        cat(paste0("The minimum ", data_type, " length should be more than 3, including at least tss, tes, and a gene body bin.\n"))
        mess <- paste0("Filter out ", sum(! torf), "/", nrow(gene_body_ID), " (", round(sum(! torf)/nrow(gene_body_ID)*100, 2), 
                       "%) ", data_type, "s, whose length is less than 3 bins.\n")
        cat(mess)
        
        gene_body_ID <- gene_body_ID[torf, ]
        
        rm(torf, mess)
      }
      
      # save the data
      if (T) {
        basic_file <- paste0("data/saved_data/11.", data_type, "_level_Body_region_length_ge3_bins.qs")
        if (! file.exists(basic_file)) {
          qsave(gene_body_ID, file = basic_file, nthreads = 6)
        }
      }
      
      # manual
      if (F) {
        # detailed logic see: E:\OneDrive\2021级 杨柳 人类染色质状态的特征\Yangliu with Yen\Lab meeting\03 小组会\20240815-纠正GTF区分正负链+pol2数据质控+基因体等分算法纠正并明确算法思路.pptx
        ## name:
        ##     bin_id<up_bin_num+1: paste0("U", up_bin_num+1-bin_id)
        ##     bin_id==up_bin_num+1: "TSS"
        ##     bin_id<up_bin_num+body_bin_num+2: paste0("B", bin_id-(up_bin_num+1))
        ##     bin_id==up_bin_num+body_bin_num+2: "TES"
        ##     else:
        ##         paste0("D", bin_id-(up_bin_num+body_bin_num+2))
        ## 
        ## gene_body_mat[, name]:
        ##     rep(bin_id<up_bin_num+1, nrow(gene_body_mat)):
        ##         "+": gene_body_ID$tss_BinID-(up_bin_num+1-bin_id)
        ##         "-": gene_body_ID$tss_BinID+(up_bin_num+1-bin_id)
        ##     rep(bin_id==up_bin_num+1, nrow(gene_body_mat)):
        ##         gene_body_ID$tss_BinID
        ##     rep(bin_id<up_bin_num+body_bin_num+2, nrow(gene_body_mat)):
        ##         "+": 
        ##             gene_body_mat$unit>=1: paste0(
        ##                 pos_strand_rounding((gene_body_ID$tss_BinID+1) + (bin_id-(up_bin_num+1) -1) * gene_body_mat$unit), 
        ##                 "-", 
        ##                 pos_strand_rounding((gene_body_ID$tss_BinID+1) + (bin_id-(up_bin_num+1)) * gene_body_mat$unit - 1)
        ##                 )
        ##             gene_body_mat$unit<1: paste0(
        ##                 floor((gene_body_ID$tss_BinID+1) + (bin_id-(up_bin_num+1) -1) * gene_body_mat$unit + 1e-5), 
        ##                 "-", 
        ##                 floor((gene_body_ID$tss_BinID+1) + (bin_id-(up_bin_num+1)) * gene_body_mat$unit - 1e-5)
        ##                 )
        ##         "-": 
        ##             gene_body_mat$unit>=1: paste0(
        ##                 neg_strand_rounding((gene_body_ID$tes_BinID+1) + (body_bin_num+1-(bin_id-(up_bin_num+1)) -1) * gene_body_mat$unit + 1), 
        ##                 "-", 
        ##                 neg_strand_rounding((gene_body_ID$tes_BinID+1) + (body_bin_num+1-(bin_id-(up_bin_num+1))) * gene_body_mat$unit)
        ##                 )
        ##             gene_body_mat$unit<1: paste0(
        ##                 floor((gene_body_ID$tes_BinID+1) + (body_bin_num+1-(bin_id-(up_bin_num+1)) -1) * gene_body_mat$unit + 1e-5), 
        ##                 "-", 
        ##                 floor((gene_body_ID$tes_BinID+1) + (body_bin_num+1-(bin_id-(up_bin_num+1))) * gene_body_mat$unit - 1e-5)
        ##                 )
        ##     rep(bin_id==up_bin_num+body_bin_num+2, nrow(gene_body_mat)):
        ##         gene_body_ID$tes_BinID
        ##     else
        ##         "+": gene_body_ID$tes_BinID+(bin_id-(up_bin_num+body_bin_num+2))
        ##         "-": gene_body_ID$tes_BinID-(bin_id-(up_bin_num+body_bin_num+2))
      }
      
      # Tx up down matrix
      if (T) {
        total_bin_num <- up_bin_num+body_bin_num+down_bin_num
        
        head(gene_body_ID)
        gene_body_mat <- data.frame(gene_id = gene_body_ID$tx_id)
        gene_body_mat$strand <- gene_body_ID$strand
        gene_body_mat$unit <- ifelse(gene_body_ID$strand == "+", 
                                     ((gene_body_ID$tes_BinID-1) - (gene_body_ID$tss_BinID+1) +1) / body_bin_num, 
                                     ((gene_body_ID$tss_BinID-1) - (gene_body_ID$tes_BinID+1) +1) / body_bin_num)
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
        
        for (bin_id in 1:(total_bin_num+2)) {
          name <- ifelse(bin_id<up_bin_num+1, paste0("U", up_bin_num+1-bin_id), 
                         ifelse(bin_id==up_bin_num+1, "TSS", 
                                ifelse(bin_id<up_bin_num+body_bin_num+2, paste0("B", bin_id-(up_bin_num+1)), 
                                       ifelse(bin_id==up_bin_num+body_bin_num+2, "TES", paste0("D", bin_id-(up_bin_num+body_bin_num+2))))))
          print(name)
          
          gene_body_mat[, name] <- ifelse(rep(bin_id<up_bin_num+1, nrow(gene_body_mat)), 
                                          ifelse(gene_body_ID$strand == "+", 
                                                 gene_body_ID$tss_BinID-(up_bin_num+1-bin_id), 
                                                 gene_body_ID$tss_BinID+(up_bin_num+1-bin_id)), 
                                          
                                          ifelse(rep(bin_id==up_bin_num+1, nrow(gene_body_mat)), 
                                                 gene_body_ID$tss_BinID, 
                                                 
                                                 ifelse(rep(bin_id<up_bin_num+body_bin_num+2, nrow(gene_body_mat)), 
                                                        ifelse(gene_body_ID$strand == "+", 
                                                               ifelse(gene_body_mat$unit>=1, 
                                                                      paste0(
                                                                        pos_strand_rounding((gene_body_ID$tss_BinID+1) + (bin_id-(up_bin_num+1) -1) * gene_body_mat$unit), 
                                                                        "-", 
                                                                        pos_strand_rounding((gene_body_ID$tss_BinID+1) + (bin_id-(up_bin_num+1)) * gene_body_mat$unit - 1)
                                                                      ), 
                                                                      paste0(
                                                                        floor((gene_body_ID$tss_BinID+1) + (bin_id-(up_bin_num+1) -1) * gene_body_mat$unit + 1e-5), 
                                                                        "-", 
                                                                        floor((gene_body_ID$tss_BinID+1) + (bin_id-(up_bin_num+1)) * gene_body_mat$unit - 1e-5)
                                                                      )), 
                                                               ifelse(gene_body_mat$unit>=1, 
                                                                      paste0(
                                                                        neg_strand_rounding((gene_body_ID$tes_BinID+1) + (body_bin_num+1-(bin_id-(up_bin_num+1)) -1) * gene_body_mat$unit + 1), 
                                                                        "-", 
                                                                        neg_strand_rounding((gene_body_ID$tes_BinID+1) + (body_bin_num+1-(bin_id-(up_bin_num+1))) * gene_body_mat$unit)
                                                                      ), 
                                                                      paste0(
                                                                        floor((gene_body_ID$tes_BinID+1) + (body_bin_num+1-(bin_id-(up_bin_num+1)) -1) * gene_body_mat$unit + 1e-5), 
                                                                        "-", 
                                                                        floor((gene_body_ID$tes_BinID+1) + (body_bin_num+1-(bin_id-(up_bin_num+1))) * gene_body_mat$unit - 1e-5)
                                                                      ))),
                                                        
                                                        ifelse(rep(bin_id==up_bin_num+body_bin_num+2, nrow(gene_body_mat)), 
                                                               gene_body_ID$tes_BinID,
                                                               
                                                               ifelse(gene_body_ID$strand == "+", 
                                                                      gene_body_ID$tes_BinID+(bin_id-(up_bin_num+body_bin_num+2)), 
                                                                      gene_body_ID$tes_BinID-(bin_id-(up_bin_num+body_bin_num+2)))))))
        }
        
        rm(total_bin_num, bin_id, name, pos_strand_rounding, neg_strand_rounding)
      }
      
      # remove the end point
      if (T) {
        bin$Chr <- seqnames(bin)
        ends <- do.call(rbind, tapply(bin$ID, bin$Chr, range))
        ends_bins <- c(ends)
        
        # the ends locate in the tx
        torf <- apply(gene_body_mat, 1, function(x) {
          # x <- unlist(gene_body_mat[1, ])
          up <- as.numeric(x[paste0("U", up_bin_num-1)])
          down <- as.numeric(x[paste0("D", down_bin_num-1)])
          
          torf <- sum(up:down %in% ends_bins)>0
          
          return(torf)
        })
        
        # the number of tx location in chromatin ends (gene:2 tx:6)
        sum(torf)
        mess <- paste0("There are ", sum(torf), " ", data_type, " location in chromatin ends, these ", data_type, "s will be removed.\n")
        cat(mess)
        
        gene_body_mat <- gene_body_mat[!torf, ]
      }
      
      qsave(gene_body_mat, file = file, nthreads = 6)
      
      mess <- paste0("The result of ", up_bin_num, "bins-[TSS, ", body_bin_num, "bins_Body , TES]-", down_bin_num, "bins has been done.\n")
      cat(mess)
    }
    
    cat(paste0("|", paste(rep("-", 50), collapse = ""), " Done ", paste(rep("-", 50), collapse = ""), "|\n"))
  }
}



# line1-20: 读入binID的坐标信息

# line21-274: 以转录本为单位进行区域分割
#   1) line22-32: 提取chr1-22,XY上的所有转录本
#   2) line33-107: 针对所有转录本而言，以TSS为ref点，根据指定上下游区域，获取相应binID，并去除处于染色质端点的转录本
#   3) line108-274: 
#               1. 针对所有转录本而言，先进行过滤，需要保证同时满足以下2点方可进行后续分析:
#                   - TSS和TES均不与blacklist有关
#                   - 转录本长度>3，从而满足TSS, TES, 以及gene body (该区域可以进一步在body_bin_num的设定下进行分割)至少可以分配1个bin进行后续分析
#               2. 对于过滤后的转录本而言，以TSS和TES为ref点:
#                   - bin_id<up_bin_num+1（TSS上游区域）: gene_body_ID$tss_BinID-(up_bin_num+1-bin_id)
#                   - bin_id==up_bin_num+1（TSS ref点）: gene_body_ID$tss_BinID
#                   - bin_id<up_bin_num+body_bin_num+2（TSS+1 ~ TES-1）: 基因body区域等分为特定个区块
#                     起点：【(gene_body_ID$tss_BinID+1) + ((gene_body_ID$tes_BinID-1)-(gene_body_ID$tss_BinID+1)) / body_bin_num * (bin_id-(up_bin_num+1) - 1)】
#                     终点：【(gene_body_ID$tss_BinID+1) + ((gene_body_ID$tes_BinID-1)-(gene_body_ID$tss_BinID+1)) / body_bin_num * (bin_id-(up_bin_num+1))】
#                   - bin_id==up_bin_num+body_bin_num+2（TES ref点）: gene_body_ID$tes_BinID
#                   - 其余区域（TSS下游区域）: gene_body_ID$tes_BinID+(bin_id-(up_bin_num+body_bin_num+2))
#                   - **注意**：为避免重复计算问题，仅最后一个区间为首尾均闭合，其他区间均为【）模式
#               3. 去除处于染色质端点的转录本

# line278-531: 以基因为单位进行区域分割
#   1) line279-289: 提取chr1-22,XY上的所有基因
#   2) line290-364: 针对所有基因而言，以TSS为ref点，根据指定上下游区域，获取相应binID，并去除处于染色质端点的基因
#   3) line365-531: 
#               1. 针对所有基因而言，先进行过滤，需要保证同时满足以下2点方可进行后续分析:
#                   - TSS和TES均不与blacklist有关
#                   - 基因长度>3，从而满足TSS, TES, 以及gene body (该区域可以进一步在body_bin_num的设定下进行分割)至少可以分配1个bin进行后续分析
#               2. 对于过滤后的基因而言，以TSS和TES为ref点:
#                   - bin_id<up_bin_num+1（TSS上游区域）: gene_body_ID$tss_BinID-(up_bin_num+1-bin_id)
#                   - bin_id==up_bin_num+1（TSS ref点）: gene_body_ID$tss_BinID
#                   - bin_id<up_bin_num+body_bin_num+2（TSS+1 ~ TES-1）: 基因body区域等分为特定个区块
#                     起点：【(gene_body_ID$tss_BinID+1) + ((gene_body_ID$tes_BinID-1)-(gene_body_ID$tss_BinID+1)) / body_bin_num * (bin_id-(up_bin_num+1) - 1)】
#                     终点：【(gene_body_ID$tss_BinID+1) + ((gene_body_ID$tes_BinID-1)-(gene_body_ID$tss_BinID+1)) / body_bin_num * (bin_id-(up_bin_num+1))】
#                   - bin_id==up_bin_num+body_bin_num+2（TES ref点）: gene_body_ID$tes_BinID
#                   - 其余区域（TSS下游区域）: gene_body_ID$tes_BinID+(bin_id-(up_bin_num+body_bin_num+2))
#                   - **注意**：为避免重复计算问题，仅最后一个区间为首位均闭合，其他区间均为【）模式
#               3. 去除处于染色质端点的基因