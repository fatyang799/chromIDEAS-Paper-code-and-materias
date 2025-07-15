# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(reshape2)
  library(stringr)
  library(qs)
}

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
body_bin_num <- 10
rna_levels <- 10
length_leveles <- 10

# get gene body metadata: gene_body_ID
if (T) {
  file <- "data/saved_data/20.geneBody_metadata_tx_level.qs"
  gene_body_ID <- qread(file)
}

# length filter > body_bin_num
if (T) {
  # 12846
  mess <- paste0("There are ", sum(gene_body_ID$Len >= body_bin_num), " (", 
                 round(sum(gene_body_ID$Len >= body_bin_num)/nrow(gene_body_ID)*100, 2), 
                 "%) genes passing the body length filter.\n")
  cat(mess)
  gene_body_ID <- gene_body_ID[gene_body_ID$Len >= body_bin_num, ]
  
  rm(mess)
}

# gene body length statistics
if (T) {
  head(gene_body_ID)
  
  gene_body_ID$Len_type <- NA
  for (q in 1:length_leveles) {
    # q <- 1
    cutoff1 <- quantile(gene_body_ID$Len, (q-1)/10)
    cutoff2 <- quantile(gene_body_ID$Len, q/10)
    
    if (q<length_leveles) {
      gene_body_ID$Len_type <- ifelse(gene_body_ID$Len>=cutoff1 & gene_body_ID$Len<cutoff2, 
                                      paste0("Q", q), 
                                      gene_body_ID$Len_type)
    } else {
      gene_body_ID$Len_type <- ifelse(gene_body_ID$Len>=cutoff1 & gene_body_ID$Len<=cutoff2, 
                                      paste0("Q", q), 
                                      gene_body_ID$Len_type)
    }
    
  }
  
  rm(cutoff1, cutoff2, q)
}

# stratified random sampling based on the length and the expression level: target_tx
if (T) {
  # statistics cell 1
  if (T) {
    tab <- as.data.frame.table(table(paste0(gene_body_ID$rna_type_cell1, "_", gene_body_ID$Len_type)))
    tab$rna <- gsub("_Q[0-9]{1,3}", "", tab$Var1)
    tab$rna <- factor(tab$rna, levels = c("subdat0", paste0("Q", seq(1, body_bin_num))))
    tab$len <- stringr::str_extract(tab$Var1, "Q[0-9]{1,3}$")
    tab$len <- factor(tab$len, levels = paste0("Q", seq(1, body_bin_num)))
    
    ggplot(data = tab) +
      geom_bar(aes(x=len, y=Freq, fill=rna), stat = "identity") +
      facet_grid(rna  ~ ., scales = "fixed") + 
      ggtitle(cell1) +
      theme(legend.position = "none")
  }
  
  # statistics cell 2
  if (T) {
    tab <- as.data.frame.table(table(paste0(gene_body_ID$rna_type_cell2, "_", gene_body_ID$Len_type)))
    tab$rna <- gsub("_Q[0-9]{1,3}", "", tab$Var1)
    tab$rna <- factor(tab$rna, levels = c("subdat0", paste0("Q", seq(1, body_bin_num))))
    tab$len <- stringr::str_extract(tab$Var1, "Q[0-9]{1,3}$")
    tab$len <- factor(tab$len, levels = paste0("Q", seq(1, body_bin_num)))
    
    ggplot(data = tab) +
      geom_bar(aes(x=len, y=Freq, fill=rna), stat = "identity") +
      facet_grid(rna  ~ ., scales = "fixed") + 
      ggtitle(cell2) +
      theme(legend.position = "none")
  }
  rm(tab)
  
  # stratified random sampling
  if (T) {
    file <- "data/saved_data/22.stratified_random_sampling.qs"
    if (file.exists(file)) {
      target_tx <- qread(file, nthreads = 6)
    }
    if (! file.exists(file)) {
      # 每个组别中选20个tx进行后续分析，相同组别的tx按sd最大+表达量最高进行挑选
      target_tx <- lapply(seq(1:2), function(c) {
        # c <- 1
        cell <- ifelse(c == 1, cell1, cell2)
        cat(paste0("Now process ", cell, " dat:\n"))
        
        # get gene body part chromatin state percentage mat
        if (T) {
          gene_body_mat <- qread(paste0("data/saved_data/21.tx_Body(", body_bin_num, "parts)_percentage_based_on_Body(", length_leveles, 
                                        "parts)_RNA(", rna_levels, "parts).bs", bin_size, ".", cell, ".qs"), nthreads = 6)
        }
        
        # get consistent id
        if (T) {
          rownames_mat <- rownames(gene_body_mat)
          rownames_mat <- str_split(rownames_mat, "@", simplify = T)[, 1]
          id_mat <- unique(rownames_mat)
          
          overid <- intersect(id_mat, gene_body_ID$tx_id)
          
          gene_body_ID <- gene_body_ID[match(overid, gene_body_ID$tx_id), ]
          gene_body_mat <- gene_body_mat[rownames_mat %in% overid, ]
          
          rm(overid, id_mat, rownames_mat)
        }
        
        target_tx_type <- paste0("Q", 1:length_leveles)
        
        # common value
        if (T) {
          all_states <- colnames(gene_body_mat)
          all_states_n <- as.numeric(gsub("S", "", all_states))
          
          # default select 2000 genes
          total <- 2000
          n_length_state <- ceiling(
            (total/length(all_states)) / length_leveles
          )
        }
        
        # get group specific target txs
        if (T) {
          target_tx <- lapply(target_tx_type, function(type) {
            # type <- target_tx_type[1]
            # mess
            if (T) {
              n <- which(target_tx_type == type)
              total <- length(target_tx_type)
              
              cat(paste0("\t(", n, "/", total, ") ", type, "\n"))
            }
            
            # filter data
            if (T) {
              subid <- gene_body_ID$tx_id[gene_body_ID$Len_type == type]
              submat <- gene_body_mat[str_split(rownames(gene_body_mat), "@", simplify = T)[, 1] %in% subid, ]
            }
            
            # get statistics value for all genes
            if (T) {
              stat_ids <- data.frame(t(
                sapply(subid, function(id) {
                  # id <- subid[1]
                  
                  id_mat <- submat[str_split(rownames(submat), "@", simplify = T)[, 1] == id, ]
                  sapply(id_mat, sum)
                })
              ))
              stat_ids$geneid <- rownames(stat_ids)
            }
            
            # for each state, get genes with most specific state percentage
            if (T) {
              target_gene <- sapply(all_states, function(state) {
                # state <- all_states[1]
                
                # get state specific genes
                if (T) {
                  stat_ids <- stat_ids[order(stat_ids[, state], decreasing = T), ]
                  state_specific <- stat_ids$geneid
                }
                
                res <- head(state_specific, n_length_state)
                
                return(res)
              })
              
              if (n_length_state>1) {
                target_gene <- as.data.frame(target_gene)
              } else if (n_length_state == 1) {
                target_gene <- data.frame(matrix(target_gene, nrow = 1, dimnames = list("1", names(target_gene))))
              }
            }
            
            # format the res
            if (T) {
              target_gene <- melt(target_gene, id.vars = NULL, measure.vars = all_states, variable.name = "State", value.name = "target_genes")
              target_gene$Cell <- cell
              target_gene$Len <- type
            }
            
            return(target_gene)
          })
        }
        
        # format the data
        if (T) {
          target_tx <- data.frame(do.call(rbind, target_tx))
        }
        
        return(target_tx)
      })
      target_tx <- data.frame(do.call(rbind, target_tx))
      
      qsave(target_tx, file = file, nthreads = 6)
    }
  }
}

# get sampled txs
if (T) {
  # 2623
  length(unique(target_tx$target_genes))
  
  # 2623
  gene_body_ID <- gene_body_ID[gene_body_ID$tx %in% target_tx$target_genes, ]
  sort(table(gene_body_ID$rna_type_cell1)/nrow(gene_body_ID)*100)
  sort(table(gene_body_ID$rna_type_cell2)/nrow(gene_body_ID)*100)
  
  table(gene_body_ID$rna_type_cell1, gene_body_ID$Len_type)
  table(gene_body_ID$rna_type_cell2, gene_body_ID$Len_type)
}

# info merge and summary
if (T) {
  head(gene_body_ID)
  head(target_tx)
  
  # get cell info
  gene_body_ID$cell1 <- cell1
  gene_body_ID$cell2 <- cell2
  
  # for cell1
  genes <- target_tx[target_tx$Cell == cell1, ]
  gene_body_ID[match(genes$target_genes, gene_body_ID$tx), "cell1_info"] <- genes$Len
  
  # for cell2
  genes <- target_tx[target_tx$Cell == cell2, ]
  gene_body_ID[match(genes$target_genes, gene_body_ID$tx), "cell2_info"] <- genes$Len
}

# save the data
if (T) {
  file <- "data/saved_data/22.geneBody_metadata_tx_level_stratified_random_sampling.qs"
  qsave(gene_body_ID, file = file, nthreads = 6)
}

# line1-22: 读入数据，gene_body_ID信息，里面包含了tss、tes的信息以及mRNA的信息
# line23-34: 对gene_body_ID数据进行过滤：基因所占bin数一定要>=body_bin_num，以保证每个bin至少独立享受1个bin
#               - 过滤：后续分层抽样时从这里进行抽样
#               - 未过滤：当所有chromatin states分完群后，进行差异分析时使用该矩阵
# line35-50: 在过滤完成后，统计剩余基因的基因body长度，并进行分级。
# line51-115: 这里需要根据基因长度以及基因表达量进行**分层抽样**。**分层抽样**相关考虑：
#              - 基因表达量：所有基因分为subdat0以及Q1-Q10一共11级
#              - 基因body长度：所有基因长度分为Q1-Q10一共10级
#              - 每层抽取基因数量：单细胞中2000个高变基因即可完成细胞准确分群，故这里考虑每层一共抽取20个基因，这样分层后共会抽取20*10*11=2200个基因
#              - 每层抽样选择：在相同表达水平、相同基因长度程度中，挑选细胞中表达量最高的20个基因。如果数量不够20，则有多少取多少
# line116-122: 在gene_body_ID中仅保留抽样得到的基因信息
# line123-140: 在gene_body_ID变量中加入RNAseq的信息


