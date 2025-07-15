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
mode <- "mode1"

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
    file <- paste0("data/saved_data/26.stratified_random_sampling.", mode, ".qs")
    if (file.exists(file)) {
      target_tx <- qread(file, nthreads = 6)
    }
    if (! file.exists(file)) {
      # 每个组别中选20个tx进行后续分析，相同组别的tx按sd最大+表达量最高进行挑选
      target_tx_type <- paste0("Q", 1:length_leveles)
      target_tx <- lapply(target_tx_type, function(type) {
        # type <- target_tx_type[2]
        # mess
        if (T) {
          n <- which(target_tx_type == type)
          total <- length(target_tx_type)
          
          cat(paste0("\t(", n, "/", total, ") ", type, "\n"))
        }
        
        # get gene body part chromatin state percentage mat
        if (T) {
          gene_body_mat1 <- qread(paste0("data/saved_data/21.tx_Body(", body_bin_num, "parts)_percentage_based_on_Body(", length_leveles, 
                                         "parts)_RNA(", rna_levels, "parts).bs", bin_size, ".", cell1, ".qs"), nthreads = 6)
          gene_body_mat2 <- qread(paste0("data/saved_data/21.tx_Body(", body_bin_num, "parts)_percentage_based_on_Body(", length_leveles, 
                                         "parts)_RNA(", rna_levels, "parts).bs", bin_size, ".", cell2, ".qs"), nthreads = 6)
        }
        
        # get consistent id
        if (T) {
          common_dat <- function(gene_body_mat, gene_body_ID) {
            rownames_mat <- rownames(gene_body_mat)
            rownames_mat <- str_split(rownames_mat, "@", simplify = T)[, 1]
            id_mat <- unique(rownames_mat)
            
            overid <- intersect(id_mat, gene_body_ID$tx_id)
            
            gene_body_ID <- gene_body_ID[match(overid, gene_body_ID$tx_id), ]
            gene_body_mat <- gene_body_mat[rownames_mat %in% overid, ]
            
            return(gene_body_mat)
          }
          gene_body_mat1 <- common_dat(gene_body_mat1, gene_body_ID)
          gene_body_mat2 <- common_dat(gene_body_mat2, gene_body_ID)
        }
        
        # common value
        if (identical(rownames(gene_body_mat1), rownames(gene_body_mat2)) & identical(colnames(gene_body_mat1), colnames(gene_body_mat2))) {
          all_states <- colnames(gene_body_mat1)
          all_states_n <- as.numeric(gsub("S", "", all_states))
          
          # default select 2000 genes
          total <- 2000
          n_length_state <- ceiling(
            (total/length(all_states)) / length_leveles
          )
        } else {
          stop("The order of cell1 data and cell2 data are not identical")
        }
        
        # filter data
        if (T) {
          subid <- gene_body_ID$tx_id[gene_body_ID$Len_type == type]
          submat1 <- gene_body_mat1[str_split(rownames(gene_body_mat1), "@", simplify = T)[, 1] %in% subid, ]
          submat2 <- gene_body_mat2[str_split(rownames(gene_body_mat2), "@", simplify = T)[, 1] %in% subid, ]
        }
        
        # get statistics value for all genes
        if (T) {
          stat_state_percentage <- function(subid, submat) {
            res <- data.frame(t(
              sapply(subid, function(id) {
                # id <- subid[1]
                
                id_mat <- submat[str_split(rownames(submat), "@", simplify = T)[, 1] == id, ]
                sapply(id_mat, sum)
              })
            ))
            res$geneid <- rownames(res)
            rownames(res) <- 1:nrow(res)
            
            return(res)
          }
          stat_ids1 <- stat_state_percentage(subid, submat1)
          stat_ids2 <- stat_state_percentage(subid, submat2)
        }
        
        # for each state, get genes with most specific state percentage
        if (T) {
          target_gene <- lapply(all_states, function(state) {
            # state <- all_states[1]
            
            # get state specific genes function
            get_cell_specific_genes <- function(stat_ids, state, n_length_state) {
              stat_ids <- stat_ids[order(stat_ids[, state], decreasing = T), ]
              state_specific <- stat_ids$geneid
              
              res <- head(state_specific, n_length_state)
              
              return(res)
            }
            
            res1 <- get_cell_specific_genes(stat_ids1, state, n_length_state)
            res2 <- get_cell_specific_genes(stat_ids2, state, n_length_state)
            
            res <- data.frame(cell1 = res1, 
                              cell2 = res2, 
                              len = type, 
                              state = state)
            
            return(res)
          })
        }
        
        # merge the res
        if (T) {
          target_gene <- data.frame(do.call(rbind, target_gene))
        }
        
        return(target_gene)
      })
      
      # merge the data
      if (T) {
        target_tx <- data.frame(do.call(rbind, target_tx))
      }
      
      qsave(target_tx, file = file, nthreads = 6)
    }
  }
}

# get sampled txs
if (T) {
  head(target_tx)
  target_tx <- melt(target_tx, id.vars = c("len", "state"), measure.vars = c("cell1", "cell2"), variable.name = "cell", value.name = "target_genes")
  
  # 2623
  length(unique(target_tx$target_genes))
  table(target_tx$cell, target_tx$state)
  
  # 2623
  gene_body_ID <- gene_body_ID[gene_body_ID$tx %in% target_tx$target_genes, ]
  
  # add additional info about cell
  gene_body_ID$cell1_info <- ifelse(gene_body_ID$tx_id %in% target_tx$target_genes[target_tx$cell == "cell1"], T, F)
  gene_body_ID$cell2_info <- ifelse(gene_body_ID$tx_id %in% target_tx$target_genes[target_tx$cell == "cell2"], T, F)
  
  # 1780
  length(unique(target_tx$target_genes[target_tx$cell == "cell1"]))
  table(gene_body_ID$cell1_info)
  # 1777
  length(unique(target_tx$target_genes[target_tx$cell == "cell2"]))
  table(gene_body_ID$cell2_info)
  
  sort(table(gene_body_ID$rna_type_cell1)/nrow(gene_body_ID)*100)
  sort(table(gene_body_ID$rna_type_cell2)/nrow(gene_body_ID)*100)
  
  table(gene_body_ID$rna_type_cell1, gene_body_ID$Len_type)
  table(gene_body_ID$rna_type_cell2, gene_body_ID$Len_type)
}

# save the data
if (T) {
  file <- paste0("data/saved_data/26.geneBody_metadata_tx_level_stratified_random_sampling.", mode, ".qs")
  qsave(gene_body_ID, file = file, nthreads = 6)
}

# 合并后的抽样原则
# 分别对2种细胞数据进行挑选，按照state比例由高到低进行挑选，每种state在2种细胞中分别挑选n_length_state个基因，相当于对单独挑选的结果进行直接合并mode1：
#     优势：更凸显每种染色质状态的特点，分群效果更好。且没有细胞bias

# 针对每个基因，计算每种染色质状态在2种细胞中的平均百分比，按照state比例由高到低进行挑选，每种state在2种细胞中一共挑选n_length_state个基因mode2：
#     优势：没有细胞bias
#     劣势：部分染色质状态特征可能会被平均

# 合并2种细胞数据挑选时，直接合并数据，在合并后的数据中按照state比例进行挑选，每种state在2种细胞中一共挑选n_length_state个基因mode3：
#     优势：更凸显每种染色质状态的特点，分群效果更好
#     劣势：可能某种细胞数据会占主体，需检查细胞比例bias

