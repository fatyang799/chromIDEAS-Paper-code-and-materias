# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(reshape2)
  library(stringr)
  library(qs)
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(circlize))
  suppressPackageStartupMessages(library(ComplexHeatmap))
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
  gene_body_ID <- gene_body_ID[gene_body_ID$Len >= body_bin_num, ]
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

total <- c(seq(100, 400, 100), 
           seq(500, 4000, 500))

# stratified random sampling
if (T) {
  file <- "data/saved_data/34.stratified_random_sampling_various_number_highly_informative_txs.qs"
  if (file.exists(file)) {
    target_tx <- qread(file, nthreads = 6)
  }
  if (! file.exists(file)) {
    # get sampling selection
    if (T) {
      target_tx_max <- lapply(seq(1:2), function(c) {
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
          
          n_length_state <- ceiling(
            (max(total)/length(all_states)) / length_leveles
          )
          n_length_state <- 20
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
              target_gene$Cell <- cell
              target_gene$Len <- type
              target_gene$Rank <- 1:n_length_state
            }
            
            return(target_gene)
          })
        }
        
        # merge the res
        if (T) {
          target_tx <- data.frame(do.call(rbind, target_tx))
        }
        
        return(target_tx)
      })
    }
    
    # format the data
    if (T) {
      target_tx <- data.frame(do.call(rbind, target_tx_max))
    }
    
    qsave(target_tx, file = file, nthreads = 6)
  }
  
  rm(gene_body_ID)
}

resolutions <- c(seq(0.1, 0.9, 0.1), 
                 seq(1, 1.8, 0.2), 
                 seq(2, 5, 1))

# clustering res
if (T) {
  file <- "data/saved_data/34.various_highly_informative_txs_mat.qs"
  if (file.exists(file)) {
    hit_dat <- qread(file, nthreads = 6)
  }
  if (! file.exists(file)) {
    res <- lapply(c(cell1, cell2), function(cell) {
      # cell <- "thp1"
      # cell <- "cd34"
      
      # info
      if (T) {
        cat(paste0("Now process ", cell, ":\n"))
      }
      
      # read the raw data
      if (T) {
        # emission data
        emission <- read.table("data/raw_data/1.emission/chromIDEAS.emission.txt", header = T, sep = "\t")
      }
      
      subres <- lapply(total, function(sub_total) {
        # sub_total <- total[1]
        
        # info
        if (T) {
          cat(paste0("\t", sub_total, " highly informative txs\n"))
        }
        
        # gene body state percentage
        if (T) {
          n_length_state <- ceiling(
            (sub_total/nrow(emission)) / length_leveles
          )
          cell_gene_pat <- target_tx[target_tx$Cell == cell, ]
          target_gene_pat <- cell_gene_pat[cell_gene_pat$Rank <= n_length_state, emission$State]
          target_gene_pat <- unique(c(as.matrix(target_gene_pat)))
          
          set.seed(799)
          remain <- sub_total - length(target_gene_pat)
          remain_pool <- unique(c(as.matrix(cell_gene_pat[cell_gene_pat$Rank > n_length_state, emission$State])))
          if (remain>0) {
            target_gene_pat <- c(target_gene_pat, 
                                 sample(remain_pool, remain, replace = F))
          } else {
            target_gene_pat <- sample(target_gene_pat, sub_total, replace = F)
          }
          
          # CS percentage input matrix
          gene <- qread(paste0("data/saved_data/21.tx_Body(", body_bin_num, "parts)_percentage_based_on_Body(", 
                               length_leveles, "parts)_RNA(", rna_levels, "parts).bs", bin_size, ".", cell, ".qs"))
          rowname <- rownames(gene)
          rowname <- str_split(rowname, "@", simplify = T)[, 1]
          torf <- rowname %in% target_gene_pat
          sum(torf)
          gene <- gene[torf, ]
          
          rm(target_gene_pat, n_length_state, rowname, torf)
        }
        
        # Create Seurat Object
        if (T) {
          # gene parts' percentage
          if (T) {
            # get the sorted matrix
            gene <- gene[, paste0("S", seq(min(as.numeric(gsub("S", "", colnames(gene)))), max(as.numeric(gsub("S", "", colnames(gene))))))]
            gene[1:4, 1:4]
            sapply(gene, function(x) {
              summary(x[x>0])
            })
            gene <- gene*1e4
            sapply(gene, function(x) {
              summary(x[x>0])
            })
            gene <- round(gene)
            gene <- log2(gene+1)
            
            gene <- CreateSeuratObject(counts = gene,
                                       assay = "gene",
                                       names.delim = "+",
                                       min.cells = 0, min.features = 0)
            gene@meta.data$orig.ident <- colnames(gene)
          }
          
          # emission
          if (T) {
            emission[1:4, 1:4]
            rownames(emission) <- emission[, 1]
            emission <- emission[, -c(1:2)]
            emission <- data.frame(t(emission))
            emission <- emission[, colnames(gene)]
            
            sapply(emission, function(x) {
              summary(x[x>0])
            })
            emission <- emission*1e2
            sapply(emission, function(x) {
              summary(x[x>0])
            })
            emission <- round(emission)
            emission <- log2(emission+1)
            
            emission <- CreateSeuratObject(counts = emission,
                                           assay = "emission",
                                           names.delim = "+",
                                           min.cells = 0, min.features = 0)
          }
          
          # merge the data
          if (T) {
            dat <- gene
            dat@assays$emission <- emission@assays$emission
            
            Assays(dat)
          }
        }
        
        # pre-qc for gene
        if (T) {
          DefaultAssay(dat) <- 'gene'
          
          dat[["gene"]]$data <- dat[["gene"]]$counts
          VariableFeatures(dat) <- rownames(dat[["gene"]])
          dat <- ScaleData(dat, features = rownames(dat))
          dat <- RunPCA(dat, features = VariableFeatures(dat), npcs = ncol(dat), verbose=F, approx=F)
          
          n_pc_gene <- ncol(dat@reductions$pca)
        }
        
        # pre-qc for emission
        if (T) {
          DefaultAssay(dat) <- 'emission'
          
          dat[["emission"]]$data <- dat[["emission"]]$counts
          VariableFeatures(dat) <- rownames(dat[["emission"]])
          dat <- ScaleData(dat, features = rownames(dat))
          dat <- RunPCA(dat, features = VariableFeatures(dat), npcs = nrow(dat), reduction.name = 'apca', verbose=F, approx=F)
          
          n_pc_emission <- ncol(dat@reductions$apca)
        }
        
        # WNN
        if (T) {
          dat <- FindMultiModalNeighbors(
            dat, 
            reduction.list = list("pca", "apca"), 
            dims.list = list(1:n_pc_gene, 1:n_pc_emission),
            k.nn = round(ncol(dat)*0.1),
            knn.range = round(ncol(dat)*0.5),
          )
        }
        
        # find cluster
        if (T) {
          dat <- FindClusters(dat, graph.name = "wsnn", algorithm = 3, resolution = resolutions, verbose = FALSE)
          
          library(clustree)
          p <- clustree(dat@meta.data, prefix = "wsnn_res.")
          
          cutoff <- ifelse(cell == "thp1", 3, 4)
          
          p + geom_hline(yintercept = cutoff)
          ggsave(paste0("results/2.pic/34.test_various_highly_informative_txs/", cell, "_hit_clustree_", sub_total, "_txs.png"), width = 7, height = 10)
          
          resolution <- ifelse(cell == "thp1", 2, 1.8)
          dat <- FindClusters(dat, graph.name = "wsnn", algorithm = 3, resolution = resolution, verbose = FALSE)
        }
        
        # format the res
        if (T) {
          cell_specific_group <- dat@meta.data
          cell_specific_group$cell <- cell
          cell_specific_group$total_hit <- sub_total
        }
        
        return(cell_specific_group)
      })
      subres <- do.call(rbind, subres)
      
      return(subres)
    })
    hit_dat <- do.call(rbind,res)
    
    qsave(hit_dat, file, nthreads = 6)
  }
  
  rm(target_tx)
}

# calculate the ari value compared with standard res
if (T) {
  ggdat <- lapply(c(cell1, cell2), function(cell) {
    # cell <- cell1
    
    # read now used dat
    if (T) {
      standard <- read.table(paste0("results/1.tab/24.", cell, "_states_group_info.csv"), header = T, sep = ",", fill = T, comment.char = "")
    }
    
    # select best resolution (和标准分群相同分辨率 / **和标准分群一样分5群** / 标准流程挑选最佳分辨率)
    if (T) {
      subdat <- hit_dat[hit_dat$cell == cell, ]
      compared <- lapply(total, function(n) {
        # n <- total[1]
        
        x <- subdat[subdat$total_hit == n, ]
        
        # 和标准分群一样分5群
        if (T) {
          n_cluster <- sapply(x[, paste0("wsnn_res.", resolutions)], function(c) {
            # c <- x$wsnn_res.0.9
            length(unique(c))
          })
          names(n_cluster) <- paste0("wsnn_res.", resolutions)
          
          resol <- resolutions[which((n_cluster == 5))]
          max_res <- max(resol)
          
          cluster <- x[, paste0("wsnn_res.", max_res)]
          cluster <- data.frame(
            total_hit = n, 
            cluster = cluster
          )
        }
        
        # 和标准分群相同分辨率
        if (F) {
          if (cell == cell1) {
            cluster <- data.frame(
              total_hit = n, 
              cluster = x$wsnn_res.2
            )
          }
          if (cell == cell2) {
            cluster <- data.frame(
              total_hit = n, 
              cluster = x$wsnn_res.1.8
            )
          }
        }
        
        # 标准流程挑选最佳分辨率
        if (F) {
          if (cell == cell1) {
            res <- ifelse(n %in% c(100, 200, 300, 1000), 1.8, 
                          ifelse(n %in% c(400, 500), 1.2, 2))
          }
          if (cell == cell2) {
            res <- ifelse(n %in% c(100, 200, 300, 3000, 3500, 4000), 2, 1.8)
          }
          
          cluster <- data.frame(
            total_hit = n, 
            cluster = x[, paste0("wsnn_res.", res)]
          )
        }
        
        return(cluster)
      })
      
      compared <- data.frame(do.call(rbind, compared))
    }
    
    # calculate the ari mat
    if (T) {
      res <- lapply(total, function(sub_total) {
        # sub_total <- total[1]
        
        comp <- compared$cluster[compared$total_hit == sub_total]
        ri <- flexclust::randIndex(table(standard$seurat_clusters, comp),correct = F)
        ari <- flexclust::randIndex(table(standard$seurat_clusters, comp),correct = T)
        
        res <- data.frame(
          total_hit = sub_total, 
          cell = cell, 
          ri = ri, 
          ari = ari
        )
        
        return(res)
      })
      res <- do.call(rbind, res)
      rownames(res) <- 1:nrow(res)
    }
    
    return(res)
  })
  
  # format the data
  if (T) {
    ggdat <- do.call(rbind, ggdat)
    qsave(ggdat, file = "data/saved_data/34.cell_specific_ari_matrix.qs", nthreads = 6)
    
    head(ggdat)
    
    mat <- dcast(ggdat, cell~total_hit, value.var = "ari")
    rownames(mat) <- mat[, 1]
    mat <- mat[, -1]
  }
  
  # heatmap
  if (T) {
    # color setting
    if (T) {
      library(RColorBrewer)
      colors <- brewer.pal(9, "Set1") 
      col_fun <- colorRamp2(c(0, 0.5, 1), c(colors[2], "white", colors[1]))
    }
    
    p <- Heatmap(as.matrix(mat), name = "ARI", col = col_fun,
             
             border_gp = gpar(col = "black"), rect_gp = gpar(col = "black"),
             
             show_heatmap_legend = TRUE, heatmap_legend_param = list(border = T),
             
             column_title = "column_title", row_title = "row_title",
             row_title_gp = gpar(fontsize = 13.2), column_title_gp = gpar(fontsize = 13.2),
             
             row_names_max_width = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), row_names_rot = 0,
             column_names_max_height = unit(6, "cm"), column_names_gp = gpar(fontsize = 10), column_names_rot = 60,
             
             cluster_rows = F, cluster_columns = F)
    p <- ggplotify::as.ggplot(p)
    print(p)
  }
}
