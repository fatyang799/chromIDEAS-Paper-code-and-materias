# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(stringr)
  library(reshape2)
  library(qs)
}

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
up_bin_num <- 5
down_bin_num <- 5
body_bin_num <- 10
rna_levels <- 10
n_ct <- 5

# define functions
if (T) {
  # common value
  if (T) {
    order <- c(
      paste0("U", up_bin_num:1), 
      "TSS", 
      paste0("B", 1:body_bin_num), 
      "TES", 
      paste0("D", 1:down_bin_num)
    )
  }
  
  # data prepare (script36)
  quantitle_mRNA <- function(cell_mrna) {
    # geneid: Select target genes. If 0, then all genes.
    
    # test data
    if (F) {
      cell_mrna <- mrna[[cell1]]
      cell_mrna <- cell_mrna[[1]]
    }
    
    # get subdat0 and Q1-Q10
    if (T) {
      cell_mrna_0 <- names(cell_mrna)[cell_mrna == 0]
      cell_mrna_non0 <- cell_mrna[cell_mrna > 0]
      
      cell_mrna_non0 <- lapply(seq(1, rna_levels), function(q) {
        min <- quantile(cell_mrna_non0, (q-1)/10)
        max <- quantile(cell_mrna_non0, q/10)
        
        if (max == max(cell_mrna_non0)) {
          target <- cell_mrna_non0[cell_mrna_non0>=min & cell_mrna_non0<=max]
        }
        if (max < max(cell_mrna_non0)) {
          target <- cell_mrna_non0[cell_mrna_non0>=min & cell_mrna_non0<max]
        }
        
        return(names(target))
      })
      names(cell_mrna_non0) <- paste0("Q", 1:rna_levels)
      cell_mrna_non0[["subdat0"]] <- cell_mrna_0
      
      cell_mrna <- cell_mrna_non0
      
      rm(cell_mrna_0, cell_mrna_non0)
    }
    
    return(cell_mrna)
  }
  
  # gpat pattern seq stat
  gene_body_pattern_seq_num_stat <- function(gpat_cs, ids) {
    # test dat
    if (F) {
      ids <- sample(gpat_cs$gene_id, 1000)
    }
    
    # get specific id dat
    if (T) {
      subdat <- gpat_cs[gpat_cs$gene_id %in% ids, ]
    }
    
    # calculate the gene body pattern type
    if (T) {
      stat <- apply(subdat[, order], 1, paste, collapse="+")
      names(stat) <- subdat$gene_id
    }
    
    # statistics the number of type for gene body pattern type
    if (T) {
      pattern_type_num <- as.data.frame.table(table(stat), stringsAsFactors = F)
      colnames(pattern_type_num) <- c("Pattern", "Number")
    }
    
    # sorting
    if (T) {
      pattern_type_num <- pattern_type_num[order(pattern_type_num$Number, decreasing=T), ]
      pattern_type_num$Top <- 1:nrow(pattern_type_num)
      rownames(pattern_type_num) <- pattern_type_num$Top
    }
    
    return(pattern_type_num)
  }
  gene_body_pattern_seq_num_batch_stat <- function(gpat_cs, IDs) {
    dat <- lapply(names(IDs), function(mrna) {
      # mrna <- names(IDs)[1]
      ids <- IDs[[mrna]]
      
      # calcluate the number of each gene body pattern seq
      if (T) {
        subdat <- gene_body_pattern_seq_num_stat(gpat_cs, ids)
        subdat$mRNA <- mrna
        subdat$Percentage <- subdat$Number/sum(subdat$Number)*100
      }
      
      return(subdat)
    })
    dat <- do.call(rbind, dat)
    return(dat)
  }
  
  # gpat pattern motif stat
  gene_body_ppm <- function(gpat_cs, ids, all_values) {
    # test dat
    if (F) {
      ids <- sample(gpat_cs$gene_id, 1000)
      all_values <- sort(as.numeric(unique(gpat_cs$TSS)))
    }
    
    # get specific data
    if (T) {
      # get specific data
      subdat <- gpat_cs[gpat_cs$gene_id %in% ids, ]
      
      # stat the ppm
      if (T) {
        ppm <- data.frame(sapply(subdat[, order], function(x) {
          # x <- subdat$U5
          stat <- table(x)/length(x)*100
          stat <- sapply(as.character(all_values), function(x) {
            dat <- ifelse(is.na(stat[x]), 0, stat[x])
            dat <- as.numeric(dat)
            return(dat)
          })
          
          return(stat)
        }))
      }
      
      # format the data
      if (T) {
        ppm$state <- rownames(ppm)
      }
    }
    
    # format the data
    if (T) {
      rownames(ppm) <- 1:nrow(ppm)
    }
    
    return(ppm)
  }
  gene_body_batch_ppm <- function(gpat_cs, IDs) {
    all_values <- sort(unique(as.numeric(gpat_cs$TSS)))
    
    dat <- lapply(names(IDs), function(mrna) {
      # mrna <- "Q1"
      ids <- IDs[[mrna]]
      
      ppm <- gene_body_ppm(gpat_cs, ids, all_values)
      ppm$mRNA <- mrna
      
      return(ppm)
    })
    
    # merge the data
    if (T) {
      dat <- do.call(rbind, dat)
    }
    
    return(dat)
  }
  gene_body_ppm_barplot <- function(gpat_motif_cs, state_prefix, title, legend_title, legend_nrow, color) {
    # test data
    if (F) {
      state_prefix <- "S"
      title <- "title"
      legend_title <- "chromatin states"
      legend_nrow <- 2
    }
    
    # format the data
    if (T) {
      ggdat <- melt(gpat_motif_cs, id.vars = c("state", "mRNA"), variable.name = "Loc", value.name = "Percentage")
      ggdat$state <- factor(paste0(state_prefix, ggdat$state), levels = paste0(state_prefix, sort(unique(as.numeric(ggdat$state)))))
      ggdat$mRNA <- factor(ggdat$mRNA, levels = c("subdat0", paste0("Q", 1:rna_levels)))
      ggdat$Loc <- factor(ggdat$Loc, levels = order)
    }
    
    # color setting
    if (F) {
      n_total <- length(levels(ggdat$state))
      n <- ceiling(n_total / 4)
      
      set.seed(799)
      color <- c(
        circlize::rand_color(n, hue = "red", luminosity = "bright"), 
        circlize::rand_color(n, hue = "orange", luminosity = "bright"), 
        circlize::rand_color(n, hue = "green", luminosity = "bright"), 
        circlize::rand_color(n, hue = "purple", luminosity = "bright")
      )
      color <- color[1:n_total]
      names(color) <- levels(ggdat$state)
    }
    
    # ggplot
    if (T) {
      head(ggdat)
      p <- ggplot(ggdat) +
        geom_bar(aes(x=Loc, y=Percentage, fill=state), color="black", position="stack", stat="identity") +
        ggtitle(title) +
        scale_fill_manual(values = color) +
        scale_x_discrete(name=NULL, breaks = c(paste0("U", up_bin_num), "TSS", "TES", paste0("D", down_bin_num))) +
        facet_grid(mRNA~.) +
        cowplot::theme_cowplot() +
        theme(panel.border = element_rect(color="black"), 
              strip.background = element_rect(fill=NA, color=NA), 
              legend.position = "bottom") +
        guides(fill = guide_legend(title=legend_title, nrow=legend_nrow, title.hjust=0.5))
    }
    
    return(p)
  }
  
  color_palette <- function(n, colors) {
    # test data
    if (F) {
      n <- 10
      colors <- c("#669BBC", "#C1121F", "#003049", "#FDF0D5")
    }
    
    # continuous colors
    cols <- colorRampPalette(c(colors))(n)
    
    return(cols)
  }
}

# calculate the number of txs with identical GSTpat seq
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir <- paste0("results/2.pic/43.GSPat_motif/", data_type)
      if (! dir.exists(dir)) {
        dir.create(dir, showWarnings = F, recursive = T)
      }
    }
    
    for (type in c("merged", "single")) {
      cat(paste0("\t", type, ": \n"))
      
      # cell specific data
      for (cell in c(cell1, cell2)) {
        cat(paste0("\t\t", cell, ": \n"))
        
        # test data
        if (F) {
          data_type <- "tx"
          type <- "merged"
          cell <- cell1
        }
        
        # get non_overlap-txs ID
        if (T) {
          file <- "data/saved_data/20.txs_with_nonoverlap_region.qs"
          
          non_overlap_txs <- qread(file, nthreads = 6)
          non_overlap_txs <- non_overlap_txs$tx_id
        }
        
        # load the gpat data: only non_overlap_txs
        if (T) {
          # GSpat
          if (T) {
            cell_dat <- function(cell, non_overlap_txs) {
              file <- paste0("data/saved_data/42.gene_state_pattern/all_txs/42.", data_type, "_level_profile_Body_", 
                             up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".state.", cell, ".qs")
              
              gpat_cs <- qread(file, nthreads = 6)
              gpat_cs <- gpat_cs[match(non_overlap_txs, gpat_cs$gene_id), ]
              
              return(gpat_cs)
            }
            gpat_cs <- cell_dat(cell, non_overlap_txs)
          }
          
          # GTpat
          if (T) {
            cell_dat <- function(type, cluster, cell) {
              file <- paste0("data/saved_data/42.gene_state_pattern/all_txs/42.", data_type, "_Body_", up_bin_num, "UP_", down_bin_num, 
                             "DW_", body_bin_num, "Body.bs", bin_size, ".", type, ".clusters(", cluster, ").", cell, ".qs")
              dat <- qread(file, nthreads = 6)
              dat <- dat[match(non_overlap_txs, dat$gene_id), ]
              
              variable <- ifelse(cluster %in% as.character(c(1:n_ct)), 
                                 paste0("gpat_ct", cluster), 
                                 paste0("gpat_", tolower(cluster)))
              
              return(list(variable, dat))
            }
            
            for (cluster in c("chromIDEAS", "kmeans", 1:n_ct)) {
              # cluster <- "chromIDEAS"
              
              dat <- cell_dat(type, cluster, cell)
              assign(dat[[1]], dat[[2]])
            }
          }
          
          rm(dat, cluster)
          
          # merge the data
          if (T) {
            gpat_dat <- list(
              "cs" = gpat_cs,
              "chromideas" = gpat_chromideas,
              "kmeans" = gpat_kmeans,
              "ct1" = gpat_ct1,
              "ct2" = gpat_ct2,
              "ct3" = gpat_ct3,
              "ct4" = gpat_ct4,
              "ct5" = gpat_ct5
            )
            rm(gpat_chromideas,gpat_kmeans,gpat_ct1,gpat_ct2,gpat_ct3,gpat_ct4,gpat_ct5)
          }
        }
        
        # get mRNA level data: IDs
        if (T) {
          file <- ifelse(data_type == "gene", 
                         paste0("data/saved_data/6.4_type_gene_expression_level_classification_", cell, ".qs"), 
                         paste0("data/saved_data/7.4_type_tx_expression_level_classification_", cell, ".qs"))
          
          cell_mrna <- qread(file, nthreads = 6)
          cell_mrna <- cell_mrna[[1]]
          
          cell_mrna <- cell_mrna[names(cell_mrna) %in% gpat_cs$gene_id]
          IDs <- quantitle_mRNA(cell_mrna)
          
          rm(cell_mrna)
        }
        
        # statistics the gpat pattern seq: non_overlap_txs
        if (T) {
          # statistics the number of gene body pattern seq
          if (T) {
            gpat_seq <- lapply(names(gpat_dat), function(method) {
              # method <- names(gpat_dat)[1]
              
              subdat <- gpat_dat[[method]]
              
              res <- gene_body_pattern_seq_num_batch_stat(subdat, IDs=IDs)
              
              res$method <- ifelse(method == "chromideas", "chromIDEAS", 
                                   ifelse(method == "kmeans", "kmeans", toupper(method)))
              
              return(res)
            })
          }
          
          # merge the results
          if (T) {
            gpat_seq_summary <- data.frame(do.call(rbind, gpat_seq))
          }
          
          # get only top5 results to visualization
          if (T) {
            cutoff <- 5
            
            ggdat <- gpat_seq_summary[gpat_seq_summary$Top<=cutoff+1, ]
            
            for (mrna in unique(ggdat$mRNA)) {
              for (method in unique(ggdat$method)) {
                # test data
                if (F) {
                  mrna <- unique(ggdat$mRNA)[1]
                  method <- unique(ggdat$method)[1]
                }
                
                subdat <- gpat_seq_summary[gpat_seq_summary$mRNA == mrna & gpat_seq_summary$method == method, ]
                subdat <- subdat[subdat$Top>cutoff, ]
                
                ggdat$Pattern[(ggdat$mRNA == mrna) & (ggdat$method == method) & (ggdat$Top == cutoff+1)] <- NA
                ggdat$Number[(ggdat$mRNA == mrna) & (ggdat$method == method) & (ggdat$Top == cutoff+1)] <- sum(subdat$Number)
                ggdat$Percentage[(ggdat$mRNA == mrna) & (ggdat$method == method) & (ggdat$Top == cutoff+1)] <- sum(subdat$Percentage)
              }
            }
            
            ggdat$Top[ggdat$Top == cutoff+1] <- paste0(">=", cutoff+1)
            
            rm(subdat, mrna, method)
          }
          
          write.table(ggdat, file = paste0("results/1.tab/43.", data_type, "_GSTPat_motif.", cell, ".", type, ".csv"), 
                      quote = F, sep = ",", col.names = T, row.names = F)
          
          # format the data
          if (T) {
            head(ggdat)
            
            ggdat$mRNA <- factor(ggdat$mRNA, levels = c("subdat0", paste0("Q", 1:rna_levels)))
            ggdat$method <- factor(ggdat$method, levels = c("CS", "chromIDEAS", "kmeans", paste0("CT", 1:5)))
            ggdat$Top <- factor(ggdat$Top, levels = c(paste0(">=", cutoff+1), cutoff:1))
            ggdat$group <- paste0(ggdat$method, "@", ggdat$mRNA)
          }
          
          # color setting
          if (T) {
            n <- ifelse(cutoff %% 3 == 0, (cutoff+0)/3, 
                        ifelse(cutoff %% 3 == 1, (cutoff+2)/3, (cutoff+1)/3))
            
            set.seed(799)
            color <- c(
              circlize::rand_color(n, hue = "red", luminosity = "bright"), 
              circlize::rand_color(n, hue = "orange", luminosity = "bright"), 
              circlize::rand_color(n, hue = "blue", luminosity = "bright")
            )
            color <- color[1:cutoff]
            color[cutoff+1] <- "#F2F2F2"
            
            names(color) <- rev(levels(ggdat$Top))
          }
          
          # ggplot
          if (T) {
            head(ggdat)
            p <- ggplot(ggdat) +
              geom_bar(aes(x=mRNA, y=Percentage, fill=Top), colour="black", position="stack", stat="identity") +
              scale_fill_manual(values=color) +
              ggtitle(paste0(data_type, ": ", cell)) +
              facet_grid(method~.) +
              cowplot::theme_cowplot() +
              theme(legend.position = "bottom",
                    strip.background = element_rect(fill=NA, color=NA)) +
              guides(fill = guide_legend(title = "Top Rank", nrow = 1, reverse = TRUE))
          }
          
          file <- paste0(dir, "/43.top", cutoff, "_enriched_GSTPat_seq_percentage_in_", data_type, "s_Body.", cell, ".", type, ".pdf")
          ggsave(filename = file, plot = p, width = 12, height = 10)
          
          rm(p, n ,color, cutoff, gpat_seq_summary)
        }
      }
      
      # merged data
      if (type == "merged") {
        cat(paste0("\t\tmerged: \n"))
        
        # test data
        if (F) {
          data_type <- "tx"
          type <- "merged"
        }
        
        # get non_overlap-txs ID
        if (T) {
          file <- "data/saved_data/20.txs_with_nonoverlap_region.qs"
          
          non_overlap_txs <- qread(file, nthreads = 6)
          non_overlap_txs <- non_overlap_txs$tx_id
        }
        
        # load the gpat data: only non_overlap_txs
        if (T) {
          # GSpat
          if (T) {
            cell_dat <- function(cell, non_overlap_txs) {
              file <- paste0("data/saved_data/42.gene_state_pattern/all_txs/42.", data_type, "_level_profile_Body_", 
                             up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".state.", cell, ".qs")
              
              gpat_cs <- qread(file, nthreads = 6)
              gpat_cs <- gpat_cs[match(non_overlap_txs, gpat_cs$gene_id), ]
              
              return(gpat_cs)
            }
            dat1 <- cell_dat(cell1, non_overlap_txs)
            dat2 <- cell_dat(cell2, non_overlap_txs)
            dat1$gene_id <- paste0(dat1$gene_id, "@", cell1)
            dat2$gene_id <- paste0(dat2$gene_id, "@", cell2)
            gpat_cs <- rbind(dat1, dat2)
            
            rm(dat1, dat2)
          }
          
          # GTpat
          if (T) {
            cell_dat <- function(type, cluster, cell) {
              file <- paste0("data/saved_data/42.gene_state_pattern/all_txs/42.", data_type, "_Body_", up_bin_num, "UP_", down_bin_num, 
                             "DW_", body_bin_num, "Body.bs", bin_size, ".", type, ".clusters(", cluster, ").", cell, ".qs")
              dat <- qread(file, nthreads = 6)
              dat <- dat[match(non_overlap_txs, dat$gene_id), ]
              
              variable <- ifelse(cluster %in% as.character(c(1:n_ct)), 
                                 paste0("gpat_ct", cluster), 
                                 paste0("gpat_", tolower(cluster)))
              
              return(list(variable, dat))
            }
            
            for (cluster in c("chromIDEAS", "kmeans", 1:n_ct)) {
              # cluster <- "chromIDEAS"
              
              variable <- cell_dat(type, cluster, cell1)
              dat1 <- variable[[2]]
              dat2 <- cell_dat(type, cluster, cell2)[[2]]
              variable <- variable[[1]]
              
              dat1$gene_id <- paste0(dat1$gene_id, "@", cell1)
              dat2$gene_id <- paste0(dat2$gene_id, "@", cell2)
              dat <- rbind(dat1, dat2)
              
              assign(variable, dat)
            }
          }
          
          rm(variable, dat1, dat2, dat, cluster)
          
          # merge the data
          if (T) {
            gpat_dat <- list(
              "cs" = gpat_cs,
              "chromideas" = gpat_chromideas,
              "kmeans" = gpat_kmeans,
              "ct1" = gpat_ct1,
              "ct2" = gpat_ct2,
              "ct3" = gpat_ct3,
              "ct4" = gpat_ct4,
              "ct5" = gpat_ct5
            )
            rm(gpat_cs, gpat_chromideas,gpat_kmeans,gpat_ct1,gpat_ct2,gpat_ct3,gpat_ct4,gpat_ct5)
          }
        }
        
        # get mRNA level data: IDs
        if (T) {
          rna_Q <- function(cell, non_overlap_txs) {
            # cell <- cell1
            file <- ifelse(data_type == "gene", 
                           paste0("data/saved_data/6.4_type_gene_expression_level_classification_", cell, ".qs"), 
                           paste0("data/saved_data/7.4_type_tx_expression_level_classification_", cell, ".qs"))
            
            cell_mrna <- qread(file, nthreads = 6)
            cell_mrna <- cell_mrna[[1]]
            
            cell_mrna <- cell_mrna[names(cell_mrna) %in% non_overlap_txs]
            names(cell_mrna) <- paste0(names(cell_mrna), "@", cell)
            IDs <- quantitle_mRNA(cell_mrna)
            
            return(IDs)
          }
          
          dat1 <- rna_Q(cell1, non_overlap_txs)
          dat2 <- rna_Q(cell2, non_overlap_txs)
          
          IDs <- lapply(names(dat1), function(x) {
            # x <- names(dat1)[1]
            
            subdat1 <- dat1[[x]]
            subdat2 <- dat2[[x]]
            
            return(c(subdat1, subdat2))
          })
          names(IDs) <- names(dat1)
          
          rm(dat1, dat2)
        }
        
        # statistics the gpat pattern seq: non_overlap_txs
        if (T) {
          # statistics the number of gene body pattern seq
          if (T) {
            gpat_seq <- lapply(names(gpat_dat), function(method) {
              # method <- names(gpat_dat)[1]
              
              subdat <- gpat_dat[[method]]
              
              res <- gene_body_pattern_seq_num_batch_stat(subdat, IDs=IDs)
              
              res$method <- ifelse(method == "chromideas", "chromIDEAS", 
                                   ifelse(method == "kmeans", "kmeans", toupper(method)))
              
              return(res)
            })
          }
          
          # merge the results
          if (T) {
            gpat_seq_summary <- data.frame(do.call(rbind, gpat_seq))
          }
          
          # get only top5 results to visualization
          if (T) {
            cutoff <- 5
            
            ggdat <- gpat_seq_summary[gpat_seq_summary$Top<=cutoff+1, ]
            
            for (mrna in unique(ggdat$mRNA)) {
              for (method in unique(ggdat$method)) {
                # test data
                if (F) {
                  mrna <- unique(ggdat$mRNA)[1]
                  method <- unique(ggdat$method)[1]
                }
                
                subdat <- gpat_seq_summary[gpat_seq_summary$mRNA == mrna & gpat_seq_summary$method == method, ]
                subdat <- subdat[subdat$Top>cutoff, ]
                
                ggdat$Pattern[(ggdat$mRNA == mrna) & (ggdat$method == method) & (ggdat$Top == cutoff+1)] <- NA
                ggdat$Number[(ggdat$mRNA == mrna) & (ggdat$method == method) & (ggdat$Top == cutoff+1)] <- sum(subdat$Number)
                ggdat$Percentage[(ggdat$mRNA == mrna) & (ggdat$method == method) & (ggdat$Top == cutoff+1)] <- sum(subdat$Percentage)
              }
            }
            
            ggdat$Top[ggdat$Top == cutoff+1] <- paste0(">=", cutoff+1)
            
            rm(subdat, mrna, method)
          }
          
          write.table(ggdat, file = paste0("results/1.tab/43.", data_type, "_GSTPat_motif.merged.", type, ".csv"), 
                      quote = F, sep = ",", col.names = T, row.names = F)
          
          # format the data
          if (T) {
            head(ggdat)
            
            ggdat$mRNA <- factor(ggdat$mRNA, levels = c("subdat0", paste0("Q", 1:rna_levels)))
            ggdat$method <- factor(ggdat$method, levels = c("CS", "chromIDEAS", "kmeans", paste0("CT", 1:5)))
            ggdat$Top <- factor(ggdat$Top, levels = c(paste0(">=", cutoff+1), cutoff:1))
            ggdat$group <- paste0(ggdat$method, "@", ggdat$mRNA)
          }
          
          # color setting
          if (T) {
            n <- ifelse(cutoff %% 3 == 0, (cutoff+0)/3, 
                        ifelse(cutoff %% 3 == 1, (cutoff+2)/3, (cutoff+1)/3))
            
            set.seed(799)
            color <- c(
              circlize::rand_color(n, hue = "red", luminosity = "bright"), 
              circlize::rand_color(n, hue = "orange", luminosity = "bright"), 
              circlize::rand_color(n, hue = "blue", luminosity = "bright")
            )
            color <- color[1:cutoff]
            color[cutoff+1] <- "#F2F2F2"
            
            names(color) <- rev(levels(ggdat$Top))
          }
          
          # ggplot
          if (T) {
            head(ggdat)
            p <- ggplot(ggdat) +
              geom_bar(aes(x=mRNA, y=Percentage, fill=Top), colour="black", position="stack", stat="identity") +
              scale_fill_manual(values=color) +
              ggtitle(paste0(data_type, ": merged")) +
              facet_grid(method~.) +
              cowplot::theme_cowplot() +
              theme(legend.position = "bottom",
                    strip.background = element_rect(fill=NA, color=NA)) +
              guides(fill = guide_legend(title = "Top Rank", nrow = 1, reverse = TRUE))
          }
          
          file <- paste0(dir, "/43.top", cutoff, "_enriched_GSTPat_seq_percentage_in_", data_type, "s_Body.merged.", type, ".pdf")
          ggsave(filename = file, plot = p, width = 12, height = 10)
          
          rm(p, n ,color, cutoff, gpat_seq_summary)
        }
      }
    }
  }
}

# state motif for all txs
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir <- paste0("results/2.pic/43.GSPat_motif/", data_type)
      if (! dir.exists(dir)) {
        dir.create(dir, showWarnings = F, recursive = T)
      }
    }
    
    for (type in c("merged", "single")) {
      cat(paste0("\t", type, ": \n"))
      
      # cell specific data
      for (cell in c(cell1, cell2)) {
        cat(paste0("\t\t", cell, ": \n"))
        
        # test data
        if (F) {
          data_type <- "tx"
          type <- "merged"
          cell <- cell1
        }
        
        # load the gpat data: all txs
        if (T) {
          # GSpat
          if (T) {
            cell_dat <- function(cell) {
              file <- paste0("data/saved_data/42.gene_state_pattern/all_txs/42.", data_type, "_level_profile_Body_", 
                             up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".state.", cell, ".qs")
              
              gpat_cs <- qread(file, nthreads = 6)
              
              return(gpat_cs)
            }
            gpat_cs <- cell_dat(cell)
          }
          
          # GTpat
          if (T) {
            cell_dat <- function(type, cluster, cell) {
              file <- paste0("data/saved_data/42.gene_state_pattern/all_txs/42.", data_type, "_Body_", up_bin_num, "UP_", down_bin_num, 
                             "DW_", body_bin_num, "Body.bs", bin_size, ".", type, ".clusters(", cluster, ").", cell, ".qs")
              dat <- qread(file, nthreads = 6)
              
              variable <- ifelse(cluster %in% as.character(c(1:n_ct)), 
                                 paste0("gpat_ct", cluster), 
                                 paste0("gpat_", tolower(cluster)))
              
              return(list(variable, dat))
            }
            
            for (cluster in c("chromIDEAS", "kmeans", 1:n_ct)) {
              # cluster <- "chromIDEAS"
              
              dat <- cell_dat(type, cluster, cell)
              assign(dat[[1]], dat[[2]])
            }
          }
          
          rm(dat, cluster)
          
          # merge the data
          if (T) {
            gpat_dat <- list(
              "cs" = gpat_cs,
              "chromideas" = gpat_chromideas,
              "kmeans" = gpat_kmeans,
              "ct1" = gpat_ct1,
              "ct2" = gpat_ct2,
              "ct3" = gpat_ct3,
              "ct4" = gpat_ct4,
              "ct5" = gpat_ct5
            )
            rm(gpat_chromideas,gpat_kmeans,gpat_ct1,gpat_ct2,gpat_ct3,gpat_ct4,gpat_ct5)
          }
        }
        
        # get mRNA level data: IDs
        if (T) {
          file <- ifelse(data_type == "gene", 
                         paste0("data/saved_data/6.4_type_gene_expression_level_classification_", cell, ".qs"), 
                         paste0("data/saved_data/7.4_type_tx_expression_level_classification_", cell, ".qs"))
          
          cell_mrna <- qread(file, nthreads = 6)
          cell_mrna <- cell_mrna[[1]]
          
          cell_mrna <- cell_mrna[names(cell_mrna) %in% gpat_cs$gene_id]
          IDs <- quantitle_mRNA(cell_mrna)
          
          rm(cell_mrna)
        }
        
        # statistics the gpat pattern motif: all txs
        if (T) {
          # calculate the ppm for each mRNA level
          if (T) {
            gpat_motif_cs <- gene_body_batch_ppm(gpat_cs, IDs)
            gpat_motif <- lapply(gpat_dat, gene_body_batch_ppm, IDs=IDs)
            names(gpat_motif)
          }
          
          # ggplot
          if (T) {
            # color setting
            if (T) {
              col_state <- color_palette(length(unique(gpat_motif$cs$state)), c("#6a994e", "#fca311", "#C1121F", "#d4a373", "#3a0ca3"))
              col_cluster <- color_palette(length(unique(gpat_motif$chromideas$state)), c("#C1121F", "#fca311", "#6a994e", "#3a0ca3", "#d4a373"))
            }
            
            # gpat_motif_cs
            if (T) {
              p <- gene_body_ppm_barplot(gpat_motif[[1]], 
                                         state_prefix="S", 
                                         title=paste0(data_type, ": ", cell, " (cs)"), 
                                         legend_title="Chromatin\nStates", 
                                         legend_nrow=2, 
                                         col_state)
              file <- paste0(dir, "/43.ppm_barplot_chromatinStates_in_", data_type, "s_Body.", cell, ".", type, ".pdf")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
            
            # gpat_motif_chromideas
            if (T) {
              p <- gene_body_ppm_barplot(gpat_motif[[2]], 
                                         state_prefix="C", 
                                         title=paste0(data_type, ": ", cell, " (chromideas)"), 
                                         legend_title="chromIDEAS cluster", 
                                         legend_nrow=1, 
                                         col_cluster)
              file <- paste0(dir, "/43.ppm_barplot_chromIDEAS_in_", data_type, "s_Body.", cell, ".", type, ".pdf")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
            
            # gpat_motif_kmeans
            if (T) {
              p <- gene_body_ppm_barplot(gpat_motif[[3]], 
                                         state_prefix="C", 
                                         title=paste0(data_type, ": ", cell, " (kmeans)"), 
                                         legend_title="kmeans cluster", 
                                         legend_nrow=1, 
                                         col_cluster)
              file <- paste0(dir, "/43.ppm_barplot_kmeas_in_", data_type, "s_Body.", cell, ".", type, ".pdf")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
            
            # gpat_motif_ct1
            if (T) {
              p <- gene_body_ppm_barplot(gpat_motif[[4]], 
                                         state_prefix="C", 
                                         title=paste0(data_type, ": ", cell, " (ct1)"), 
                                         legend_title="ct1", 
                                         legend_nrow=1, 
                                         col_cluster)
              file <- paste0(dir, "/43.ppm_barplot_ct1_in_", data_type, "s_Body.", cell, ".", type, ".pdf")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
            
            # gpat_motif_ct2
            if (T) {
              p <- gene_body_ppm_barplot(gpat_motif[[5]], 
                                         state_prefix="C", 
                                         title=paste0(data_type, ": ", cell, " (ct2)"), 
                                         legend_title="ct2", 
                                         legend_nrow=1, 
                                         col_cluster)
              file <- paste0(dir, "/43.ppm_barplot_ct2_in_", data_type, "s_Body.", cell, ".", type, ".pdf")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
            
            # gpat_motif_ct3
            if (T) {
              p <- gene_body_ppm_barplot(gpat_motif[[6]], 
                                         state_prefix="C", 
                                         title=paste0(data_type, ": ", cell, " (ct3)"), 
                                         legend_title="ct3", 
                                         legend_nrow=1, 
                                         col_cluster)
              file <- paste0(dir, "/43.ppm_barplot_ct3_in_", data_type, "s_Body.", cell, ".", type, ".pdf")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
            
            # gpat_motif_ct4
            if (T) {
              p <- gene_body_ppm_barplot(gpat_motif[[7]], 
                                         state_prefix="C", 
                                         title=paste0(data_type, ": ", cell, " (ct4)"), 
                                         legend_title="ct4", 
                                         legend_nrow=1, 
                                         col_cluster)
              file <- paste0(dir, "/43.ppm_barplot_ct4_in_", data_type, "s_Body.", cell, ".", type, ".pdf")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
            
            # gpat_motif_ct5
            if (T) {
              p <- gene_body_ppm_barplot(gpat_motif[[8]], 
                                         state_prefix="C", 
                                         title=paste0(data_type, ": ", cell, " (ct5)"), 
                                         legend_title="ct5", 
                                         legend_nrow=1, 
                                         col_cluster)
              file <- paste0(dir, "/43.ppm_barplot_ct5_in_", data_type, "s_Body.", cell, ".", type, ".pdf")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
          }
        }
      }
      
      # merged data
      if (type == "merged") {
        cat(paste0("\t\tmerged: \n"))
        
        # test data
        if (F) {
          data_type <- "tx"
          type <- "merged"
        }
        
        # load the gpat data: all txs
        if (T) {
          # GSpat
          if (T) {
            cell_dat <- function(cell) {
              file <- paste0("data/saved_data/42.gene_state_pattern/all_txs/42.", data_type, "_level_profile_Body_", 
                             up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".state.", cell, ".qs")
              
              gpat_cs <- qread(file, nthreads = 6)
              
              return(gpat_cs)
            }
            dat1 <- cell_dat(cell1)
            dat2 <- cell_dat(cell2)
            dat1$gene_id <- paste0(dat1$gene_id, "@", cell1)
            dat2$gene_id <- paste0(dat2$gene_id, "@", cell2)
            gpat_cs <- rbind(dat1, dat2)
            
            rm(dat1, dat2)
          }
          
          # GTpat
          if (T) {
            cell_dat <- function(type, cluster, cell) {
              file <- paste0("data/saved_data/42.gene_state_pattern/all_txs/42.", data_type, "_Body_", up_bin_num, "UP_", down_bin_num, 
                             "DW_", body_bin_num, "Body.bs", bin_size, ".", type, ".clusters(", cluster, ").", cell, ".qs")
              dat <- qread(file, nthreads = 6)
              
              variable <- ifelse(cluster %in% as.character(c(1:n_ct)), 
                                 paste0("gpat_ct", cluster), 
                                 paste0("gpat_", tolower(cluster)))
              
              return(list(variable, dat))
            }
            
            for (cluster in c("chromIDEAS", "kmeans", 1:n_ct)) {
              # cluster <- "chromIDEAS"
              
              variable <- cell_dat(type, cluster, cell1)
              dat1 <- variable[[2]]
              dat2 <- cell_dat(type, cluster, cell2)[[2]]
              variable <- variable[[1]]
              
              dat1$gene_id <- paste0(dat1$gene_id, "@", cell1)
              dat2$gene_id <- paste0(dat2$gene_id, "@", cell2)
              dat <- rbind(dat1, dat2)
              
              assign(variable, dat)
            }
          }
          
          rm(variable, dat1, dat2, dat, cluster)
          
          # merge the data
          if (T) {
            gpat_dat <- list(
              "cs" = gpat_cs,
              "chromideas" = gpat_chromideas,
              "kmeans" = gpat_kmeans,
              "ct1" = gpat_ct1,
              "ct2" = gpat_ct2,
              "ct3" = gpat_ct3,
              "ct4" = gpat_ct4,
              "ct5" = gpat_ct5
            )
            rm(gpat_chromideas,gpat_kmeans,gpat_ct1,gpat_ct2,gpat_ct3,gpat_ct4,gpat_ct5)
          }
        }
        
        # get mRNA level data: IDs
        if (T) {
          rna_Q <- function(cell) {
            # cell <- cell1
            file <- ifelse(data_type == "gene", 
                           paste0("data/saved_data/6.4_type_gene_expression_level_classification_", cell, ".qs"), 
                           paste0("data/saved_data/7.4_type_tx_expression_level_classification_", cell, ".qs"))
            
            cell_mrna <- qread(file, nthreads = 6)
            cell_mrna <- cell_mrna[[1]]
            
            names(cell_mrna) <- paste0(names(cell_mrna), "@", cell)
            cell_mrna <- cell_mrna[names(cell_mrna) %in% gpat_dat$cs$gene_id]
            IDs <- quantitle_mRNA(cell_mrna)
            
            return(IDs)
          }
          
          dat1 <- rna_Q(cell1)
          dat2 <- rna_Q(cell2)
          
          IDs <- lapply(names(dat1), function(x) {
            # x <- names(dat1)[1]
            
            subdat1 <- dat1[[x]]
            subdat2 <- dat2[[x]]
            
            return(c(subdat1, subdat2))
          })
          names(IDs) <- names(dat1)
          
          rm(dat1, dat2)
        }
        
        # statistics the gpat pattern motif: all txs
        if (T) {
          # calculate the ppm for each mRNA level
          if (T) {
            gpat_motif_cs <- gene_body_batch_ppm(gpat_cs, IDs)
            gpat_motif <- lapply(gpat_dat, gene_body_batch_ppm, IDs=IDs)
            names(gpat_motif)
          }
          
          # ggplot
          if (T) {
            # color setting
            if (T) {
              col_state <- color_palette(length(unique(gpat_motif$cs$state)), c("#6a994e", "#fca311", "#C1121F", "#d4a373", "#3a0ca3"))
              col_cluster <- color_palette(length(unique(gpat_motif$chromideas$state)), c("#C1121F", "#fca311", "#6a994e", "#3a0ca3", "#d4a373"))
            }
            
            # gpat_motif_cs
            if (T) {
              p <- gene_body_ppm_barplot(gpat_motif[[1]], 
                                         state_prefix="S", 
                                         title=paste0(data_type, ": merged (cs)"), 
                                         legend_title="Chromatin\nStates", 
                                         legend_nrow=2, 
                                         col_state)
              file <- paste0(dir, "/43.ppm_barplot_chromatinStates_in_", data_type, "s_Body.merged.", type, ".pdf")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
            
            # gpat_motif_chromideas
            if (T) {
              p <- gene_body_ppm_barplot(gpat_motif[[2]], 
                                         state_prefix="C", 
                                         title=paste0(data_type, ": merged (chromideas)"), 
                                         legend_title="chromIDEAS cluster", 
                                         legend_nrow=1, 
                                         col_cluster)
              file <- paste0(dir, "/43.ppm_barplot_chromIDEAS_in_", data_type, "s_Body.merged.", type, ".pdf")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
            
            # gpat_motif_kmeans
            if (T) {
              p <- gene_body_ppm_barplot(gpat_motif[[3]], 
                                         state_prefix="C", 
                                         title=paste0(data_type, ": merged (kmeans)"), 
                                         legend_title="kmeans cluster", 
                                         legend_nrow=1, 
                                         col_cluster)
              file <- paste0(dir, "/43.ppm_barplot_kmeas_in_", data_type, "s_Body.merged.", type, ".pdf")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
            
            # gpat_motif_ct1
            if (T) {
              p <- gene_body_ppm_barplot(gpat_motif[[4]], 
                                         state_prefix="C", 
                                         title=paste0(data_type, ": merged (ct1)"), 
                                         legend_title="ct1", 
                                         legend_nrow=1, 
                                         col_cluster)
              file <- paste0(dir, "/43.ppm_barplot_ct1_in_", data_type, "s_Body.merged.", type, ".pdf")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
            
            # gpat_motif_ct2
            if (T) {
              p <- gene_body_ppm_barplot(gpat_motif[[5]], 
                                         state_prefix="C", 
                                         title=paste0(data_type, ": merged (ct2)"), 
                                         legend_title="ct2", 
                                         legend_nrow=1, 
                                         col_cluster)
              file <- paste0(dir, "/43.ppm_barplot_ct2_in_", data_type, "s_Body.merged.", type, ".pdf")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
            
            # gpat_motif_ct3
            if (T) {
              p <- gene_body_ppm_barplot(gpat_motif[[6]], 
                                         state_prefix="C", 
                                         title=paste0(data_type, ": merged (ct3)"), 
                                         legend_title="ct3", 
                                         legend_nrow=1, 
                                         col_cluster)
              file <- paste0(dir, "/43.ppm_barplot_ct3_in_", data_type, "s_Body.merged.", type, ".pdf")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
            
            # gpat_motif_ct4
            if (T) {
              p <- gene_body_ppm_barplot(gpat_motif[[7]], 
                                         state_prefix="C", 
                                         title=paste0(data_type, ": merged (ct4)"), 
                                         legend_title="ct4", 
                                         legend_nrow=1, 
                                         col_cluster)
              file <- paste0(dir, "/43.ppm_barplot_ct4_in_", data_type, "s_Body.merged.", type, ".pdf")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
            
            # gpat_motif_ct5
            if (T) {
              p <- gene_body_ppm_barplot(gpat_motif[[8]], 
                                         state_prefix="C", 
                                         title=paste0(data_type, ": merged (ct5)"), 
                                         legend_title="ct5", 
                                         legend_nrow=1, 
                                         col_cluster)
              file <- paste0(dir, "/43.ppm_barplot_ct5_in_", data_type, "s_Body.merged.", type, ".pdf")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
          }
        }
      }
    }
  }
}
