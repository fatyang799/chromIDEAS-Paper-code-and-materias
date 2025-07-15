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
    
    # all mRNA levels
    mrna <- c("subdat0", paste0("Q", 1:rna_levels))
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
  
  # gpat pattern seq stat (script 43)
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
  
  calculate_cor_with_target <- function(gpat_motif_cs) {
    # calculate the cor compare with target
    if (T) {
      cor_compare_with_target <- lapply(mrna, function(base_rna) {
        # base_rna <- "Q1"
        
        # get base data
        if (T) {
          base <- gpat_motif_cs[gpat_motif_cs$mRNA == base_rna, ]
        }
        
        # remain RNA levels to be compared
        remain_rna <- setdiff(mrna, base_rna)
        
        # calculate the cor
        if (T) {
          cor <- data.frame(sapply(remain_rna, function(rna) {
            # rna <- remain_rna[1]
            loc_cor <- sapply(order, function(loc) {
              # loc <- order[1]
              base_loc <- base[, loc]
              
              remain_rna_dat <- gpat_motif_cs[gpat_motif_cs$mRNA == rna, ]
              remain_rna_dat <- remain_rna_dat[match(base$state, remain_rna_dat$state), ]
              remain_rna_loc <- remain_rna_dat[, loc]
              
              cor(base_loc, remain_rna_loc)
            })
            
            return(loc_cor)
          }))
        }
        
        # format the data
        if (T) {
          cor$Loc <- rownames(cor)
          
          format_cor <- melt(cor, id.vars = "Loc", variable.name = "mRNA", value.name = "Cor", factorsAsStrings = F)
          format_cor$Loc <- as.character(format_cor$Loc)
          format_cor$mRNA <- as.character(format_cor$mRNA)
          format_cor$base <- base_rna
        }
        
        return(format_cor)
      })
    }
    
    # merge the data
    if (T) {
      cor_compare_with_target <- data.frame(do.call(rbind, cor_compare_with_target))
    }
    
    return(cor_compare_with_target)
  }
  plot_cor_with_target <- function(gpat_motif_cs, title) {
    # test data
    if (F) {
      title <- "title"
    }
    
    # get formated data
    if (T) {
      ggdat <- calculate_cor_with_target(gpat_motif_cs)
    }
    
    # format the data
    if (T) {
      ggdat$mRNA <- factor(ggdat$mRNA, levels = mrna)
      ggdat$base <- factor(ggdat$base, levels = mrna)
      ggdat$Loc <- factor(ggdat$Loc, levels = order)
    }
    
    # color setting
    if (T) {
      col_fun <- circlize::colorRamp2(c(1, length(mrna)), c("blue", "red"))
    }
    
    # ggplot
    if (T) {
      head(ggdat)
      
      p <- ggplot(ggdat) +
        geom_point(aes(x=Loc, y=Cor), size=1.5) +
        geom_line(aes(x=Loc, y=Cor, group=mRNA, color=mRNA), linewidth=1.2) +
        scale_x_discrete(name=NULL, breaks = c(paste0("U", up_bin_num), "TSS", "TES", paste0("D", down_bin_num))) +
        scale_color_manual(values = col_fun(1:length(mrna)),
                           breaks = mrna) +
        ggtitle(title) +
        facet_wrap(~base, ncol = 4) +
        cowplot::theme_cowplot() +
        theme(panel.border = element_rect(color="black"), 
              strip.background = element_rect(fill=NA, color=NA), 
              legend.position = "bottom") +
        guides(colour = guide_legend(nrow=1))
    }
    
    return(p)
  }
  plot_cor_similarity <- function(gpat_motif_cs, 
                                  gpat_motif_others, 
                                  names, title, type="rain_plot", method="pearson") {
    # test data
    if (F) {
      gpat_motif_others <- gpat_motif[2:8]
      names <- c("chromIDEAS", "kmeans", paste0("CT", 1:n_ct))
      title <- "title"
      type <- "barplot"
      method <- "pearson"
      type <- "rain_plot"
    }
    
    # get data
    if (T) {
      dat1 <- calculate_cor_with_target(gpat_motif_cs)
      dat2 <- lapply(gpat_motif_others, calculate_cor_with_target)
      names(dat2) <- names
    }
    
    # calculate the similarity of correlation
    if (T) {
      ggdat <- lapply(names, function(subname) {
        # subname <- names[1]
        subdat <- dat2[[subname]]
        
        base_rna_similarity <- lapply(mrna, function(base_rna) {
          # base_rna <- mrna[1]
          
          # remain RNA levels to be compared
          remain_rna <- setdiff(mrna, base_rna)
          
          # calculate the cor
          if (T) {
            cor_similarity <- sapply(remain_rna, function(r_rna) {
              # r_rna <- remain_rna[1]
              
              subdat1 <- dat1[dat1$base == base_rna & dat1$mRNA == r_rna, ]
              subdat2 <- subdat[subdat$base == base_rna & subdat$mRNA == r_rna, ]
              
              subdat2 <- subdat2[match(subdat1$Loc, subdat2$Loc), ]
              
              similarity <- cor(subdat1$Cor, subdat2$Cor, method = method)
              
              return(similarity)
            })
          }
          
          # format the cor similarity
          if (T) {
            cor_similarity <- data.frame(similarity = cor_similarity)
            cor_similarity$method <- subname
            cor_similarity$base <- base_rna
            cor_similarity$mRNA <- remain_rna
          }
          
          return(cor_similarity)
        })
        base_rna_similarity <- do.call(rbind, base_rna_similarity)
        
        return(base_rna_similarity)
      })
    }
    
    # format the data
    if (T) {
      ggdat <- data.frame(do.call(rbind, ggdat))
      rownames(ggdat) <- 1:nrow(ggdat)
      head(ggdat)
      
      ggdat$mRNA <- factor(ggdat$mRNA, levels = mrna)
      ggdat$base <- factor(ggdat$base, levels = mrna)
      ggdat$method <- factor(ggdat$method, levels = names)
      ggdat$similarity <- as.numeric(ggdat$similarity)
    }
    
    # statistis group
    if (T) {
      ref <- "chromIDEAS"
      remain <- setdiff(levels(ggdat$method), ref)
      my_comparisons <- lapply(remain, function(x) {
        c(ref, x)
      })
      
      print(do.call(rbind, tapply(ggdat$similarity, ggdat$method, summary)))
    }
    
    # ggplot
    if (type == "barplot") {
      head(ggdat)
      
      p <- ggplot(ggdat, aes(x=method, y=similarity, group=method)) +
        geom_violin(aes(fill=method, group=method), scale = "width") +
        geom_boxplot(aes(group=method), fill="NA", color="black", width=0.5) +
        ggpubr::stat_compare_means(method="wilcox.test", paired = T, label="p.signif", comparisons = my_comparisons) +
        ggpubr::stat_compare_means(method="anova", label.y=2.2) +
        scale_x_discrete(name=NULL) +
        scale_y_continuous(name="Data Fidelity Level", limits = c(NA, 2.3)) +
        ggtitle(title) +
        cowplot::theme_cowplot() +
        theme(panel.border = element_rect(color="black"), 
              strip.background = element_rect(fill=NA, color=NA), 
              legend.position = "bottom") +
        guides(fill = guide_legend(nrow=1))
    }
    
    # ggplot
    if (type == "rain_plot") {
      head(ggdat)
      
      # ggplot
      p <- ggplot(ggdat, aes(x=method, y=similarity, group=method)) +
        gghalves::geom_half_violin(aes(fill=method, group=method), side = "r", color="black", alpha=0.5) +
        gghalves::geom_half_boxplot(aes(fill=method, group=method), side = "r", errorbar.draw = FALSE, width=0.2, linewidth=0.5) +
        gghalves::geom_half_point_panel(aes(fill=mRNA), side = "l", shape=21, size=3, color="white") +
        ggpubr::stat_compare_means(method="wilcox.test", paired = T, label="p.signif", comparisons = my_comparisons) +
        ggpubr::stat_compare_means(method="anova", label.y=2.2) +
        scale_x_discrete(name=NULL) +
        scale_y_continuous(name="Data Fidelity Level") +
        scale_fill_manual(values=color) +
        ggtitle(title) +
        cowplot::theme_cowplot() +
        theme(panel.border = element_rect(color="black"), 
              strip.background = element_rect(fill=NA, color=NA), 
              legend.position = "bottom") +
        guides(fill = guide_legend(nrow=1, title="Compared_object"))
    }
    
    res <- list("figure"=p, "dat"=ggdat)
    
    return(res)
  }
  
  move_file <- function(file, dir) {
    stat <- file.copy(from = file, to = paste0(dir, "/"), overwrite = T, copy.date = T)
    if (stat) {
      stat <- file.remove(file)
    }
  }
}

# calculate the data fidelity
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir <- paste0("results/2.pic/44.cluster_dat_fidelity/", data_type)
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
        
        # calculate the ppm for each mRNA level
        if (T) {
          gpat_motif_cs <- gene_body_batch_ppm(gpat_cs, IDs)
          gpat_motif <- lapply(gpat_dat, gene_body_batch_ppm, IDs=IDs)
          names(gpat_motif)
        }
        
        # statistics the gpat pattern motif similarity
        if (T) {
          # raw chromatin states
          if (T) {
            cluster_type <- "CS"
            title <- paste0(data_type, ": ", cell, " ", cluster_type)
            p <- plot_cor_with_target(gpat_motif[[1]], title)
            
            file <- paste0("44.ppm_cor_in_all_", data_type, "s_Body.", cluster_type, ".", cell, ".", type, ".pdf")
            ggsave(filename = file, plot = p, width = 12, height = 8)
            move_file(file, dir)
          }
          
          # chromIDEAS
          if (T) {
            cluster_type <- "chromideas"
            title <- paste0(data_type, ": ", cell, " ", cluster_type)
            p <- plot_cor_with_target(gpat_motif[[2]], title)
            
            file <- paste0("44.ppm_cor_in_all_", data_type, "s_Body.", cluster_type, ".", cell, ".", type, ".pdf")
            ggsave(filename = file, plot = p, width = 12, height = 8)
            move_file(file, dir)
          }
          
          # kmeans
          if (T) {
            cluster_type <- "kmeans"
            title <- paste0(data_type, ": ", cell, " ", cluster_type)
            p <- plot_cor_with_target(gpat_motif[[3]], title)
            
            file <- paste0("44.ppm_cor_in_all_", data_type, "s_Body.", cluster_type, ".", cell, ".", type, ".pdf")
            ggsave(filename = file, plot = p, width = 12, height = 8)
            move_file(file, dir)
          }
          
          # ct1
          if (T) {
            cluster_type <- "ct1"
            title <- paste0(data_type, ": ", cell, " ", cluster_type)
            p <- plot_cor_with_target(gpat_motif[[4]], title)
            
            file <- paste0("44.ppm_cor_in_all_", data_type, "s_Body.", cluster_type, ".", cell, ".", type, ".pdf")
            ggsave(filename = file, plot = p, width = 12, height = 8)
            move_file(file, dir)
          }
          
          # ct2
          if (T) {
            cluster_type <- "ct2"
            title <- paste0(data_type, ": ", cell, " ", cluster_type)
            p <- plot_cor_with_target(gpat_motif[[5]], title)
            
            file <- paste0("44.ppm_cor_in_all_", data_type, "s_Body.", cluster_type, ".", cell, ".", type, ".pdf")
            ggsave(filename = file, plot = p, width = 12, height = 8)
            move_file(file, dir)
          }
          
          # ct3
          if (T) {
            cluster_type <- "ct3"
            title <- paste0(data_type, ": ", cell, " ", cluster_type)
            p <- plot_cor_with_target(gpat_motif[[6]], title)
            
            file <- paste0("44.ppm_cor_in_all_", data_type, "s_Body.", cluster_type, ".", cell, ".", type, ".pdf")
            ggsave(filename = file, plot = p, width = 12, height = 8)
            move_file(file, dir)
          }
          
          # ct4
          if (T) {
            cluster_type <- "ct4"
            title <- paste0(data_type, ": ", cell, " ", cluster_type)
            p <- plot_cor_with_target(gpat_motif[[7]], title)
            
            file <- paste0("44.ppm_cor_in_all_", data_type, "s_Body.", cluster_type, ".", cell, ".", type, ".pdf")
            ggsave(filename = file, plot = p, width = 12, height = 8)
            move_file(file, dir)
          }
          
          # ct5
          if (T) {
            cluster_type <- "ct5"
            title <- paste0(data_type, ": ", cell, " ", cluster_type)
            p <- plot_cor_with_target(gpat_motif[[8]], title)
            
            file <- paste0("44.ppm_cor_in_all_", data_type, "s_Body.", cluster_type, ".", cell, ".", type, ".pdf")
            ggsave(filename = file, plot = p, width = 12, height = 8)
            move_file(file, dir)
          }
          
          # summary
          if (T) {
            res <- plot_cor_similarity(gpat_motif_cs = gpat_motif[[1]], 
                                       gpat_motif_others = gpat_motif[2:8], 
                                       names = c("chromIDEAS", "kmeans", paste0("CT", 1:n_ct)), 
                                       title = paste0(data_type, ": ", cell, " Correlation Similarity Summary"), 
                                       type="barplot", 
                                       method="pearson")
            
            file <- paste0("44.ppm_cor_in_all_", data_type, "s_Body.summary.", cell, ".", type, ".pdf")
            ggsave(filename = file, plot = res[["figure"]], width = 12, height = 8)
            move_file(file, dir)
            
            file <- paste0("results/1.tab/44.ppm_cor_in_all_", data_type, "s_Body.summary.", cell, ".", type, ".csv")
            write.table(res[["dat"]], file = file, quote = F, sep = ",", col.names = T, row.names = F)
          }
          
          rm(title, cluster_type, res)
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
        
        # calculate the ppm for each mRNA level
        if (T) {
          gpat_motif_cs <- gene_body_batch_ppm(gpat_cs, IDs)
          gpat_motif <- lapply(gpat_dat, gene_body_batch_ppm, IDs=IDs)
          names(gpat_motif)
        }
        
        # statistics the gpat pattern motif similarity
        if (T) {
          # raw chromatin states
          if (T) {
            cluster_type <- "CS"
            title <- paste0(data_type, ": merged ", cluster_type)
            p <- plot_cor_with_target(gpat_motif[[1]], title)
            
            file <- paste0("44.ppm_cor_in_all_", data_type, "s_Body.", cluster_type, ".merged.", type, ".pdf")
            ggsave(filename = file, plot = p, width = 12, height = 8)
            move_file(file, dir)
          }
          
          # chromIDEAS
          if (T) {
            cluster_type <- "chromideas"
            title <- paste0(data_type, ": merged ", cluster_type)
            p <- plot_cor_with_target(gpat_motif[[2]], title)
            
            file <- paste0("44.ppm_cor_in_all_", data_type, "s_Body.", cluster_type, ".merged.", type, ".pdf")
            ggsave(filename = file, plot = p, width = 12, height = 8)
            move_file(file, dir)
          }
          
          # kmeans
          if (T) {
            cluster_type <- "kmeans"
            title <- paste0(data_type, ": merged ", cluster_type)
            p <- plot_cor_with_target(gpat_motif[[3]], title)
            
            file <- paste0("44.ppm_cor_in_all_", data_type, "s_Body.", cluster_type, ".merged.", type, ".pdf")
            ggsave(filename = file, plot = p, width = 12, height = 8)
            move_file(file, dir)
          }
          
          # ct1
          if (T) {
            cluster_type <- "ct1"
            title <- paste0(data_type, ": merged ", cluster_type)
            p <- plot_cor_with_target(gpat_motif[[4]], title)
            
            file <- paste0("44.ppm_cor_in_all_", data_type, "s_Body.", cluster_type, ".merged.", type, ".pdf")
            ggsave(filename = file, plot = p, width = 12, height = 8)
            move_file(file, dir)
          }
          
          # ct2
          if (T) {
            cluster_type <- "ct2"
            title <- paste0(data_type, ": merged ", cluster_type)
            p <- plot_cor_with_target(gpat_motif[[5]], title)
            
            file <- paste0("44.ppm_cor_in_all_", data_type, "s_Body.", cluster_type, ".merged.", type, ".pdf")
            ggsave(filename = file, plot = p, width = 12, height = 8)
            move_file(file, dir)
          }
          
          # ct3
          if (T) {
            cluster_type <- "ct3"
            title <- paste0(data_type, ": merged ", cluster_type)
            p <- plot_cor_with_target(gpat_motif[[6]], title)
            
            file <- paste0("44.ppm_cor_in_all_", data_type, "s_Body.", cluster_type, ".merged.", type, ".pdf")
            ggsave(filename = file, plot = p, width = 12, height = 8)
            move_file(file, dir)
          }
          
          # ct4
          if (T) {
            cluster_type <- "ct4"
            title <- paste0(data_type, ": merged ", cluster_type)
            p <- plot_cor_with_target(gpat_motif[[7]], title)
            
            file <- paste0("44.ppm_cor_in_all_", data_type, "s_Body.", cluster_type, ".merged.", type, ".pdf")
            ggsave(filename = file, plot = p, width = 12, height = 8)
            move_file(file, dir)
          }
          
          # ct5
          if (T) {
            cluster_type <- "ct5"
            title <- paste0(data_type, ": merged ", cluster_type)
            p <- plot_cor_with_target(gpat_motif[[8]], title)
            
            file <- paste0("44.ppm_cor_in_all_", data_type, "s_Body.", cluster_type, ".merged.", type, ".pdf")
            ggsave(filename = file, plot = p, width = 12, height = 8)
            move_file(file, dir)
          }
          
          # summary
          if (T) {
            res <- plot_cor_similarity(gpat_motif_cs = gpat_motif[[1]], 
                                       gpat_motif_others = gpat_motif[2:8], 
                                       names = c("chromIDEAS", "kmeans", paste0("CT", 1:n_ct)), 
                                       title = paste0(data_type, ": merged Correlation Similarity Summary"), 
                                       type="barplot", 
                                       method="pearson")
            
            file <- paste0("44.ppm_cor_in_all_", data_type, "s_Body.summary.merged.", type, ".pdf")
            ggsave(filename = file, plot = res[["figure"]], width = 12, height = 8)
            move_file(file, dir)
            
            file <- paste0("results/1.tab/44.ppm_cor_in_all_", data_type, "s_Body.summary.merged.", type, ".csv")
            write.table(res[["dat"]], file = file, quote = F, sep = ",", col.names = T, row.names = F)
          }
          
          rm(title, cluster_type, res)
        }
      }
    }
  }
}

## tx: 
##     merged: 
##         thp1: 
##                   Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
## chromIDEAS  0.44457276 0.6288829 0.7921549 0.7445477 0.8627114 0.9481134
## kmeans      0.64296501 0.7927084 0.8653699 0.8493964 0.9188587 0.9486708
## CT1         0.21538762 0.4976982 0.5898270 0.5811678 0.6997604 0.8505953
## CT2         0.31645968 0.6652507 0.7265272 0.7493558 0.9103782 0.9432345
## CT3         0.26384850 0.5504452 0.7757099 0.7156581 0.8841129 0.9438957
## CT4         0.47549319 0.7175853 0.8418604 0.8068588 0.9241596 0.9716859
## CT5        -0.03612412 0.3797089 0.5953043 0.5432563 0.7388598 0.9554623
##         cd34: 
##                 Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
## chromIDEAS 0.1262339 0.4529419 0.6714030 0.6527256 0.8325187 0.9328481
## kmeans     0.5707462 0.7776881 0.8429533 0.8449903 0.9157654 0.9761501
## CT1        0.2424244 0.5448191 0.6280259 0.6091782 0.6904326 0.8832605
## CT2        0.1376623 0.5183038 0.6401977 0.6588992 0.8548577 0.9437822
## CT3        0.2196073 0.7118153 0.8238975 0.7921050 0.9329800 0.9590919
## CT4        0.3311033 0.6294232 0.8012470 0.7560099 0.8964378 0.9780867
## CT5        0.1715529 0.6208734 0.7046759 0.7075330 0.8135408 0.9411212
##         merged: 
##                   Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
## chromIDEAS  0.32878824 0.5940473 0.7422734 0.7097691 0.8556487 0.9324792
## kmeans      0.69211010 0.8051993 0.8537958 0.8542174 0.9043336 0.9592061
## CT1         0.24080880 0.5366252 0.6162329 0.6009202 0.7055073 0.9227042
## CT2         0.22764268 0.6558108 0.7023309 0.7386204 0.9036892 0.9446441
## CT3         0.27490993 0.6406875 0.7884412 0.7452489 0.9079720 0.9432336
## CT4         0.44554942 0.6817128 0.8387880 0.7935267 0.9190061 0.9769859
## CT5        -0.08000722 0.4992816 0.6655295 0.6264194 0.7915105 0.9372893
##     single: 
##         thp1: 
##                   Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
## chromIDEAS 0.426889039 0.6360750 0.7879271 0.7362470 0.8684816 0.9420208
## kmeans     0.635039579 0.7792079 0.8569980 0.8403179 0.9190404 0.9479722
## CT1        0.348353192 0.5944706 0.6589047 0.6413135 0.6958905 0.8002764
## CT2        0.228108565 0.6050411 0.6846919 0.7008925 0.8426474 0.9180830
## CT3        0.319863384 0.5199319 0.6711514 0.6531613 0.7832884 0.9172157
## CT4        0.543827158 0.7220376 0.8242535 0.8049616 0.8972492 0.9689109
## CT5        0.002525975 0.4226718 0.6379635 0.5762154 0.7583293 0.9369608
##         cd34: 
##                   Min.    1st Qu.     Median        Mean   3rd Qu.      Max.
## chromIDEAS  0.51793491  0.7351761  0.8162522  0.79792555 0.8872047 0.9604046
## kmeans      0.57074617  0.7776881  0.8429533  0.84499029 0.9157654 0.9761501
## CT1        -0.66678519 -0.4634805 -0.2196425 -0.08329825 0.3289247 0.7867167
## CT2         0.10228245  0.6014681  0.6772209  0.68597427 0.8409425 0.9218225
## CT3         0.11759986  0.5237834  0.6218481  0.66517026 0.9368422 0.9740603
## CT4         0.41781627  0.6338447  0.7702548  0.74638060 0.8573933 0.9603536
## CT5         0.04989744  0.6322012  0.7266665  0.71326100 0.8338832 0.9351467