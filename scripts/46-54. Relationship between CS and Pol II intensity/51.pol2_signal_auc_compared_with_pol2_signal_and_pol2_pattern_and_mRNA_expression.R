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
pol2_levels <- 10

# define functions
if (T) {
  # common value
  if (T) {
    name <- c(
      paste0("U", up_bin_num:1), 
      "TSS", 
      paste0("B", 1:body_bin_num), 
      "TES", 
      paste0("D", 1:down_bin_num)
    )
  }
  
  # auc as signal for each gene (modify based on script38)
  gene_pol2_signal <- function(mat) {
    auc <- apply(mat, 1, function(row) {
      # row <- unlist(mat[1, ])
      value <- as.numeric(row[name])
      
      a <- value[-1]
      b <- value[-length(value)]
      
      sum((a+b)*1/2)
    })
    auc <- data.frame(id = mat$gene_id, 
                      signal = as.numeric(auc))
    
    return(auc)
  }
  
  # data prepare (modify based on script36)
  quantitle_pol2 <- function(mat, pol2_levels) {
    # get subdat0 and Q1-Q10
    if (T) {
      for (q in seq(1, pol2_levels)) {
        min <- quantile(mat$signal, (q-1)/10)
        max <- quantile(mat$signal, q/10)
        
        if (max == max(mat$signal)) {
          mat[mat$signal>=min & mat$signal<=max, "Quantile"] <- paste0("Q", q)
        }
        if (max < max(mat$signal)) {
          mat[mat$signal>=min & mat$signal<max, "Quantile"] <- paste0("Q", q)
        }
      }
    }
    
    return(mat)
  }
  
  move_file <- function(file, dir) {
    stat <- file.copy(from = file, to = paste0(dir, "/"), overwrite = T, copy.date = T)
    if (stat) {
      stat <- file.remove(file)
    }
  }
}

# calculate the pol2 auc value
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir1 <- paste0("data/saved_data/51.pol2_signal/", data_type)
      if (! dir.exists(dir1)) {
        dir.create(dir1, showWarnings = F, recursive = T)
      }
      
      dir2 <- paste0("results/2.pic/51.pol2_signal/", data_type)
      if (! dir.exists(dir2)) {
        dir.create(dir2, showWarnings = F, recursive = T)
      }
    }
    
    for (type in c("merged", "single")) {
      cat(paste0("\t", type, ": \n"))
      
      # cell specific data
      if (type == "single") {
        for (cell in c(cell1, cell2)) {
          cat(paste0("\t\t", cell, ": \n"))
          
          # test data
          if (F) {
            data_type <- "tx"
            type <- "single"
            cell <- cell1
          }
          
          # calculate the pol2 auc for all genes
          if (T) {
            file <- paste0(dir1, "/51.pol2_signal_auc_Q", pol2_levels, "_in_", data_type, "s_profile_Body.", cell, ".", type, ".qs")
            
            if (file.exists(file)) {
              mat <- qread(file, nthreads = 6)
            }
            if (! file.exists(file)) {
              # get pol2 matrix: mat
              if (T) {
                cell_dat <- function(cell, type) {
                  input <- paste0("data/saved_data/48.mk_signal_distribution/tx/48.", data_type, "_Body_", up_bin_num, "UP_", 
                                  down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".pol2.mean_mat.", cell, ".", type, ".qs")
                  
                  mat <- qread(input, nthreads = 6)
                  rownames(mat) <- paste0(mat$gene_id, "@", cell)
                  
                  return(mat)
                }
                mat <- cell_dat(cell, type)
              }
              
              # calculate the auc value
              if (T) {
                signal <- gene_pol2_signal(mat)
              }
              
              # merge the data
              if (T) {
                mat <- cbind(mat, 
                             signal = signal$signal)
                
                rm(signal)
              }
              
              # get pol2 signal quantile levels
              if (T) {
                mat <- quantitle_pol2(mat, pol2_levels)
              }
              
              qsave(mat, file, nthreads = 6)
            }
          }
          
          k <- 4
          
          # get pol2 cluster data
          if (T) {
            input <- paste0("data/saved_data/50.pol2_pattern/", data_type, "/50.pol2_pat_bs", bin_size, 
                            "_kmeans", k, "_", data_type, "_", cell, ".", type, ".qs")
            cluster <- qread(input)
            cluster <- cluster[[2]]
          }
          
          # add cluster info into mat
          if (T) {
            mat <- mat[cluster$ID, ]
            mat$cluster <- cluster$cluster
          }
          
          # save the data
          if (T) {
            file <- paste0(dir1, "/51.pol2_signal_auc_Q", pol2_levels, "_in_", data_type, "s_format.keamns", k, ".", cell, ".", type, ".qs")
            
            if (! file.exists(file)) {
              mat <- mat[, c("gene_id", "cluster", "signal", "Quantile")]
              
              qsave(mat, file = file, nthreads = 6)
            }
          }
        }
      }
      
      # merged data
      if (type == "merged") {
        # test data
        if (F) {
          data_type <- "tx"
          type <- "merged"
        }
        
        # calculate the pol2 auc for all genes
        if (T) {
          file <- paste0(dir1, "/51.pol2_signal_auc_Q", pol2_levels, "_in_", data_type, "s_profile_Body.merged.", type, ".qs")
          
          if (file.exists(file)) {
            mat <- qread(file, nthreads = 6)
          }
          if (! file.exists(file)) {
            # get pol2 matrix: mat
            if (T) {
              cell_dat <- function(cell, type) {
                input <- paste0("data/saved_data/48.mk_signal_distribution/tx/48.", data_type, "_Body_", up_bin_num, "UP_", 
                                down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".pol2.mean_mat.", cell, ".", type, ".qs")
                
                mat <- qread(input, nthreads = 6)
                rownames(mat) <- paste0(mat$gene_id, "@", cell)
                
                return(mat)
              }
              dat1 <- cell_dat(cell1, type="single")
              dat2 <- cell_dat(cell2, type="single")
              
              mat <- rbind(dat1, dat2)
              
              rm(dat1, dat2)
            }
            
            # calculate the auc value
            if (T) {
              signal <- gene_pol2_signal(mat)
            }
            
            # merge the data
            if (T) {
              mat <- cbind(mat, 
                           signal = signal$signal)
              
              rm(signal)
            }
            
            # get pol2 signal quantile levels
            if (T) {
              mat <- quantitle_pol2(mat, pol2_levels)
            }
            
            qsave(mat, file, nthreads = 6)
          }
        }
        
        k <- 4
        
        # get pol2 cluster data
        if (T) {
          input <- paste0("data/saved_data/50.pol2_pattern/", data_type, "/50.pol2_pat_bs", bin_size, 
                          "_kmeans", k, "_", data_type, "_merged.", type, ".qs")
          cluster <- qread(input)
          cluster <- cluster[[2]]
        }
        
        # add cluster info into mat
        if (T) {
          mat <- mat[cluster$ID, ]
          mat$cluster <- cluster$cluster
        }
        
        # save the data
        if (T) {
          file <- paste0(dir1, "/51.pol2_signal_auc_Q", pol2_levels, "_in_", data_type, "s_format.keamns", k, ".merged.", type, ".qs")
          
          if (! file.exists(file)) {
            mat <- mat[, c("gene_id", "cluster", "signal", "Quantile")]
            
            qsave(mat, file = file, nthreads = 6)
          }
        }
      }
    }
  }
}

# plot pol2 auc figure
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir1 <- paste0("data/saved_data/51.pol2_signal/", data_type)
      if (! dir.exists(dir1)) {
        dir.create(dir1, showWarnings = F, recursive = T)
      }
      
      dir2 <- paste0("results/2.pic/51.pol2_signal/", data_type)
      if (! dir.exists(dir2)) {
        dir.create(dir2, showWarnings = F, recursive = T)
      }
    }
    
    for (type in c("merged", "single")) {
      cat(paste0("\t", type, ": \n"))
      
      # cell specific data
      if (type == "single") {
        for (cell in c(cell1, cell2)) {
          cat(paste0("\t\t", cell, ": \n"))
          
          # test data
          if (F) {
            data_type <- "tx"
            type <- "single"
            cell <- cell1
          }
          
          k <- 4
          
          # read pol2 auc of genes
          if (T) {
            file <- paste0(dir1, "/51.pol2_signal_auc_Q", pol2_levels, "_in_", data_type, "s_format.keamns", k, ".", cell, ".", type, ".qs")
            
            mat <- qread(file, nthreads = 6)
          }
          
          # format the data
          if (T) {
            summary(mat$signal)
            mat$log2signal <- log2(mat$signal)
            
            mat$Quantile <- factor(mat$Quantile, levels = paste0("Q", 1:pol2_levels))
            mat$cluster <- factor(paste0("pol2_", mat$cluster), levels = paste0("pol2_", 1:k))
            
            summary(mat$log2signal)
          }
          
          # ggplot2: pol2 signal vs pol2 pat
          if (T) {
            head(mat)
            
            p <- ggplot(mat, aes(x=cluster, y=log2signal)) +
              geom_violin(aes(group=cluster, fill=cluster), scale = "width") +
              geom_boxplot(aes(group=cluster), fill=NA, width=0.2) +
              xlab(NULL) +
              ylab("log2(Pol2 signal AUC)") +
              ggtitle(paste0(cell, "+", type)) +
              cowplot::theme_cowplot() +
              theme(axis.title = element_text(size = rel(1.2)),
                    axis.text = element_text(size = rel(1)),
                    legend.position = "none",
                    panel.background = element_rect(fill = NA, color="black"), 
                    strip.background = element_rect(fill = NA, color="black"), 
                    strip.text = element_text(size = rel(1)))
            
            file <- paste0("51.pol2_auc_distribution_vs_pol2_", k, "pat_in_", data_type, "s_.", cell, ".", type, ".pdf")
            ggsave(filename = file, plot = p, width = 12, height = 8)
            move_file(file, dir2)
          }
          
          # ggplot2: number of each pol2 pat group within pol2 signal quantile level
          if (T) {
            # statistics raw data
            if (T) {
              group <- as.data.frame.table(table(paste0(mat$Quantile, "@", mat$cluster)))
            }
            
            # format the info
            if (T) {
              head(group)
              
              group$quantile <- str_split(group$Var1, "@", simplify = T)[, 1]
              group$cluster <- str_split(group$Var1, "@", simplify = T)[, 2]
              
              group$quantile <- factor(group$quantile, levels = paste0("Q", 1:pol2_levels))
              group$cluster <- factor(group$cluster, levels = paste0("pol2_", 1:k))
            }
            
            # statistics the gene percentage of each cluster
            if (T) {
              head(group)
              
              for (q in paste0("Q", 1:pol2_levels)) {
                # q <- "Q1"
                torf <- group$quantile == q
                value <- group[torf, "Freq"]
                
                group[torf, "Percentage"] <- value/sum(value)*100
                
                rm(torf, value, q)
              }
              
              group$label <- paste0(round(group$Percentage, 2), "%")
            }
            
            # percentage ggplot
            if (T) {
              head(group)
              
              p <- ggplot(group) +
                geom_bar(aes(x=quantile, y=Percentage, fill=cluster), stat="identity", position="stack", color="black") +
                scale_x_discrete(name = "Pol2 signal level") +
                scale_y_continuous(breaks = seq(0, 100, 25)) +
                ylab("Gene Percentage") +
                ggtitle(paste0(cell, ": ", type)) +
                cowplot::theme_cowplot() +
                theme(panel.border = element_rect(color="black"), 
                      strip.background = element_rect(fill=NA, color=NA), 
                      legend.position = "right")
            }
            
            file <- paste0("51.number_of_pol2_", k, "pat_within_pol2_auc_levels_in_", data_type, "s.", cell, ".", type, ".pdf")
            ggsave(filename = file, plot = p, width = 12, height = 8)
            move_file(file, dir2)
            
            rm(p, file, group)
          }
        }
      }
      
      # merged data
      if (type == "merged") {
        # test data
        if (F) {
          data_type <- "tx"
          type <- "merged"
        }
        
        k <- 4
        
        # read pol2 auc of genes
        if (T) {
          file <- paste0(dir1, "/51.pol2_signal_auc_Q", pol2_levels, "_in_", data_type, "s_format.keamns", k, ".merged.", type, ".qs")
          
          mat <- qread(file, nthreads = 6)
        }
        
        # format the data
        if (T) {
          summary(mat$signal)
          mat$log2signal <- log2(mat$signal)
          
          mat$Quantile <- factor(mat$Quantile, levels = paste0("Q", 1:pol2_levels))
          mat$cluster <- factor(paste0("pol2_", mat$cluster), levels = paste0("pol2_", 1:k))
          
          summary(mat$log2signal)
        }
        
        # ggplot2: pol2 signal vs pol2 pat
        if (T) {
          head(mat)
          
          p <- ggplot(mat, aes(x=cluster, y=log2signal)) +
            geom_violin(aes(group=cluster, fill=cluster), scale = "width") +
            geom_boxplot(aes(group=cluster), fill=NA, width=0.2) +
            xlab(NULL) +
            ylab("log2(Pol2 signal AUC)") +
            ggtitle(paste0("merged+", type)) +
            cowplot::theme_cowplot() +
            theme(axis.title = element_text(size = rel(1.2)),
                  axis.text = element_text(size = rel(1)),
                  legend.position = "none",
                  panel.background = element_rect(fill = NA, color="black"), 
                  strip.background = element_rect(fill = NA, color="black"), 
                  strip.text = element_text(size = rel(1)))
          
          file <- paste0("51.pol2_auc_distribution_vs_pol2_", k, "pat_in_", data_type, "s_.merged.", type, ".pdf")
          ggsave(filename = file, plot = p, width = 12, height = 8)
          move_file(file, dir2)
        }
        
        # ggplot2: number of each pol2 pat group within pol2 signal quantile level
        if (T) {
          # statistics raw data
          if (T) {
            group <- as.data.frame.table(table(paste0(mat$Quantile, "@", mat$cluster)))
          }
          
          # format the info
          if (T) {
            head(group)
            
            group$quantile <- str_split(group$Var1, "@", simplify = T)[, 1]
            group$cluster <- str_split(group$Var1, "@", simplify = T)[, 2]
            
            group$quantile <- factor(group$quantile, levels = paste0("Q", 1:pol2_levels))
            group$cluster <- factor(group$cluster, levels = paste0("pol2_", 1:k))
          }
          
          # statistics the gene percentage of each cluster
          if (T) {
            head(group)
            
            for (q in paste0("Q", 1:pol2_levels)) {
              # q <- "Q1"
              torf <- group$quantile == q
              value <- group[torf, "Freq"]
              
              group[torf, "Percentage"] <- value/sum(value)*100
              
              rm(torf, value, q)
            }
            
            group$label <- paste0(round(group$Percentage, 2), "%")
          }
          
          # percentage ggplot
          if (T) {
            head(group)
            
            p <- ggplot(group) +
              geom_bar(aes(x=quantile, y=Percentage, fill=cluster), stat="identity", position="stack", color="black") +
              scale_x_discrete(name = "Pol2 signal level") +
              scale_y_continuous(breaks = seq(0, 100, 25)) +
              ylab("Gene Percentage") +
              ggtitle(paste0("merged: ", type)) +
              cowplot::theme_cowplot() +
              theme(panel.border = element_rect(color="black"), 
                    strip.background = element_rect(fill=NA, color=NA), 
                    legend.position = "right")
          }
          
          file <- paste0("51.number_of_pol2_", k, "pat_within_pol2_auc_levels_in_", data_type, "s.merged.", type, ".pdf")
          ggsave(filename = file, plot = p, width = 12, height = 8)
          move_file(file, dir2)
          
          rm(p, file, group)
        }
      }
    }
  }
}

# the relationship between pol2 auc and mRNA level
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir1 <- paste0("data/saved_data/51.pol2_signal/", data_type)
      if (! dir.exists(dir1)) {
        dir.create(dir1, showWarnings = F, recursive = T)
      }
      
      dir2 <- paste0("results/2.pic/51.pol2_signal/", data_type)
      if (! dir.exists(dir2)) {
        dir.create(dir2, showWarnings = F, recursive = T)
      }
    }
    
    for (type in c("merged", "single")) {
      cat(paste0("\t", type, ": \n"))
      
      # cell specific data
      if (type == "single") {
        for (cell in c(cell1, cell2)) {
          cat(paste0("\t\t", cell, ": \n"))
          
          # test data
          if (F) {
            data_type <- "tx"
            type <- "single"
            cell <- cell1
          }
          
          k <- 4
          
          # read pol2 auc of genes
          if (T) {
            file <- paste0(dir1, "/51.pol2_signal_auc_Q", pol2_levels, "_in_", data_type, "s_format.keamns", k, ".", cell, ".", type, ".qs")
            
            mat <- qread(file, nthreads = 6)
          }
          
          # format the data
          if (T) {
            summary(mat$signal)
            mat$log2signal <- log2(mat$signal)
            
            mat$Quantile <- factor(mat$Quantile, levels = paste0("Q", 1:pol2_levels))
            mat$cluster <- factor(paste0("pol2_", mat$cluster), levels = paste0("pol2_", 1:k))
            
            summary(mat$log2signal)
          }
          
          # get mRNA level data
          if (T) {
            file <- ifelse(data_type == "gene", 
                           paste0("data/saved_data/6.4_type_gene_expression_level_classification_", cell, ".qs"), 
                           paste0("data/saved_data/7.4_type_tx_expression_level_classification_", cell, ".qs"))
            
            cell_mrna <- qread(file, nthreads = 6)
            cell_mrna <- cell_mrna[[1]]
            
            cell_mrna <- cell_mrna[names(cell_mrna) %in% mat$gene_id]
          }
          
          # add mRNA into mat
          if (T) {
            mat$mRNA <- cell_mrna[match(mat$gene_id, names(cell_mrna))]
            
            non0min <- min(cell_mrna[cell_mrna>0])
            mat$log2mRNA <- log2(mat$mRNA+non0min)
            
            rm(cell_mrna, non0min)
          }
          
          # color setting
          if (T) {
            col_fun <- circlize::colorRamp2(c(1, pol2_levels), c("blue", "red"))
          }
          
          # statistics group
          if (T) {
            max <- max(mat$log2mRNA[mat$Quantile=="Q5"])
            top <- max*2
            
            my_comparisons_signal <- lapply(2:pol2_levels, function(b) {
              # b <- 2
              a <- (b-1)
              
              paste0("Q", a:b)
            })
          }
          
          # pol2 sig auc vs mRNA (not split by pol2 pat)
          if (T) {
            p <- ggplot(mat, aes(x=Quantile, y=log2mRNA)) +
              geom_violin(aes(group=Quantile, fill=Quantile), scale = "width") +
              geom_boxplot(aes(group=Quantile), fill=NA, width=0.2) +
              scale_fill_manual(values = col_fun(1:pol2_levels),breaks = paste0("Q", 1:pol2_levels)) +
              ggpubr::stat_compare_means(method="t.test", label="p.signif", 
                                         comparisons = my_comparisons_signal, vjust=0.5, label.y=seq(max, top, length.out=pol2_levels-1)) +
              xlab("Pol2 Signal Strength") +
              scale_y_continuous(name="log2(RPKM)") +
              ggtitle(paste0(cell, "+", type)) +
              cowplot::theme_cowplot() +
              theme(axis.title = element_text(size = rel(1.2)),
                    axis.text = element_text(size = rel(1)),
                    legend.position = "none",
                    panel.background = element_rect(fill = NA, color="black"), 
                    strip.background = element_rect(fill = NA, color="black"), 
                    strip.text = element_text(size = rel(1)))
            
            file <- paste0("51.pol2_auc_vs_mRNA_in_", data_type, "s.", cell, ".", type, ".pdf")
            ggsave(filename = file, plot = p, width = 12, height = 8)
            
            move_file(file, dir2)
          }
          
          # pol2 sig auc vs mRNA (split by pol2 pat)
          if (T) {
            p <- p + facet_grid(cluster~.)
            
            file <- paste0("51.pol2_auc_vs_mRNA_under_pol2_", k, "pat_in_", data_type, "s.", cell, ".", type, ".pdf")
            ggsave(filename = file, plot = p, width = 12, height = 8)
            
            move_file(file, dir2)
          }
        }
      }
      
      # merged data
      if (type == "merged") {
        # test data
        if (F) {
          data_type <- "tx"
          type <- "merged"
        }
        
        k <- 4
        
        # read pol2 auc of genes
        if (T) {
          file <- paste0(dir1, "/51.pol2_signal_auc_Q", pol2_levels, "_in_", data_type, "s_format.keamns", k, ".merged.", type, ".qs")
          
          mat <- qread(file, nthreads = 6)
        }
        
        # format the data
        if (T) {
          summary(mat$signal)
          mat$log2signal <- log2(mat$signal)
          
          mat$Quantile <- factor(mat$Quantile, levels = paste0("Q", 1:pol2_levels))
          mat$cluster <- factor(paste0("pol2_", mat$cluster), levels = paste0("pol2_", 1:k))
          
          summary(mat$log2signal)
        }
        
        # get mRNA level data
        if (T) {
          rna_Q <- function(cell) {
            # cell <- cell1
            file <- ifelse(data_type == "gene", 
                           paste0("data/saved_data/6.4_type_gene_expression_level_classification_", cell, ".qs"), 
                           paste0("data/saved_data/7.4_type_tx_expression_level_classification_", cell, ".qs"))
            
            cell_mrna <- qread(file, nthreads = 6)
            cell_mrna <- cell_mrna[[1]]
            
            names(cell_mrna) <- paste0(names(cell_mrna), "@", cell)
            cell_mrna <- cell_mrna[names(cell_mrna) %in% rownames(mat)]
            
            return(cell_mrna)
          }
          
          dat1 <- rna_Q(cell1)
          dat2 <- rna_Q(cell2)
          cell_mrna <- c(dat1, dat2)
          
          rm(dat1, dat2)
        }
        
        # add mRNA into mat
        if (T) {
          mat$mRNA <- cell_mrna[match(rownames(mat), names(cell_mrna))]
          
          non0min <- min(cell_mrna[cell_mrna>0])
          mat$log2mRNA <- log2(mat$mRNA+non0min)
          
          rm(cell_mrna, non0min)
        }
        
        # color setting
        if (T) {
          col_fun <- circlize::colorRamp2(c(1, pol2_levels), c("blue", "red"))
        }
        
        # statistics group
        if (T) {
          max <- max(mat$log2mRNA[mat$Quantile=="Q5"])
          top <- max*2
          
          my_comparisons_signal <- lapply(2:pol2_levels, function(b) {
            # b <- 2
            a <- (b-1)
            
            paste0("Q", a:b)
          })
        }
        
        # pol2 sig auc vs mRNA (not split by pol2 pat)
        if (T) {
          p <- ggplot(mat, aes(x=Quantile, y=log2mRNA)) +
            geom_violin(aes(group=Quantile, fill=Quantile), scale = "width") +
            geom_boxplot(aes(group=Quantile), fill=NA, width=0.2) +
            scale_fill_manual(values = col_fun(1:pol2_levels),breaks = paste0("Q", 1:pol2_levels)) +
            ggpubr::stat_compare_means(method="t.test", label="p.signif", 
                                       comparisons = my_comparisons_signal, vjust=0.5, label.y=seq(max, top, length.out=pol2_levels-1)) +
            xlab("Pol2 Signal Strength") +
            scale_y_continuous(name="log2(RPKM)") +
            ggtitle(paste0("merged+", type)) +
            cowplot::theme_cowplot() +
            theme(axis.title = element_text(size = rel(1.2)),
                  axis.text = element_text(size = rel(1)),
                  legend.position = "none",
                  panel.background = element_rect(fill = NA, color="black"), 
                  strip.background = element_rect(fill = NA, color="black"), 
                  strip.text = element_text(size = rel(1)))
          
          file <- paste0("51.pol2_auc_vs_mRNA_in_", data_type, "s.merged.", type, ".pdf")
          ggsave(filename = file, plot = p, width = 12, height = 8)
          
          move_file(file, dir2)
        }
        
        # pol2 sig auc vs mRNA (split by pol2 pat)
        if (T) {
          p <- p + facet_grid(cluster~.)
          
          file <- paste0("51.pol2_auc_vs_mRNA_under_pol2_", k, "pat_in_", data_type, "s.merged.", type, ".pdf")
          ggsave(filename = file, plot = p, width = 12, height = 8)
          
          move_file(file, dir2)
        }
      }
    }
  }
}


