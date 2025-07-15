# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(stringr)
  library(reshape2)
  library(qs)
  suppressPackageStartupMessages(library(circlize))
  suppressPackageStartupMessages(library(ComplexHeatmap))
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
    order <- c(
      paste0("U", up_bin_num:1), 
      "TSS", 
      paste0("B", 1:body_bin_num), 
      "TES", 
      paste0("D", 1:down_bin_num)
    )
    
    labels <- c(
      paste0("U", up_bin_num), 
      "TSS", 
      "TES", 
      paste0("D", down_bin_num)
    )
  }
  
  # gpat pattern seq stat (script43)
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
  
  # auc as signal for each gene (modify based on script51)
  state_auc <- function(gpat_ppm, order) {
    # format the data
    if (T) {
      auc_dat <- gpat_ppm[, order]
    }
    
    # calculate the auc value
    if (T) {
      auc <- apply(auc_dat, 1, function(row) {
        # row <- unlist(auc_dat[1, ])
        value <- as.numeric(row[order])
        
        a <- value[-1]
        b <- value[-length(value)]
        
        sum((a+b)*1/2)
      })
    }
    
    return(auc)
  }
  
  # data prepare (script51)
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
  
  # color setting (script43)
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
  
  # state percentage across pol2 signal
  state_percentage_across_pol2_sig <- function(ggdat, title) {
    # test data
    if (F) {
      title <- "title"
    }
    
    # color setting
    if (T) {
      col_fun <- colorRamp2(c(1, pol2_levels), c("blue", "red"))
    }
    
    # ggplot2
    if (T) {
      p <- ggplot(ggdat, aes(x=Loc, y=Percentage)) +
        geom_line(aes(group=group, color=pol2_sig), linewidth=1, alpha=0.6) +
        scale_colour_manual(values = col_fun(1:pol2_levels),
                            breaks = paste0("Q", 1:pol2_levels)) +
        scale_x_discrete(name = NULL, breaks=labels) +
        ggtitle(title) +
        facet_wrap(~state, ncol = 5, scales = "free") +
        cowplot::theme_cowplot() +
        theme(legend.position = c(0.90, 0.05), 
              panel.border = element_rect(color="black"), 
              strip.background = element_rect(color=NA, fill=NA))
    }
    
    return(p)
  }
  
  # plot state auc tendency across pol2 signal
  plot_state_auc_tendency <- function(gpat_ppm, title) {
    # test data
    if (F) {
      title <- "title"
    }
    
    p <- ggplot(gpat_ppm, aes(x=pol2_sig, y=auc)) +
      geom_line(aes(group=state, color=state), linewidth=1) +
      geom_point(size=1, alpha=0.8) +
      ggtitle(title) +
      xlab(NULL) +
      ylab("AUC (Area Under Distribution Curve)") +
      facet_wrap(~state, scale = "free", ncol = 5) +
      cowplot::theme_cowplot() +
      theme(legend.position = "none", 
            panel.border = element_rect(color="black"), 
            strip.background = element_rect(fill=NA, color=NA))
    
    return(p)
  }
}

# calculate the state auc value under pol2 sig genes
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir1 <- paste0("data/saved_data/52.pol2_ppm/", data_type)
      if (! dir.exists(dir1)) {
        dir.create(dir1, showWarnings = F, recursive = T)
      }
      
      dir2 <- paste0("results/2.pic/52.pol2_ppm/", data_type)
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
          
          # load the gpat data
          if (T) {
            file <- paste0("data/saved_data/42.gene_state_pattern/all_txs/42.", data_type, "_level_profile_Body_", up_bin_num, 
                           "UP_", down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".state.", cell, ".qs")
            gpat_cs <- qread(file)
          }
          
          k <- 4
          
          # get pol2 gene cluster: pol2cluster
          if (T) {
            file <- paste0("data/saved_data/51.pol2_signal/", data_type, "/51.pol2_signal_auc_Q", pol2_levels, 
                           "_in_", data_type, "s_format.keamns", k, ".", cell, ".", type, ".qs")
            
            pol2cluster <- qread(file)
          }
          
          # merge the data
          if (T) {
            over <- intersect(gpat_cs$gene_id, pol2cluster$gene_id)
            
            pol2cluster <- pol2cluster[match(over, pol2cluster$gene_id), ]
            gpat_cs <- gpat_cs[match(over, gpat_cs$gene_id), ]
            gpat_cs <- merge(gpat_cs, pol2cluster, by="gene_id")
            
            rm(over, pol2cluster)
          }
          
          # state auc: pol2 sig only
          if (T) {
            # calculate the ppm for each pol2 sig
            if (T) {
              gpat_ppm <- gene_body_batch_ppm(gpat_cs, split(gpat_cs$gene_id, gpat_cs$Quantile))
              colnames(gpat_ppm)[colnames(gpat_ppm) == "mRNA"] <- "Quantile"
            }
            
            # calculate the auc value
            if (T) {
              file <- paste0(dir1, "/52.cs_ppm_pol2_Q", pol2_levels, "_AUC_in_", data_type, "s.", cell, ".", type, ".qs")
              if (! file.exists(file)) {
                gpat_ppm$auc <- state_auc(gpat_ppm, order)
                
                qsave(gpat_ppm, file = file, nthreads = 6)
              }
            }
          }
          
          # state auc: pol2 sig + pol2 pat
          if (T) {
            # calculate the ppm for each pol2 sig and each pol2 pat
            if (T) {
              gpat_ppm <- gene_body_batch_ppm(gpat_cs, split(gpat_cs$gene_id, paste0(gpat_cs$Quantile, "@", gpat_cs$cluster)))
              colnames(gpat_ppm)[colnames(gpat_ppm) == "mRNA"] <- "ID"
              
              gpat_ppm$Quantile <- str_split(gpat_ppm$ID, "@", simplify = T)[, 1]
              gpat_ppm$cluster <- str_split(gpat_ppm$ID, "@", simplify = T)[, 2]
            }
            
            # calculate the auc value
            if (T) {
              file <- paste0(dir1, "/52.cs_ppm_pol2_Q", pol2_levels, "_pol2_", k, "pat_AUC_in_", data_type, "s.", cell, ".", type, ".qs")
              if (! file.exists(file)) {
                gpat_ppm$auc <- state_auc(gpat_ppm, order)
                
                qsave(gpat_ppm, file = file, nthreads = 6)
              }
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
        
        # load the gpat data
        if (T) {
          cell_dat <- function(cell) {
            file <- paste0("data/saved_data/42.gene_state_pattern/all_txs/42.", data_type, "_level_profile_Body_", up_bin_num, 
                           "UP_", down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".state.", cell, ".qs")
            dat <- qread(file)
            dat$gene_id <- paste0(dat$gene_id, "@", cell)
            
            return(dat)
          }
          dat1 <- cell_dat(cell1)
          dat2 <- cell_dat(cell2)
          
          gpat_cs <- rbind(dat1, dat2)
          
          rm(dat1, dat2)
        }
        
        k <- 4
        
        # get pol2 gene cluster: pol2cluster
        if (T) {
          cell_dat <- function(cell, type="single") {
            file <- paste0("data/saved_data/51.pol2_signal/", data_type, "/51.pol2_signal_auc_Q", pol2_levels, 
                           "_in_", data_type, "s_format.keamns", k, ".", cell, ".", type, ".qs")
            dat <- qread(file)
            
            dat$gene_id <- paste0(dat$gene_id, "@", cell)
            
            return(dat)
          }
          
          dat1 <- cell_dat(cell1, type="single")
          dat2 <- cell_dat(cell2, type="single")
          
          pol2cluster <- rbind(dat1, dat2)
          pol2cluster <- quantitle_pol2(pol2cluster, pol2_levels)
          
          rm(dat1, dat2)
        }
        
        # merge the data
        if (T) {
          over <- intersect(gpat_cs$gene_id, pol2cluster$gene_id)
          
          pol2cluster <- pol2cluster[match(over, pol2cluster$gene_id), ]
          gpat_cs <- gpat_cs[match(over, gpat_cs$gene_id), ]
          gpat_cs <- merge(gpat_cs, pol2cluster, by="gene_id")
          
          rm(over, pol2cluster)
        }
        
        # state auc: pol2 sig only
        if (T) {
          # calculate the ppm for each pol2 sig
          if (T) {
            gpat_ppm <- gene_body_batch_ppm(gpat_cs, split(gpat_cs$gene_id, gpat_cs$Quantile))
            colnames(gpat_ppm)[colnames(gpat_ppm) == "mRNA"] <- "Quantile"
          }
          
          # calculate the auc value
          if (T) {
            file <- paste0(dir1, "/52.cs_ppm_pol2_Q", pol2_levels, "_AUC_in_", data_type, "s.merged.", type, ".qs")
            if (! file.exists(file)) {
              gpat_ppm$auc <- state_auc(gpat_ppm, order)
              
              qsave(gpat_ppm, file = file, nthreads = 6)
            }
          }
        }
        
        # state auc: pol2 sig + pol2 pat
        if (T) {
          # calculate the ppm for each pol2 sig and each pol2 pat
          if (T) {
            gpat_ppm <- gene_body_batch_ppm(gpat_cs, split(gpat_cs$gene_id, paste0(gpat_cs$Quantile, "@", gpat_cs$cluster)))
            colnames(gpat_ppm)[colnames(gpat_ppm) == "mRNA"] <- "ID"
            
            gpat_ppm$Quantile <- str_split(gpat_ppm$ID, "@", simplify = T)[, 1]
            gpat_ppm$cluster <- str_split(gpat_ppm$ID, "@", simplify = T)[, 2]
          }
          
          # calculate the auc value
          if (T) {
            file <- paste0(dir1, "/52.cs_ppm_pol2_Q", pol2_levels, "_pol2_", k, "pat_AUC_in_", data_type, "s.merged.", type, ".qs")
            if (! file.exists(file)) {
              gpat_ppm$auc <- state_auc(gpat_ppm, order)
              
              qsave(gpat_ppm, file = file, nthreads = 6)
            }
          }
        }
      }
    }
  }
}

# plot the auc tendency
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir1 <- paste0("data/saved_data/52.pol2_ppm/", data_type)
      if (! dir.exists(dir1)) {
        dir.create(dir1, showWarnings = F, recursive = T)
      }
      
      dir2 <- paste0("results/2.pic/52.pol2_ppm/", data_type)
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
          
          # state auc: pol2 sig only
          if (T) {
            # load the auc value
            if (T) {
              file <- paste0(dir1, "/52.cs_ppm_pol2_Q", pol2_levels, "_AUC_in_", data_type, "s.", cell, ".", type, ".qs")
              gpat_ppm <- qread(file)
            }
            
            # common value
            if (T) {
              all_pol2_sig <- paste0("Q", 1:pol2_levels)
              all_states <- sort(unique(as.numeric(gpat_ppm$state)))
            }
            
            # format the data
            if (T) {
              head(gpat_ppm)
              
              ggdat <- melt(gpat_ppm, id.vars = c("state", "Quantile"), measure.vars = order, variable.name = "Loc", value.name = "Percentage")
              colnames(ggdat)[colnames(ggdat) == "Quantile"] <- "pol2_sig"
              
              ggdat$state <- factor(paste0("S", ggdat$state), levels = paste0("S", all_states))
              ggdat$pol2_sig <- factor(ggdat$pol2_sig, levels = all_pol2_sig)
              ggdat$Loc <- factor(ggdat$Loc, levels = order)
            }
            
            # barplot for cs at each location (motif)
            if (T) {
              # color setting
              col_state <- color_palette(length(all_states), c("#6a994e", "#fca311", "#C1121F", "#d4a373", "#3a0ca3"))
              
              head(ggdat)
              
              # ggplot (script43)
              p <- ggplot(ggdat) +
                geom_bar(aes(x=Loc, y=Percentage, fill=state), color="black", position="stack", stat="identity") +
                ggtitle(paste0(type, ": ", cell, " (cs)")) +
                scale_fill_manual(values = col_state, breaks = paste0("S", all_states)) +
                scale_x_discrete(name=NULL, breaks = labels) +
                facet_grid(pol2_sig~.) +
                cowplot::theme_cowplot() +
                theme(panel.border = element_rect(color="black"), 
                      strip.background = element_rect(fill=NA, color=NA), 
                      legend.position = "bottom") +
                guides(fill = guide_legend(title="Chromatin\nState", nrow=2, title.hjust=0.5))
              
              # plotly::ggplotly(p)
              file <- paste0(dir2, "/52.pol2_ppm_barplot_cs_in_", data_type, "s.", cell, ".", type, ".pdf")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
            
            # state percentage change across pol2 sig
            if (T) {
              head(ggdat)
              ggdat$group <- paste0(ggdat$state, "@", ggdat$pol2_sig)
              
              title <- paste0(cell, ": ", type)
              p <- state_percentage_across_pol2_sig(ggdat, title)
              
              file <- paste0(dir2, "/52.CS_percentage_across_pol2_sig_in_", data_type, "s.", cell, ".", type, ".pdf")
              ggsave(filename = file, plot = p, width = 10, height = 12)
            }
            
            # state percentage auc
            if (T) {
              colnames(gpat_ppm)[colnames(gpat_ppm) == "Quantile"] <- "pol2_sig"
              gpat_ppm$state <- factor(paste0("S", gpat_ppm$state), levels = paste0("S", all_states))
              gpat_ppm$pol2_sig <- factor(gpat_ppm$pol2_sig, levels = all_pol2_sig)
              
              title <- paste0(cell, ": ", type)
              p <- plot_state_auc_tendency(gpat_ppm, title)
              
              file <- paste0(dir2, "/52.CS_AUC_across_pol2_sig_in_", data_type, "s.", cell, ".", type, ".pdf")
              ggsave(filename = file, plot = p, width = 10, height = 12)
            }
          }
          
          k <- 4
          
          # state auc: pol2 sig + pol2 pat
          if (T) {
            # load the auc value
            if (T) {
              file <- paste0(dir1, "/52.cs_ppm_pol2_Q", pol2_levels, "_pol2_", k, "pat_AUC_in_", data_type, "s.", cell, ".", type, ".qs")
              gpat_ppm <- qread(file)
            }
            
            # common value
            if (T) {
              all_pol2_sig <- paste0("Q", 1:pol2_levels)
              all_pol2_pat <- paste0("pol2_", 1:k)
              all_states <- sort(unique(as.numeric(gpat_ppm$state)))
            }
            
            # format the data
            if (T) {
              head(gpat_ppm)
              
              ggdat <- melt(gpat_ppm, id.vars = c("state", "Quantile", "cluster"), measure.vars = order, variable.name = "Loc", value.name = "Percentage")
              colnames(ggdat)[colnames(ggdat) == "Quantile"] <- "pol2_sig"
              colnames(ggdat)[colnames(ggdat) == "cluster"] <- "pol2_pat"
              
              ggdat$state <- factor(paste0("S", ggdat$state), levels = paste0("S", all_states))
              ggdat$pol2_sig <- factor(ggdat$pol2_sig, levels = all_pol2_sig)
              ggdat$pol2_pat <- factor(paste0("pol2_", ggdat$pol2_pat), levels = all_pol2_pat)
              ggdat$Loc <- factor(ggdat$Loc, levels = order)
            }
            
            # barplot for cs at each location (motif)
            if (T) {
              # color setting
              col_state <- color_palette(length(all_states), c("#6a994e", "#fca311", "#C1121F", "#d4a373", "#3a0ca3"))
              
              head(ggdat)
              
              # ggplot (script43)
              p <- ggplot(ggdat) +
                geom_bar(aes(x=Loc, y=Percentage, fill=state), color="black", position="stack", stat="identity") +
                ggtitle(paste0(type, ": ", cell, " (cs)")) +
                scale_fill_manual(values = col_state, breaks = paste0("S", all_states)) +
                scale_x_discrete(name=NULL, breaks = labels) +
                facet_grid(pol2_sig~pol2_pat) +
                cowplot::theme_cowplot() +
                theme(panel.border = element_rect(color="black"), 
                      strip.background = element_rect(fill=NA, color=NA), 
                      legend.position = "bottom") +
                guides(fill = guide_legend(title="Chromatin\nState", nrow=2, title.hjust=0.5))
              
              # plotly::ggplotly(p)
              file <- paste0(dir2, "/52.pol2_ppm_barplot_cs_in_", data_type, "s.pol2_", k, "pat.", cell, ".", type, ".pdf")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
            
            # state percentage change across pol2 sig
            if (T) {
              head(ggdat)
              ggdat$group <- paste0(ggdat$state, "@", ggdat$pol2_sig, "@", ggdat$pol2_pat)
              
              for (pol2_pat in all_pol2_pat) {
                # pol2_pat <- all_pol2_pat[1]
                title <- paste0(cell, ": ", type, " ", pol2_pat)
                p <- state_percentage_across_pol2_sig(ggdat[ggdat$pol2_pat == pol2_pat, ], title)
                
                file <- paste0(dir2, "/52.CS_percentage_across_pol2_sig_in_", data_type, "s.", pol2_pat, ".", cell, ".", type, ".pdf")
                ggsave(filename = file, plot = p, width = 10, height = 12)
              }
            }
            
            # state percentage auc
            if (T) {
              colnames(gpat_ppm)[colnames(gpat_ppm) == "Quantile"] <- "pol2_sig"
              colnames(gpat_ppm)[colnames(gpat_ppm) == "cluster"] <- "pol2_pat"
              
              gpat_ppm$state <- factor(paste0("S", gpat_ppm$state), levels = paste0("S", all_states))
              gpat_ppm$pol2_sig <- factor(gpat_ppm$pol2_sig, levels = all_pol2_sig)
              gpat_ppm$pol2_pat <- factor(paste0("pol2_", gpat_ppm$pol2_pat), levels = all_pol2_pat)
              
              for (pol2_pat in all_pol2_pat) {
                # pol2_pat <- all_pol2_pat[1]
                title <- paste0(cell, ": ", type, " ", pol2_pat)
                p <- plot_state_auc_tendency(gpat_ppm[gpat_ppm$pol2_pat == pol2_pat, ], title)
                
                file <- paste0(dir2, "/52.CS_AUC_across_pol2_sig_in_", data_type, "s.", pol2_pat, ".", cell, ".", type, ".pdf")
                ggsave(filename = file, plot = p, width = 10, height = 12)
              }
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
        
        # state auc: pol2 sig only
        if (T) {
          # load the auc value
          if (T) {
            file <- paste0(dir1, "/52.cs_ppm_pol2_Q", pol2_levels, "_AUC_in_", data_type, "s.merged.", type, ".qs")
            gpat_ppm <- qread(file)
          }
          
          # common value
          if (T) {
            all_pol2_sig <- paste0("Q", 1:pol2_levels)
            all_states <- sort(unique(as.numeric(gpat_ppm$state)))
          }
          
          # format the data
          if (T) {
            head(gpat_ppm)
            
            ggdat <- melt(gpat_ppm, id.vars = c("state", "Quantile"), measure.vars = order, variable.name = "Loc", value.name = "Percentage")
            colnames(ggdat)[colnames(ggdat) == "Quantile"] <- "pol2_sig"
            
            ggdat$state <- factor(paste0("S", ggdat$state), levels = paste0("S", all_states))
            ggdat$pol2_sig <- factor(ggdat$pol2_sig, levels = all_pol2_sig)
            ggdat$Loc <- factor(ggdat$Loc, levels = order)
          }
          
          # barplot for cs at each location (motif)
          if (T) {
            # color setting
            col_state <- color_palette(length(all_states), c("#6a994e", "#fca311", "#C1121F", "#d4a373", "#3a0ca3"))
            
            head(ggdat)
            
            # ggplot (script43)
            p <- ggplot(ggdat) +
              geom_bar(aes(x=Loc, y=Percentage, fill=state), color="black", position="stack", stat="identity") +
              ggtitle(paste0(type, ": merged (cs)")) +
              scale_fill_manual(values = col_state, breaks = paste0("S", all_states)) +
              scale_x_discrete(name=NULL, breaks = labels) +
              facet_grid(pol2_sig~.) +
              cowplot::theme_cowplot() +
              theme(panel.border = element_rect(color="black"), 
                    strip.background = element_rect(fill=NA, color=NA), 
                    legend.position = "bottom") +
              guides(fill = guide_legend(title="Chromatin\nState", nrow=2, title.hjust=0.5))
            
            # plotly::ggplotly(p)
            file <- paste0(dir2, "/52.pol2_ppm_barplot_cs_in_", data_type, "s.merged.", type, ".pdf")
            ggsave(filename = file, plot = p, width = 15, height = 10)
          }
          
          # state percentage change across pol2 sig
          if (T) {
            head(ggdat)
            ggdat$group <- paste0(ggdat$state, "@", ggdat$pol2_sig)
            
            title <- paste0("merged: ", type)
            p <- state_percentage_across_pol2_sig(ggdat, title)
            
            file <- paste0(dir2, "/52.CS_percentage_across_pol2_sig_in_", data_type, "s.merged.", type, ".pdf")
            ggsave(filename = file, plot = p, width = 10, height = 12)
          }
          
          # state percentage auc
          if (T) {
            colnames(gpat_ppm)[colnames(gpat_ppm) == "Quantile"] <- "pol2_sig"
            gpat_ppm$state <- factor(paste0("S", gpat_ppm$state), levels = paste0("S", all_states))
            gpat_ppm$pol2_sig <- factor(gpat_ppm$pol2_sig, levels = all_pol2_sig)
            
            title <- paste0("merged: ", type)
            p <- plot_state_auc_tendency(gpat_ppm, title)
            
            file <- paste0(dir2, "/52.CS_AUC_across_pol2_sig_in_", data_type, "s.merged.", type, ".pdf")
            ggsave(filename = file, plot = p, width = 10, height = 12)
          }
        }
        
        k <- 4
        
        # state auc: pol2 sig + pol2 pat
        if (T) {
          # load the auc value
          if (T) {
            file <- paste0(dir1, "/52.cs_ppm_pol2_Q", pol2_levels, "_pol2_", k, "pat_AUC_in_", data_type, "s.merged.", type, ".qs")
            gpat_ppm <- qread(file)
          }
          
          # common value
          if (T) {
            all_pol2_sig <- paste0("Q", 1:pol2_levels)
            all_pol2_pat <- paste0("pol2_", 1:k)
            all_states <- sort(unique(as.numeric(gpat_ppm$state)))
          }
          
          # format the data
          if (T) {
            head(gpat_ppm)
            
            ggdat <- melt(gpat_ppm, id.vars = c("state", "Quantile", "cluster"), measure.vars = order, variable.name = "Loc", value.name = "Percentage")
            colnames(ggdat)[colnames(ggdat) == "Quantile"] <- "pol2_sig"
            colnames(ggdat)[colnames(ggdat) == "cluster"] <- "pol2_pat"
            
            ggdat$state <- factor(paste0("S", ggdat$state), levels = paste0("S", all_states))
            ggdat$pol2_sig <- factor(ggdat$pol2_sig, levels = all_pol2_sig)
            ggdat$pol2_pat <- factor(paste0("pol2_", ggdat$pol2_pat), levels = all_pol2_pat)
            ggdat$Loc <- factor(ggdat$Loc, levels = order)
          }
          
          # barplot for cs at each location (motif)
          if (T) {
            # color setting
            col_state <- color_palette(length(all_states), c("#6a994e", "#fca311", "#C1121F", "#d4a373", "#3a0ca3"))
            
            head(ggdat)
            
            # ggplot (script43)
            p <- ggplot(ggdat) +
              geom_bar(aes(x=Loc, y=Percentage, fill=state), color="black", position="stack", stat="identity") +
              ggtitle(paste0(type, ": merged (cs)")) +
              scale_fill_manual(values = col_state, breaks = paste0("S", all_states)) +
              scale_x_discrete(name=NULL, breaks = labels) +
              facet_grid(pol2_sig~pol2_pat) +
              cowplot::theme_cowplot() +
              theme(panel.border = element_rect(color="black"), 
                    strip.background = element_rect(fill=NA, color=NA), 
                    legend.position = "bottom") +
              guides(fill = guide_legend(title="Chromatin\nState", nrow=2, title.hjust=0.5))
            
            # plotly::ggplotly(p)
            file <- paste0(dir2, "/52.pol2_ppm_barplot_cs_in_", data_type, "s.pol2_", k, "pat.merged.", type, ".pdf")
            ggsave(filename = file, plot = p, width = 15, height = 10)
          }
          
          # state percentage change across pol2 sig
          if (T) {
            head(ggdat)
            ggdat$group <- paste0(ggdat$state, "@", ggdat$pol2_sig, "@", ggdat$pol2_pat)
            
            for (pol2_pat in all_pol2_pat) {
              # pol2_pat <- all_pol2_pat[1]
              title <- paste0("merged: ", type, " ", pol2_pat)
              p <- state_percentage_across_pol2_sig(ggdat[ggdat$pol2_pat == pol2_pat, ], title)
              
              file <- paste0(dir2, "/52.CS_percentage_across_pol2_sig_in_", data_type, "s.", pol2_pat, ".merged.", type, ".pdf")
              ggsave(filename = file, plot = p, width = 10, height = 12)
            }
          }
          
          # state percentage auc
          if (T) {
            colnames(gpat_ppm)[colnames(gpat_ppm) == "Quantile"] <- "pol2_sig"
            colnames(gpat_ppm)[colnames(gpat_ppm) == "cluster"] <- "pol2_pat"
            
            gpat_ppm$state <- factor(paste0("S", gpat_ppm$state), levels = paste0("S", all_states))
            gpat_ppm$pol2_sig <- factor(gpat_ppm$pol2_sig, levels = all_pol2_sig)
            gpat_ppm$pol2_pat <- factor(paste0("pol2_", gpat_ppm$pol2_pat), levels = all_pol2_pat)
            
            for (pol2_pat in all_pol2_pat) {
              # pol2_pat <- all_pol2_pat[1]
              title <- paste0("merged: ", type, " ", pol2_pat)
              p <- plot_state_auc_tendency(gpat_ppm[gpat_ppm$pol2_pat == pol2_pat, ], title)
              
              file <- paste0(dir2, "/52.CS_AUC_across_pol2_sig_in_", data_type, "s.", pol2_pat, ".merged.", type, ".pdf")
              ggsave(filename = file, plot = p, width = 10, height = 12)
            }
          }
        }
      }
    }
  }
}


