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
    label <- c(paste0("U", up_bin_num), "TSS", "TES", paste0("D", down_bin_num))
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
  
  # mean matrix (script48)
  get_cell_specific_body_plot_dat <- function(mat, ids) {
    # test data
    if (F) {
      ids <- IDs[[1]]
    }
    
    # get specific IDs data
    if (T) {
      mess <- sum(! ids %in% mat$gene_id)
      if (mess > 0) {
        mess <- paste0("There are ", mess, "/", length(ids), " (", round(mess/length(ids)*100, 2), "%) genes in IDs not exist in matrix\n")
        cat(mess)
      }
      rm(mess)
      
      # get existing IDs
      ids <- ids[ids %in% mat$gene_id]
      
      # make sure the order are identical
      mat <- mat[match(ids, mat$gene_id), ]
    }
    
    # get average pol2 signal
    if (T) {
      colnames(mat)
      dat <- sapply(mat[, grep("^U|TSS|^B|TES|^D", colnames(mat))], mean)
    }
    
    # format the data
    if (T) {
      dat <- data.frame(Loc = names(dat), 
                        Value = dat)
      
      dat$Loc <- factor(dat$Loc, levels = name)
      dat$Assistant <- ifelse(dat$Loc %in% label, dat$Value, NA)
    }
    
    return(dat)
  }
  
  move_file <- function(file, dir) {
    stat <- file.copy(from = file, to = paste0(dir, "/"), overwrite = T, copy.date = T)
    if (stat) {
      stat <- file.remove(file)
    }
  }
}

# data prepare
if (T) {
  mks <- c("pol2", "ATAC", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K79me2", "H3K9me3")
  
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir <- paste0("results/2.pic/49.mk_signal_distribution_under_mRNA_levels/", data_type)
      if (! dir.exists(dir)) {
        dir.create(dir, showWarnings = F, recursive = T)
      }
    }
    
    for (type in c("merged", "single")) {
      cat(paste0("\t", type, ": \n"))
      
      # cell specific data
      if (type == "single") {
        for (cell in c(cell1, cell2)) {
          cat(paste0("\t\t", cell, ": \n"))
          
          for (mk in mks) {
            cat(paste0("\t\t\t", mk, ": \n"))
            
            # test data
            if (F) {
              data_type <- "tx"
              type <- "single"
              cell <- cell1
              mk <- "H3K27ac"
            }
            
            # get pol2 matrix
            if (T) {
              read_dat <- function(cell, mk, type) {
                input <- paste0("data/saved_data/48.mk_signal_distribution/", data_type, "/48.", data_type, "_Body_", up_bin_num, "UP_", 
                                down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".", mk, ".mean_mat.", cell, ".", type, ".qs")
                mat <- qread(input, nthreads = 6)
                
                return(mat)
              }
              mat <- read_dat(cell, mk, type)
            }
            
            # get mRNA level data: IDs
            if (T) {
              file <- ifelse(data_type == "gene", 
                             paste0("data/saved_data/6.4_type_gene_expression_level_classification_", cell, ".qs"), 
                             paste0("data/saved_data/7.4_type_tx_expression_level_classification_", cell, ".qs"))
              
              cell_mrna <- qread(file, nthreads = 6)
              cell_mrna <- cell_mrna[[1]]
              
              cell_mrna <- cell_mrna[names(cell_mrna) %in% mat$gene_id]
              IDs <- quantitle_mRNA(cell_mrna)
              
              rm(cell_mrna)
            }
            
            # get mRNA specific profile
            if (T) {
              dat <- lapply(IDs, function(x) {
                get_cell_specific_body_plot_dat(mat, x)
              })
            }
            
            # format the profile
            if (T) {
              profile <- lapply(names(dat), function(x) {
                # x <- names(dat)[1]
                ggdat <- dat[[x]]
                ggdat$mRNA <- x
                
                return(ggdat)
              })
              
              profile <- data.frame(do.call(rbind, profile))
              profile$mRNA <- factor(profile$mRNA, levels = c("subdat0", paste0("Q", 1:rna_levels)))
            }
            
            # ggplot
            if (T) {
              head(profile)
              col_fun <- circlize::colorRamp2(c(1, (rna_levels+1)), c("blue", "red"))
              
              p <- ggplot(profile) +
                geom_segment(aes(x=Loc, y=min(Value)*0.8, xend=Loc, yend=Assistant), linetype=2, linewidth=0.8, color="grey", na.rm=T, alpha=0.6) +
                geom_line(aes(x=Loc, y=Value, group=mRNA, color=mRNA), linewidth=1, alpha=1) +
                scale_x_discrete(name = NULL, breaks = label) +
                scale_y_continuous(name = paste0(mk, " Signal"), limits = c(min(profile$Value)*0.8, NA)) +
                ggtitle(cell) +
                scale_colour_manual(values = col_fun(1:(rna_levels+1)), breaks = levels(profile$mRNA)) +
                cowplot::theme_cowplot() +
                theme(panel.border = element_rect(color="black"), 
                      strip.background = element_rect(fill=NA, color=NA))
              
              file <- paste0("49.", mk, "_sig_across_mRNA_", data_type, "s_Body.", cell, ".", type, ".pdf")
              ggsave(filename = file, plot = p, width = 8, height = 6)
              move_file(file, dir)
              
              # plotly::ggplotly(p)
            }
          }
        }
      }
      
      # merged data
      if (type == "merged") {
        cat(paste0("\t\tmerged: \n"))
        
        for (mk in mks) {
          cat(paste0("\t\t\t", mk, ": \n"))
          
          # test data
          if (F) {
            data_type <- "tx"
            type <- "merged"
            mk <- "H3K27ac"
          }
          
          # get pol2 matrix
          if (T) {
            read_dat <- function(cell, mk, type) {
              input <- paste0("data/saved_data/48.mk_signal_distribution/", data_type, "/48.", data_type, "_Body_", up_bin_num, "UP_", 
                              down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".", mk, ".mean_mat.", cell, ".", type, ".qs")
              mat <- qread(input, nthreads = 6)
              
              return(mat)
            }
            
            dat1 <- read_dat(cell1, mk, type = "single")
            dat2 <- read_dat(cell2, mk, type = "single")
            ids <- dat1$gene_id
            dat1$gene_id <- paste0(dat1$gene_id, "@", cell1)
            dat2$gene_id <- paste0(dat2$gene_id, "@", cell2)
            
            mat <- rbind(dat1, dat2)
          }
          
          # get mRNA level data: IDs (script43)
          if (T) {
            rna_Q <- function(cell, ids) {
              # cell <- cell1
              file <- ifelse(data_type == "gene", 
                             paste0("data/saved_data/6.4_type_gene_expression_level_classification_", cell, ".qs"), 
                             paste0("data/saved_data/7.4_type_tx_expression_level_classification_", cell, ".qs"))
              
              cell_mrna <- qread(file, nthreads = 6)
              cell_mrna <- cell_mrna[[1]]
              
              cell_mrna <- cell_mrna[names(cell_mrna) %in% ids]
              names(cell_mrna) <- paste0(names(cell_mrna), "@", cell)
              IDs <- quantitle_mRNA(cell_mrna)
              
              return(IDs)
            }
            
            dat1 <- rna_Q(cell1, ids)
            dat2 <- rna_Q(cell2, ids)
            
            IDs <- lapply(names(dat1), function(x) {
              # x <- names(dat1)[1]
              
              subdat1 <- dat1[[x]]
              subdat2 <- dat2[[x]]
              
              return(c(subdat1, subdat2))
            })
            names(IDs) <- names(dat1)
            
            rm(dat1, dat2, ids)
          }
          
          # get mRNA specific profile
          if (T) {
            dat <- lapply(IDs, function(x) {
              get_cell_specific_body_plot_dat(mat, x)
            })
          }
          
          # format the profile
          if (T) {
            profile <- lapply(names(dat), function(x) {
              # x <- names(dat)[1]
              ggdat <- dat[[x]]
              ggdat$mRNA <- x
              
              return(ggdat)
            })
            
            profile <- data.frame(do.call(rbind, profile))
            profile$mRNA <- factor(profile$mRNA, levels = c("subdat0", paste0("Q", 1:rna_levels)))
          }
          
          # ggplot
          if (T) {
            head(profile)
            col_fun <- circlize::colorRamp2(c(1, (rna_levels+1)), c("blue", "red"))
            
            p <- ggplot(profile) +
              geom_segment(aes(x=Loc, y=min(Value)*0.8, xend=Loc, yend=Assistant), linetype=2, linewidth=0.8, color="grey", na.rm=T, alpha=0.6) +
              geom_line(aes(x=Loc, y=Value, group=mRNA, color=mRNA), linewidth=1, alpha=1) +
              scale_x_discrete(name = NULL, breaks = label) +
              scale_y_continuous(name = paste0(mk, " Signal"), limits = c(min(profile$Value)*0.8, NA)) +
              ggtitle("merged") +
              scale_colour_manual(values = col_fun(1:(rna_levels+1)), breaks = levels(profile$mRNA)) +
              cowplot::theme_cowplot() +
              theme(panel.border = element_rect(color="black"), 
                    strip.background = element_rect(fill=NA, color=NA))
            
            file <- paste0("49.", mk, "_sig_across_mRNA_", data_type, "s_Body.merged.", type, ".pdf")
            ggsave(filename = file, plot = p, width = 8, height = 6)
            move_file(file, dir)
            
            # plotly::ggplotly(p)
          }
        }
      }
    }
  }
}

