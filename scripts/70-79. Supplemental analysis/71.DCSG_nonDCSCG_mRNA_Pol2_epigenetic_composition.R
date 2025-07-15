# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(stringr)
  library(reshape2)
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(circlize))
  suppressPackageStartupMessages(library(clusterProfiler))
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  library(qs)
}

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
up_bin_num <- 5
down_bin_num <- 5
body_bin_num <- 10

# define function
if (T) {
  # script 48
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
    
    # data prepare
    get_mk_sig_across_body_matrix <- function(dat, gene_body_mat) {
      # prepare hello info
      if (T) {
        start_mess <- paste0("|", paste(rep("-", 100), collapse = ""), "|\n")
        cat(start_mess)
        
        gene_body_mat$num <- 1:nrow(gene_body_mat)
        breaks <- round(seq(1, nrow(gene_body_mat), length.out=100))
      }
      
      # get cell specific pol2 matrix
      if (T) {
        mat <- data.frame(t(
          apply(gene_body_mat, 1, function(x) {
            # x <- unlist(gene_body_mat[1, ])
            # print process
            if (as.numeric(x[length(x)]) %in% breaks) {
              if (which(as.numeric(x[length(x)]) == breaks) == 1) {
                cat("|*")
              }
              if (which(as.numeric(x[length(x)]) == breaks) == 100) {
                cat("*|\n")
              }
              if (! which(as.numeric(x[length(x)]) == breaks) %in% c(1, 100)) {
                cat("*")
              }
            }
            
            # nonbody
            if (T) {
              nonbody <- x[grepl("^U|TSS|TES|^D", colnames(gene_body_mat))]
              nonbody <- as.numeric(nonbody)
              nonbody <- dat[nonbody, "signal"]
              names(nonbody) <- grep("^U|TSS|TES|^D", colnames(gene_body_mat), value = T)
            }
            
            # genebody
            if (T) {
              genebody_bin <- x[grepl("^B", colnames(gene_body_mat))]
              genebody <- sapply(genebody_bin, function(gb) {
                # gb <- genebody_bin[1]
                start <- as.numeric(strsplit(gb, "-")[[1]][1])
                end <- as.numeric(strsplit(gb, "-")[[1]][2])
                pol_dat <- dat[start:end, "signal"]
                mean <- mean(pol_dat)
                
                return(mean)
              })
            }
            
            gene_dat <- c(nonbody, genebody)
            gene_dat <- gene_dat[name]
            
            # return result
            return(gene_dat)
          })
        ))
        
        mat$gene_id <- gene_body_mat$gene_id
      }
      
      return(mat)
    }
  }
  
  # script 51
  if (T) {
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
  }
  
  mkdir_fun <- function(dir) {
    for (d in dir) {
      if (! dir.exists(d)) {
        dir.create(d, showWarnings = F, recursive = T)
      }
    }
  }
}

# data prepare: cell specific body matrix for each gene
if (T) {
  mks <- c("pol2", "ATAC", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K79me2", "H3K9me3")
  
  for (data_type in c("gene")) {
    cat(paste0(data_type, ": \n"))
    
    for (type in c("single")) {
      cat(paste0("\t", type, ": \n"))
      
      # cell specific data
      for (cell in c(cell1, cell2)) {
        cat(paste0("\t\t", cell, ": \n"))
        
        for (mk in mks) {
          cat(paste0("\t\t\t", mk, ": \n"))
          
          # test data
          if (F) {
            data_type <- "gene"
            type <- "single"
            cell <- cell1
            mk <- "H3K27ac"
          }
          
          # mkdir
          if (T) {
            dir1 <- paste0("data/saved_data/71.mk_signal_distribution/", data_type)
            if (! dir.exists(dir1)) {
              dir.create(dir1, showWarnings = F, recursive = T)
            }
            
            dir2 <- paste0("results/2.pic/71.mk_signal_distribution/", data_type)
            if (! dir.exists(dir2)) {
              dir.create(dir2, showWarnings = F, recursive = T)
            }
          }
          
          # read gene body state matrix
          if (T) {
            mat_file <- paste0("data/saved_data/11.", data_type, "_level_profile_Body_", up_bin_num, "UP_", 
                               down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".bin.qs")
            gene_body_mat <- qread(mat_file, nthreads = 6)
            
            rm(mat_file)
          }
          
          # read cell specific mk signal
          if (T) {
            read_mk <- function(mk, cell) {
              path <- ifelse(mk == "pol2", 
                             paste0("data/raw_data/6.pol2/"), 
                             paste0("data/raw_data/7.histone/"))
              file <- paste0(path, "bs", bin_size, "/", cell, ".", mk, ".S3V2.bedgraph.NBP.txt")
              file <- ifelse(mk == "pol2", file, paste0(file, ".gz"))
              
              dat <- data.table::fread(file, header = F, sep = "\t", data.table = F)
              dat <- data.frame(ID = 1:nrow(dat), 
                                signal = dat[, 1])
              
              return(dat)
            }
            
            dat <- read_mk(mk, cell)
          }
          
          # get mk specific signal (mean value in body region)
          if (T) {
            file <- paste0(dir1, "/71.", data_type, "_Body_", up_bin_num, "UP_", down_bin_num, "DW_", 
                           body_bin_num, "Body.bs", bin_size, ".", mk, ".mean_mat.", cell, ".", type, ".qs")
            
            if (file.exists(file)) {
              mat <- qread(file, nthreads = 6)
            }
            if (! file.exists(file)) {
              mat <- get_mk_sig_across_body_matrix(dat, gene_body_mat)
              
              qsave(mat, file, nthreads = 6)
            }
          }
          
          # calculate the AUC value for each ID
          if (T) {
            file <- paste0(dir1, "/71.signal_auc_of_", mk, "_in_", data_type, "s_profile_Body.", cell, ".", type, ".qs")
            
            if (! file.exists(file)) {
              signal <- gene_pol2_signal(mat)
              
              qsave(signal, file, nthreads = 6)
            }
          }
        }
      }
    }
  }
}

# statistics of mRNA level, Pol2 signal, epigenetic signal composition for DCSG but non-DCSCG
if (T) {
  for (data_type in c("gene")) {
    cat(paste0(data_type, ": \n"))
    
    for (method in c("chromIDEAS")) {
      cat(paste0("\t", method, ": \n"))
      
      # test data
      if (F) {
        data_type <- "gene"
        method <- "chromIDEAS"
      }
      
      # mkdir
      dir <- paste0("results/2.pic/71.DCSG_nonDCSCG_comparison/", data_type)
      mkdir_fun(dir)
      
      # load DCSG but non-DCSCG
      if (T) {
        file <- paste0("results/1.tab/58.DCSG_conversion/", data_type, "/58.DCSCG_geneID.", method, ".txt")
        dcscg <- read.table(file, header = F)[, 1]
        
        file <- paste0("results/1.tab/58.DCSG_conversion/", data_type, "/58.DCSG_nonDCSCG_geneID.", method, ".txt")
        dcsg_nondcscg <- read.table(file, header = F)[, 1]
      }
      
      # mRNA dat
      if (T) {
        # load target gene expression
        if (T) {
          get_rna <- function(file, id) {
            dat <- qread(file, nthreads = 6)$sorted_decreasing_tx
            dat <- dat[id]
            
            return(dat)
          }
          
          thp1_dcscg <- get_rna("data/saved_data/6.4_type_gene_expression_level_classification_thp1.qs", dcscg)
          cd34_dcscg <- get_rna("data/saved_data/6.4_type_gene_expression_level_classification_cd34.qs", dcscg)
          
          thp1_dcsg_nondcscg <- get_rna("data/saved_data/6.4_type_gene_expression_level_classification_thp1.qs", dcsg_nondcscg)
          cd34_dcsg_nondcscg <- get_rna("data/saved_data/6.4_type_gene_expression_level_classification_cd34.qs", dcsg_nondcscg)
          
          mrna <- data.frame(id = c(dcscg, dcsg_nondcscg), 
                             type = c(rep("dcscg", length(dcscg)), rep("dcsg_nondcscg", length(dcsg_nondcscg))), 
                             thp1 = c(thp1_dcscg, thp1_dcsg_nondcscg), 
                             cd34 = c(cd34_dcscg, cd34_dcsg_nondcscg))
          
          non0min <- c(mrna$thp1, mrna$cd34)
          non0min <- min(non0min[non0min>0])
          
          mrna$log2fc <- log2((mrna$thp1+non0min) / (mrna$cd34+non0min))
          mrna$log2fc_abs <- abs(mrna$log2fc)
        }
        
        # statistics
        if (T) {
          # dcscg
          if (T) {
            wilcox.test(thp1_dcscg, cd34_dcscg, paired = T)
            t.test(thp1_dcscg, cd34_dcscg, paired = T)
          }
          
          # dcsg_nondcscg
          if (T) {
            wilcox.test(thp1_dcsg_nondcscg, cd34_dcsg_nondcscg, paired = T)
            t.test(thp1_dcsg_nondcscg, cd34_dcsg_nondcscg, paired = T)
          }
        }
        
        # format the data
        if (T) {
          head(mrna)
          mrna_format <- melt(mrna, id.vars = c("id", "type"), measure.vars = c(cell1, cell2), variable.name = "cell", value.name = "signal")
          
          rm(thp1_dcscg, thp1_dcsg_nondcscg, cd34_dcscg, cd34_dcsg_nondcscg, non0min, get_rna)
        }
      }
      
      # epi and pol2 dat
      if (T) {
        mks <- c("pol2", "ATAC", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K79me2", "H3K9me3")
        
        # load target gene epi signal
        if (T) {
          get_mk <- function(mk, cell, id) {
            # test data
            if (F) {
              mk <- "pol2"
              cell <- cell1
            }
            
            file <- paste0("data/saved_data/71.mk_signal_distribution/", data_type, "/71.signal_auc_of_", mk, "_in_", data_type, "s_profile_Body.", cell, ".single.qs")
            
            dat <- qread(file, nthreads = 6)
            signal <- dat$signal
            names(signal) <- dat$id
            
            dat <- signal[id]
            
            return(dat)
          }
          mk_stat <- function(mk, dcscg, dcsg_nondcscg) {
            # mk <- "pol2"
            
            # load raw data
            if (T) {
              thp1_dcscg <- get_mk(mk, "thp1", dcscg)
              thp1_dcscg <- thp1_dcscg[! is.na(thp1_dcscg)]
              cd34_dcscg <- get_mk(mk, "cd34", dcscg)
              cd34_dcscg <- cd34_dcscg[! is.na(cd34_dcscg)]
              
              thp1_dcsg_nondcscg <- get_mk(mk, "thp1", dcsg_nondcscg)
              cd34_dcsg_nondcscg <- get_mk(mk, "cd34", dcsg_nondcscg)
              
              dat <- data.frame(id = names(c(thp1_dcscg, thp1_dcsg_nondcscg)), 
                                type = c(rep("dcscg", length(thp1_dcscg)), rep("dcsg_nondcscg", length(thp1_dcsg_nondcscg))), 
                                thp1 = c(thp1_dcscg, thp1_dcsg_nondcscg), 
                                cd34 = c(cd34_dcscg, cd34_dcsg_nondcscg))
              
              dat$log2fc <- log2(dat$thp1/dat$cd34)
              dat$log2fc_abs <- abs(dat$log2fc)
            }
            
            # statistics
            if (T) {
              # dcscg
              if (T) {
                dcscg_w <- wilcox.test(thp1_dcscg, cd34_dcscg, paired = T)$p.value
                dcscg_t <- t.test(thp1_dcscg, cd34_dcscg, paired = T)$p.value
              }
              
              # dcsg_nondcscg
              if (T) {
                dcsg_nondcscg_w <- wilcox.test(thp1_dcsg_nondcscg, cd34_dcsg_nondcscg, paired = T)$p.value
                dcsg_nondcscg_t <- t.test(thp1_dcsg_nondcscg, cd34_dcsg_nondcscg, paired = T)$p.value
              }
              
              # log2
              if (T) {
                log2_abs_w <- wilcox.test(dat$log2fc_abs[dat$type == "dcscg"], dat$log2fc_abs[dat$type == "dcsg_nondcscg"])$p.value
                log2_abs_t <- t.test(dat$log2fc_abs[dat$type == "dcscg"], dat$log2fc_abs[dat$type == "dcsg_nondcscg"])$p.value
              }
            }
            
            # mess
            if (T) {
              mess <- paste0("\t", mk, ":\ndcscg\t", round(dcscg_w, 4), "\t", round(dcscg_t, 4), 
                             "\ndcsg_nondcscg\t", round(dcsg_nondcscg_w, 4), "\t", round(dcsg_nondcscg_t, 4), 
                             "\nlog2fc_abs\t", round(log2_abs_w, 4), "\t", round(log2_abs_t, 4), "\n")
              cat(mess)
            }
            
            return(dat)
          }
          
          dat <- lapply(mks, function(mk) {
            res <- mk_stat(mk, dcscg, dcsg_nondcscg)
            res$mk <- mk
            
            return(res)
          })
          dat <- do.call(rbind, dat)
          
          rm(get_mk, mk_stat)
        }
        
        # format the data
        if (T) {
          head(dat)
          dat_format <- melt(dat, id.vars = c("id", "type", "mk"), measure.vars = c(cell1, cell2), variable.name = "cell", value.name = "signal")
        }
      }
      
      # merge the data
      if (T) {
        head(dat)
        head(mrna)
        mrna$mk <- "mRNA"
        dat <- rbind(dat, mrna)
        rm(mrna)
        
        head(dat_format)
        head(mrna_format)
        mrna_format$mk <- "mRNA"
        mrna_format <- mrna_format[, colnames(dat_format)]
        dat_format <- rbind(dat_format, mrna_format)
        rm(mrna_format)
      }
      
      # format the data
      if (T) {
        dat$mk <- factor(dat$mk, levels = c("mRNA", mks))
        dat$type <- factor(dat$type, levels = c("dcscg", "dcsg_nondcscg"))
        
        dat_format$mk <- factor(dat_format$mk, levels = c("mRNA", mks))
        dat_format$cell <- factor(dat_format$cell, levels = c(cell1, cell2))
        dat_format$type <- factor(dat_format$type, levels = c("dcscg", "dcsg_nondcscg"))
      }
      
      # ggplot2
      if (T) {
        # signal distribution
        if (T) {
          head(dat_format)
          p <- ggplot(dat_format, aes(x=cell, y=log2(signal+1), group=cell, fill=cell)) +
            geom_violin(scale = "width") +
            geom_boxplot(width=0.1, fill=NA) +
            ggpubr::stat_compare_means(method="wilcox.test", paired = T, label="p.signif", label.y.npc=0.8) +
            xlab(NULL) +
            ylab("log2(Signal+1)") +
            facet_grid(mk~type, scales="free") +
            cowplot::theme_cowplot() +
            theme(axis.title = element_text(size = rel(1.2)),
                  axis.text = element_text(size = rel(1.2)),
                  panel.border = element_rect(color="black"), 
                  strip.background = element_rect(fill=NA, color=NA), 
                  legend.position = "none")
          
          file <- paste0(dir, "/71.", data_type, "_level_DCSG_nonDCSCG_mk_absolute_signal_between_cells_", method, ".pdf")
          ggsave(filename = file, plot = p, width = 8, height = 15)
        }
        
        # log2fc abs
        if (T) {
          head(dat)
          p <- ggplot(dat, aes(x=type, y=log2fc_abs, group=type, fill=type)) +
            geom_violin(scale = "width") +
            geom_boxplot(width=0.1, fill=NA) +
            ggpubr::stat_compare_means(method="wilcox.test", paired = F, label="p.signif", label.y.npc=0.8) +
            xlab(NULL) +
            ylab("Absolute log2FC") +
            facet_grid(mk~., scales="free") +
            cowplot::theme_cowplot() +
            theme(axis.title = element_text(size = rel(1.2)),
                  axis.text = element_text(size = rel(1.2)),
                  panel.border = element_rect(color="black"), 
                  strip.background = element_rect(fill=NA, color=NA), 
                  legend.position = "none")
          
          file <- paste0(dir, "/71.", data_type, "_level_DCSG_nonDCSCG_mk_absolute_log2fc_between_gene_types_", method, ".pdf")
          ggsave(filename = file, plot = p, width = 8, height = 15)
        }
      }
    }
  }
}
