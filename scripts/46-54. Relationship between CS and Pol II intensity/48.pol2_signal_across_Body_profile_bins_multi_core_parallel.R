library(stringr)
library(reshape2)
library(qs)
library(doParallel)
library(parallel)
library(foreach)
cl <- makeCluster(20, outfile="log48.txt")
registerDoParallel(cl)

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
up_bin_num <- 5
down_bin_num <- 5
body_bin_num <- 10
rna_levels <- 10

# define the funtion
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
  
  # mean matrix
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
  plot_body_distribution <- function(dat, title, label_y, group=0, group_col=NA, hide_legend=T) {
    # test data
    if (F) {
      title <- "test"
      label_y <- "pol2 sig"
      group <- 0
      group_col <- NA
      hide_legend <- T
    }
    
    # ggplot
    if (group == 0) {
      p <- ggplot(dat) +
        geom_segment(aes(x=Loc, y=min(Value)*0.8, xend=Loc, yend=Assistant), linetype=2, linewidth=0.8, color="grey", na.rm=T, alpha=0.6) +
        geom_line(aes(x=Loc, y=Value), group=1, linewidth=1, alpha=1) +
        scale_x_discrete(name = NULL, breaks = label) +
        scale_y_continuous(name = label_y, limits = c(min(dat$Value)*0.8, NA)) +
        ggtitle(title) +
        cowplot::theme_cowplot() +
        theme(axis.text.x = element_text(size = rel(0.9)))
    }
    
    # facet by group
    if (group > 0) {
      colnames(dat)[group] <- "group"
      
      p <- ggplot(dat) +
        geom_segment(aes(x=Loc, y=min(Value)*0.8, xend=Loc, yend=Assistant, group=group), linetype=2, linewidth=0.8, color="grey", na.rm=T, alpha=0.6) +
        geom_line(aes(x=Loc, y=Value, group=group, color=group), linewidth=1, alpha=1) +
        scale_x_discrete(name = NULL, breaks = label) +
        scale_y_continuous(name = label_y, limits = c(min(dat$Value)*0.8, NA)) +
        ggtitle(title) +
        cowplot::theme_cowplot() +
        theme(axis.text.x = element_text(size = rel(0.9)))
      
      # change the color of group
      if (sum(is.na(group_col)) == 0 & length(group_col) >= 1) {
        p <- p +
          scale_color_manual(values = group_col)
      }
    }
    
    if (hide_legend) {
      p <- p + theme(legend.position = "none")
    }
    
    return(p)
  }
}

# data prepare: cell specific tss/body matrix for each gene
single_dat <- function(mk) {
  # load the environment
  if (T) {
    library(stringr)
    library(reshape2)
    library(qs)
  }
  
  # common value
  if (T) {
    cell1 <- "thp1"
    cell2 <- "cd34"
    bin_size <- 200
    up_bin_num <- 5
    down_bin_num <- 5
    body_bin_num <- 10
    rna_levels <- 10
  }
  
  # define the funtion
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
    
    # mean matrix
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
    plot_body_distribution <- function(dat, title, label_y, group=0, group_col=NA, hide_legend=T) {
      # test data
      if (F) {
        title <- "test"
        label_y <- "pol2 sig"
        group <- 0
        group_col <- NA
        hide_legend <- T
      }
      
      # ggplot
      if (group == 0) {
        p <- ggplot(dat) +
          geom_segment(aes(x=Loc, y=min(Value)*0.8, xend=Loc, yend=Assistant), linetype=2, linewidth=0.8, color="grey", na.rm=T, alpha=0.6) +
          geom_line(aes(x=Loc, y=Value), group=1, linewidth=1, alpha=1) +
          scale_x_discrete(name = NULL, breaks = label) +
          scale_y_continuous(name = label_y, limits = c(min(dat$Value)*0.8, NA)) +
          ggtitle(title) +
          cowplot::theme_cowplot() +
          theme(axis.text.x = element_text(size = rel(0.9)))
      }
      
      # facet by group
      if (group > 0) {
        colnames(dat)[group] <- "group"
        
        p <- ggplot(dat) +
          geom_segment(aes(x=Loc, y=min(Value)*0.8, xend=Loc, yend=Assistant, group=group), linetype=2, linewidth=0.8, color="grey", na.rm=T, alpha=0.6) +
          geom_line(aes(x=Loc, y=Value, group=group, color=group), linewidth=1, alpha=1) +
          scale_x_discrete(name = NULL, breaks = label) +
          scale_y_continuous(name = label_y, limits = c(min(dat$Value)*0.8, NA)) +
          ggtitle(title) +
          cowplot::theme_cowplot() +
          theme(axis.text.x = element_text(size = rel(0.9)))
        
        # change the color of group
        if (sum(is.na(group_col)) == 0 & length(group_col) >= 1) {
          p <- p +
            scale_color_manual(values = group_col)
        }
      }
      
      if (hide_legend) {
        p <- p + theme(legend.position = "none")
      }
      
      return(p)
    }
  }
  
  # loop body
  if (T) {
    # test data
    if (F) {
      data_type <- "tx"
      type <- "single"
      cell <- cell1
      mk <- "H3K27ac"
    }
    
    # read gene body state matrix
    if (T) {
      # mat_file <- paste0("data/saved_data/41.", data_type, "s_with_nonoverlap_region_matrix_bs", bin_size, ".qs")
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
      file <- paste0(dir1, "/48.", data_type, "_Body_", up_bin_num, "UP_", down_bin_num, "DW_", 
                     body_bin_num, "Body.bs", bin_size, ".", mk, ".mean_mat.", cell, ".", type, ".qs")
      
      if (file.exists(file)) {
        mat <- qread(file, nthreads = 6)
      }
      if (! file.exists(file)) {
        mat <- get_mk_sig_across_body_matrix(dat, gene_body_mat)
        
        qsave(mat, file, nthreads = 6)
      }
    }
  }
}
merged_dat <- function(mk) {
  # load the environment
  if (T) {
    library(stringr)
    library(reshape2)
    library(qs)
  }
  
  # common value
  if (T) {
    cell1 <- "thp1"
    cell2 <- "cd34"
    bin_size <- 200
    up_bin_num <- 5
    down_bin_num <- 5
    body_bin_num <- 10
    rna_levels <- 10
  }
  
  # define the funtion
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
    
    # mean matrix
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
    plot_body_distribution <- function(dat, title, label_y, group=0, group_col=NA, hide_legend=T) {
      # test data
      if (F) {
        title <- "test"
        label_y <- "pol2 sig"
        group <- 0
        group_col <- NA
        hide_legend <- T
      }
      
      # ggplot
      if (group == 0) {
        p <- ggplot(dat) +
          geom_segment(aes(x=Loc, y=min(Value)*0.8, xend=Loc, yend=Assistant), linetype=2, linewidth=0.8, color="grey", na.rm=T, alpha=0.6) +
          geom_line(aes(x=Loc, y=Value), group=1, linewidth=1, alpha=1) +
          scale_x_discrete(name = NULL, breaks = label) +
          scale_y_continuous(name = label_y, limits = c(min(dat$Value)*0.8, NA)) +
          ggtitle(title) +
          cowplot::theme_cowplot() +
          theme(axis.text.x = element_text(size = rel(0.9)))
      }
      
      # facet by group
      if (group > 0) {
        colnames(dat)[group] <- "group"
        
        p <- ggplot(dat) +
          geom_segment(aes(x=Loc, y=min(Value)*0.8, xend=Loc, yend=Assistant, group=group), linetype=2, linewidth=0.8, color="grey", na.rm=T, alpha=0.6) +
          geom_line(aes(x=Loc, y=Value, group=group, color=group), linewidth=1, alpha=1) +
          scale_x_discrete(name = NULL, breaks = label) +
          scale_y_continuous(name = label_y, limits = c(min(dat$Value)*0.8, NA)) +
          ggtitle(title) +
          cowplot::theme_cowplot() +
          theme(axis.text.x = element_text(size = rel(0.9)))
        
        # change the color of group
        if (sum(is.na(group_col)) == 0 & length(group_col) >= 1) {
          p <- p +
            scale_color_manual(values = group_col)
        }
      }
      
      if (hide_legend) {
        p <- p + theme(legend.position = "none")
      }
      
      return(p)
    }
  }
  
  # loop body
  if (T) {
    # test data
    if (F) {
      data_type <- "tx"
      type <- "merged"
      mk <- "H3K27ac"
    }
    
    # read gene body state matrix
    if (T) {
      # mat_file <- paste0("data/saved_data/41.", data_type, "s_with_nonoverlap_region_matrix_bs", bin_size, ".qs")
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
      dat1 <- read_mk(mk, cell1)
      dat2 <- read_mk(mk, cell2)
      dat <- data.frame(
        ID = dat1$ID, 
        signal = sapply(1:nrow(dat1), function(i) {
          mean(c(dat1$signal[i], dat2$signal[i]))
        })
      )
    }
    
    # get mk specific signal (mean value in body region)
    if (T) {
      file <- paste0(dir1, "/48.", data_type, "_Body_", up_bin_num, "UP_", down_bin_num, "DW_", 
                     body_bin_num, "Body.bs", bin_size, ".", mk, ".mean_mat.merged.", type, ".qs")
      
      if (file.exists(file)) {
        mat <- qread(file, nthreads = 6)
      }
      if (! file.exists(file)) {
        mat <- get_mk_sig_across_body_matrix(dat, gene_body_mat)
        
        qsave(mat, file, nthreads = 6)
      }
    }
  }
}
if (T) {
  mks <- c("pol2", "ATAC", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K79me2", "H3K9me3")
  
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir1 <- paste0("data/saved_data/48.mk_signal_distribution/", data_type)
      if (! dir.exists(dir1)) {
        dir.create(dir1, showWarnings = F, recursive = T)
      }
      
      dir2 <- paste0("results/2.pic/48.mk_signal_distribution/", data_type)
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
          
          foreach(mk = mks, .export =ls()) %dopar% single_dat(mk)
        }
      }
      
      # merged data
      if (type == "merged") {
        cat(paste0("\t\tmerged: \n"))
        
        foreach(mk = mks, .export =ls()) %dopar% merged_dat(mk)
      }
    }
  }
}

stopCluster(cl)