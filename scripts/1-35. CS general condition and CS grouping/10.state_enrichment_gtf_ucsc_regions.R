# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(stringr)
  library(reshape2)
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(circlize))
  library(qs)
}

cell1 <- "thp1"
cell2 <- "cd34"

# read the states data
if (T) {
  state_file <- "data/raw_data/2.states/chromIDEAS.state"
  
  state <- data.table::fread(state_file, sep = " ", header = T, data.table = F)
  state <- state[, c(1, 5, 6)]
  colnames(state)[1] <- "ID"
  head(state)
}

# define functions
if (T) {
  # cell specific state preference for genomic regions
  state_preference_region <- function(cell, state, genomics) {
    # test data
    if (F) {
      cell <- cell1
    }
    
    # get state
    if (T) {
      cell_state <- split(state$ID, state[, cell])
    }
    
    # format the genomic regions
    if (T) {
      genomics$types <- paste0(genomics$data_source, "@", genomics$Type)
      all_types <- sort(unique(genomics$types))
    }
    
    # calculate the log2fc to represent the state preference
    if (T) {
      prefer <- lapply(names(cell_state), function(state_n) {
        # test data
        if (F) {
          state_n <- names(cell_state)[10]
        }
        
        # info print
        if (T) {
          cat(paste0("Now process S", state_n, ":\n"))
        }
        
        # calculate the preference to each region for the state
        if (T) {
          sub_state <- cell_state[[state_n]]
          
          region_prefer <- lapply(all_types, function(region) {
            # region <- all_types[1]
            mess <- paste0("\t(", which(all_types == region), "/", length(all_types), ") ", region, "\n")
            cat(mess)
            
            region_bins <- genomics$ID[genomics$types == region]
            
            # All DEGs
            n <- length(sub_state)
            
            # All genes
            N <- nrow(state)
            
            # M: All genes in TermA
            M <- length(region_bins)
            
            # k: DEGs in TermA
            k <- length(intersect(sub_state, region_bins))
            
            # enrichment P: phyper(k-1,M,N-M,n,lower.tail = F)
            p <- phyper(k-1, M, N-M, n, lower.tail = F)
            
            # fold change
            term_in_DEG <- k / n
            term_in_BG <- M / N
            fc <- term_in_DEG / term_in_BG
            
            res <- c(p, N, M, n, k, term_in_DEG, term_in_BG, fc, region)
            names(res) <- c("P_hyper", "All_genes_N", "All_genes_in_term_M", "All_DEGs_n", 
                            "DEGs_in_Term_k", "term_in_DEG", "term_in_BG", "FC", "region")
            
            return(res)
          })
        }
        
        # merge the data
        if (T) {
          region_prefer <- data.frame(do.call(rbind, region_prefer))
          region_prefer$State <- state_n
        }
        
        return(region_prefer)
      })
    }
    
    # merge and format the data
    if (T) {
      prefer <- data.frame(do.call(rbind, prefer))
      prefer[, 1:8] <- sapply(prefer[, 1:8], as.numeric)
      prefer$data_source <- str_split(prefer$region, "@", simplify = T)[, 1]
      prefer$regions <- str_split(prefer$region, "@", simplify = T)[, 2]
    }
    
    return(prefer)
  }
  
  state_preference_region2 <- function(state, genomics) {
    # get state
    if (T) {
      cell_state1 <- split(state$ID, state[, cell1])
      cell_state2 <- split(state$ID, state[, cell2])
    }
    
    # format the genomic regions
    if (T) {
      genomics$types <- paste0(genomics$data_source, "@", genomics$Type)
      all_types <- sort(unique(genomics$types))
    }
    
    # calculate the log2fc to represent the state preference
    if (T) {
      if (identical(names(cell_state1), names(cell_state2))) {
        all_states <- names(cell_state1)
      }

      prefer <- lapply(all_states, function(state_n) {
        # test data
        if (F) {
          state_n <- all_states[10]
        }
        
        # info print
        if (T) {
          cat(paste0("Now process S", state_n, ":\n"))
        }
        
        # calculate the preference to each region for the state
        if (T) {
          sub_state1 <- cell_state1[[state_n]]
          sub_state2 <- cell_state2[[state_n]]
          
          region_prefer <- lapply(all_types, function(region) {
            # region <- all_types[1]
            mess <- paste0("\t(", which(all_types == region), "/", length(all_types), ") ", region, "\n")
            cat(mess)
            
            region_bins <- genomics$ID[genomics$types == region]
            
            # All DEGs
            n <- length(sub_state1) + length(sub_state2)
            
            # All genes
            N <- nrow(state) * 2
            
            # M: All genes in TermA
            M <- length(region_bins) * 2
            
            # k: DEGs in TermA
            k <- length(intersect(sub_state1, region_bins)) + length(intersect(sub_state2, region_bins))
            
            # enrichment P: phyper(k-1,M,N-M,n,lower.tail = F)
            p <- phyper(k-1, M, N-M, n, lower.tail = F)
            
            # fold change
            term_in_DEG <- k / n
            term_in_BG <- M / N
            fc <- term_in_DEG / term_in_BG
            
            res <- c(p, N, M, n, k, term_in_DEG, term_in_BG, fc, region)
            names(res) <- c("P_hyper", "All_genes_N", "All_genes_in_term_M", "All_DEGs_n", 
                            "DEGs_in_Term_k", "term_in_DEG", "term_in_BG", "FC", "region")
            
            return(res)
          })
        }
        
        # merge the data
        if (T) {
          region_prefer <- data.frame(do.call(rbind, region_prefer))
          region_prefer$State <- state_n
        }
        
        return(region_prefer)
      })
    }
    
    # merge and format the data
    if (T) {
      prefer <- data.frame(do.call(rbind, prefer))
      prefer[, 1:8] <- sapply(prefer[, 1:8], as.numeric)
      prefer$data_source <- str_split(prefer$region, "@", simplify = T)[, 1]
      prefer$regions <- str_split(prefer$region, "@", simplify = T)[, 2]
    }
    
    return(prefer)
  }
}

# load the genomic region bin ID
if (T) {
  file <- "data/saved_data/10.genomic_region_bin_IDs.qs"
  if (! file.exists(file)) {
    gtfs <- list.files(path = "data/saved_data/", pattern = "8.common_gtf_structure.*.qs")
    ucscs <- list.files(path = "data/saved_data/", pattern = "9.common_ucsc_structure.*.qs")
    files <- c(gtfs, ucscs)
    genomics <- lapply(files, function(x) {
      # x <- files[10]
      fullpath <- paste0("data/saved_data/", x)
      dat <- qread(fullpath, nthreads = 6)
      type <- str_split_1(x, "_")[4]
      
      if (is.vector(dat)) {
        dat <- data.frame(ID = dat, 
                          Type = type)
      } else {
        dat <- data.frame(ID = sort(unique(as.numeric(dat$ID))), 
                          Type = type)
      }
      dat$data_source <- str_split_1(x, "_")[2]
      
      return(dat)
    })
    genomics <- data.frame(do.call(rbind, genomics))
    
    rm(gtfs, ucscs, files)
    
    qsave(genomics, file = file, nthreads = 6)
  }
  if (file.exists(file)) {
    genomics <- qread(file, nthreads = 6)
  }
  
  table(genomics$Type, genomics$data_source)
}

# calculate the data
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    for (cell in c(cell1, cell2)) {
      cat(paste0("\t", cell, ": \n"))
      
      # test data
      if (F) {
        cell <- cell1
        data_type <- "tx"
      }
      
      # load cell speicific state distribution
      if (T) {
        file <- paste0("data/saved_data/10.state_prefer_genomic_region_stat_", cell, ".qs")
        
        if (! file.exists(file)) {
          cell_state_prefer <- state_preference_region(cell, state, genomics)
          
          qsave(cell_state_prefer, file = file, nthreads = 6)
        }
      }
    }
    
    cat(paste0("\tMerged dat: \n"))
    
    # load cell speicific state distribution
    if (T) {
      file <- paste0("data/saved_data/10.state_prefer_genomic_region_stat_merged.qs")
      
      if (! file.exists(file)) {
        cell_state_prefer <- state_preference_region2(state, genomics)
        
        qsave(cell_state_prefer, file = file, nthreads = 6)
      }
    }
  }
}

# plot the figure
if (T) {
  for (data_type in c("tx")) {
    for (cell in c(cell1, cell2, "merged")) {
      # test data
      if (F) {
        data_type <- "tx"
        cell <- "merged"
      }
      
      # load cell speicific state distribution
      if (T) {
        file <- paste0("data/saved_data/10.state_prefer_genomic_region_stat_", cell, ".qs")
        cell_state_prefer <- qread(file, nthreads = 6)
      }
      
      # read the ordering
      if (T) {
        ord <- qread("data/saved_data/2.order_of_emission_table.qs", nthreads = 6)
        ord <- ord[[1]]
      }
      
      # format the data
      if (T) {
        head(cell_state_prefer)
        
        # convert into log2fc
        if (T) {
          non0min <- cell_state_prefer$FC
          non0min <- min(non0min[non0min>0])
          
          cell_state_prefer$Log2FC <- log2(cell_state_prefer$FC + non0min)
        }
        
        # long2wide
        if (T) {
          ggdat <- dcast(cell_state_prefer, State~region, value.var = "Log2FC")
          rownames(ggdat) <- paste0("S", ggdat[, 1])
          ggdat <- ggdat[, -1]
        }
        
        # format the colnames
        if (T) {
          colnames(ggdat) <- str_split(colnames(ggdat), "@", simplify = T)[, 2]
          colnames(ggdat) <- ifelse(colnames(ggdat) %in% c("cds", "tes", "tss", "utr3", "utr5"), toupper(colnames(ggdat)), 
                                    ifelse(colnames(ggdat) %in% c("exon", "intron", "intergenic", "repeats"), str_to_title(colnames(ggdat)), "CpG"))
          ggdat <- ggdat[, c("CpG", "TSS", "UTR5", "UTR3", "TES", "Exon", "CDS", "Intron", "Intergenic", "Repeats")]
        }
        
        # sorting
        if (T) {
          ggdat <- ggdat[ord, ]
        }
      }
      
      # heatmap
      if (T) {
        # color setting
        if (T) {
          library(RColorBrewer)
          colors <- brewer.pal(9, "Set1") 
          
          colors <- colors[1:2]
          col_fun = colorRamp2(c(ceiling(quantile(c(as.matrix(ggdat)), 0.05)), 0, ceiling(quantile(c(as.matrix(ggdat)), 0.95))),
                               c(colors[2], "white", colors[1]))
        }
        
        # heatmap
        if (T) {
          p <- Heatmap(as.matrix(ggdat), name = "Enrichment Level", col = col_fun, 
                       border_gp = gpar(col = "black"), rect_gp = gpar(col = "black"),
                       row_split = NULL, cluster_row_slices = FALSE,
                       show_heatmap_legend = TRUE, heatmap_legend_param = list(border = T, title_position = "lefttop", legend_direction = "horizontal"),
                       row_names_max_width = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), row_names_rot = 0,
                       column_names_max_height = unit(6, "cm"), column_names_gp = gpar(fontsize = 10), column_names_rot = 60,
                       cluster_rows = F, cluster_columns = F, show_row_dend = F, show_column_dend = F, 
                       clustering_distance_columns = "euclidean", clustering_distance_rows = "euclidean", 
                       clustering_method_columns = "ward.D2", clustering_method_rows = "ward.D2")
        }
        
        # save the figure
        if (T) {
          p <- ggplotify::as.ggplot(p)
          file <- paste0("results/2.pic/10.cell_specific_state_prefer_regions_heatmap.", cell, ".pdf")
          ggsave(filename = file, plot = p, width = 5.5, height = 10)
        }
      }
    }
  }
}