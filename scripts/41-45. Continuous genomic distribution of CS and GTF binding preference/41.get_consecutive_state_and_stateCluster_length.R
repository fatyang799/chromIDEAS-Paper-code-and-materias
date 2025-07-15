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
state_file <- "data/raw_data/2.states/chromIDEAS.state"

# read the states data
if (T) {
  state <- data.table::fread(state_file, sep = " ", header = T, data.table = F)
  state <- state[, c(1, 5, 6)]
  colnames(state)[1] <- "ID"
  head(state)
}

# define functions
if (T) {
  # get gene CS
  cell_geneCS <- function(genebody_stat, state) {
    # test dat
    if (F) {
      state <- state[, c("ID", cell1)]
    }
    
    # get gene specific state pattern
    if (T) {
      t_genebody_stat <- data.frame(t(genebody_stat))
      
      gene_state <- lapply(t_genebody_stat, function(stat) {
        # stat <- t_genebody_stat[, 1]
        strand <- stat[2]
        tss <- as.numeric(stat[3])
        tes <- as.numeric(stat[4])
        
        tss_up <- ifelse(strand == "+", tss - up_bin_num, 
                         ifelse(strand == "-", tss + up_bin_num, NA))
        tes_down <- ifelse(strand == "+", tes + down_bin_num, 
                           ifelse(strand == "-", tes - down_bin_num, NA))
        
        ids <- (tss_up:tes_down)
        states <- state[ids, 2]
        states <- as.numeric(states)
        
        return(states)
      })
      names(gene_state) <- genebody_stat[, 1]
    }
    
    return(gene_state)
  }
  
  # get consecutive stat
  consecutive_seq <- function(cell_geneCS_dat) {
    # prepare hello info
    if (T) {
      start_mess <- paste0("|", paste(rep("-", 100), collapse = ""), "|\n")
      cat(start_mess)
      
      breaks <- round(seq(1, length(cell_geneCS_dat), length.out=100))
      breaks <- names(cell_geneCS_dat)[breaks]
    }
    
    # get most consecutive state for each gene data
    if (T) {
      mat <- data.frame(t(
        sapply(names(cell_geneCS_dat), function(id) {
          # id <- names(cell_geneCS_dat)[1]
          # print process
          if (id %in% breaks) {
            if (which(id == breaks) == 1) {
              cat("|*")
            }
            if (which(id == breaks) == 100) {
              cat("*|\n")
            }
            if (! which(id == breaks) %in% c(1, 100)) {
              cat("*")
            }
          }
          
          # Compute the lengths for each consecutive values
          if (T) {
            gene_state <- cell_geneCS_dat[[id]]
            
            rle_values <- rle(gene_state)
            max_counts <- tapply(rle_values$lengths, rle_values$values, max)
            names(max_counts) <- paste0("S", names(max_counts))
          }
          
          # stat
          if (T) {
            state_with_max_n <- max(max_counts)
            state_with_max_name <- names(max_counts)[as.numeric(max_counts) == state_with_max_n]
            state_with_max_name <- paste(state_with_max_name, collapse = "+")
          }
          
          res <- c("most_state_name" = state_with_max_name, 
                   "most_state_number" = state_with_max_n)
          
          return(res)
        })
      ))
      mat$gene_id <- names(cell_geneCS_dat)
    }
    
    return(mat)
  }
  
  # replace CS with CScluster
  replace_CS_with_cluster <- function(cell_geneCS_dat, cluster_dat) {
    cell_geneCS_cluster <- lapply(names(cell_geneCS_dat), function(id) {
      # id <- names(cell_geneCS_dat)[1]
      states <- cell_geneCS_dat[[id]]
      cluster <- cluster_dat[match(states, cluster_dat$state), "cluster"]
      
      return(cluster)
    })
    names(cell_geneCS_cluster) <- names(cell_geneCS_dat)
    
    return(cell_geneCS_cluster)
  }
}

# get GSpat
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
      
      # get gene body location and filter by length>=3
      if (T) {
        file <- paste0("data/saved_data/11.", data_type, "_level_Body_region_length_ge3_bins.qs")
        genebody_stat <- qread(file, nthreads = 6)
        
        genebody_stat <- genebody_stat[, c(1, 4, 2, 3)]
        head(genebody_stat)
      }
      
      # mkdir
      if (T) {
        dir <- paste0("data/saved_data/41.gene_continuous_state/all_txs/", data_type)
        if (! dir.exists(dir)) {
          dir.create(dir, showWarnings = F, recursive = T)
        }
      }
      
      # get cell specific each gene state pattern: cell_geneCS_dat
      if (T) {
        file <- paste0(dir, "/41.", data_type, "_level_profile_Body_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, 
                       "Body.bs", bin_size, ".states_across_geneBody_UpDown.", cell, ".qs")
        
        if (file.exists(file)) {
          cell_geneCS_dat <- qread(file, nthreads = 6)
        }
        if (! file.exists(file)) {
          cell_geneCS_dat <- cell_geneCS(genebody_stat, state[, c("ID", cell)])
          
          qsave(cell_geneCS_dat, file, nthreads = 6)
        }
      }
      
      # each gene continuous chromatin state number statistics: gene_continuous_state
      if (T) {
        file <- paste0(dir, "/41.", data_type, "_level_profile_Body_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, 
                       "Body.bs", bin_size, ".continuous_max_state.", cell, ".qs")
        if (! file.exists(file)) {
          cat(paste0("\t\tNow get each gene continuous chromatin state number statistics: \n"))
          
          gene_continuous_state <- consecutive_seq(cell_geneCS_dat)
          qsave(gene_continuous_state, file = file, nthreads = 6)
        }
      }
    }
    
    # for merged data (in paper)
    if (T) {
      file <- paste0(dir, "/41.", data_type, "_level_profile_Body_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, 
                     "Body.bs", bin_size, ".continuous_max_state.merged.qs")
      if (! file.exists(file)) {
        # load cell specific data
        if (T) {
          file1 <- paste0(dir, "/41.", data_type, "_level_profile_Body_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, 
                          "Body.bs", bin_size, ".continuous_max_state.", cell1 ,".qs")
          dat1 <- qread(file1, nthreads = 6)
          file2 <- paste0(dir, "/41.", data_type, "_level_profile_Body_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, 
                          "Body.bs", bin_size, ".continuous_max_state.", cell2 ,".qs")
          dat2 <- qread(file2, nthreads = 6)
        }
        
        # make the order identical
        if (! identical(dat1$gene_id, dat2$gene_id)) {
          dat2 <- dat2[match(dat1$gene_id, dat2$gene_id), ]
        }
        
        # get max as representative
        if (T) {
          gene_continuous_state <- data.frame(
            most_state_name = ifelse(dat1$most_state_name == dat2$most_state_name, 
                                     dat1$most_state_name, 
                                     paste0(dat1$most_state_name, "@", dat2$most_state_name)), 
            most_state_number = sapply(1:nrow(dat1), function(x) {
              max(c(as.numeric(dat1$most_state_number[x]), as.numeric(dat2$most_state_number[x])))
            }), 
            gene_id = dat1$gene_id
          )
        }
        
        qsave(gene_continuous_state, file = file, nthreads = 6)
      }
    }
  }
}

# replace the CS with CS cluster
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir <- paste0("data/saved_data/41.gene_continuous_state/all_txs/", data_type)
      if (! dir.exists(dir)) {
        dir.create(dir, showWarnings = F, recursive = T)
      }
    }
    
    for (type in c("merged", "single")) {
      cat(paste0("\t", type, ": \n"))
      
      for (cell in c(cell1, cell2)) {
        cat(paste0("\t\t", cell, ": \n"))
        
        for (cluster in c("chromIDEAS", "kmeans", 1:n_ct)) {
          cat(paste0("\t\t\t", cluster, ": \n"))
          
          # test data
          if (F) {
            data_type <- "tx"
            type <- "merged"
            cell <- cell1
            cluster <- "kmeans"
          }
          
          # read the cluster: cluster_dat
          if (T) {
            if (type == "merged") {
              if (cluster == "chromIDEAS") {
                cluster_dat <- qread("data/saved_data/29.common_cluster_name_for_merged_2cells_CS_clusters.qs", nthreads = 6)
                
                # only retain the chromatin state and state cluster 2 columns data
                cluster_dat <- cluster_dat[, c("orig.ident", "seurat_clusters")]
                colnames(cluster_dat) <- c("state", "cluster")
              } else if (cluster %in% as.character(1:n_ct)) {
                file <- paste0("data/saved_data/40.random_cluster_control_n", n_ct, ".merged.qs")
                cluster_dat <- qread(file, nthreads = 6)
                cluster_dat <- cluster_dat[, c("state", paste0("ct", cluster))]
                colnames(cluster_dat) <- c("state", "cluster")
              }
            } else if (type == "single") {
              if (cluster == "chromIDEAS") {
                cluster_dat <- qread("data/saved_data/25.common_cluster_name_for_2_cell_CS_clusters.qs", nthreads = 6)
                cluster_dat$orig.ident <- paste0("S", cluster_dat$orig.ident)
                cluster_dat <- cluster_dat[cluster_dat$cell == cell, ]
                
                # only retain the chromatin state and state cluster 2 columns data
                cluster_dat <- cluster_dat[, c("orig.ident", "seurat_clusters")]
                colnames(cluster_dat) <- c("state", "cluster")
              } else if (cluster %in% as.character(1:n_ct)) {
                file <- paste0("data/saved_data/40.random_cluster_control_n", n_ct, ".", cell, ".qs")
                cluster_dat <- qread(file, nthreads = 6)
                cluster_dat <- cluster_dat[, c("state", paste0("ct", cluster))]
                colnames(cluster_dat) <- c("state", "cluster")
              }
            }
            if (cluster == "kmeans") {
              cluster_dat <- read.table("results/1.tab/40.emission_kmeans_cluster.csv", header = T, sep = ",")
              
              # only retain the chromatin state and state cluster 2 columns data
              cluster_dat <- cluster_dat[, c("state", "cluster5")]
              colnames(cluster_dat) <- c("state", "cluster")
            }
            
            # format the data
            if (T) {
              cluster_dat$state <- as.character(cluster_dat$state)
              cluster_dat$state <- as.numeric(gsub("S", "", cluster_dat$state))
              cluster_dat$cluster <- as.numeric(as.vector(cluster_dat$cluster))
            }
          }
          
          # read the GSpat data
          if (T) {
            file <- paste0(dir, "/41.", data_type, "_level_profile_Body_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, 
                           "Body.bs", bin_size, ".states_across_geneBody_UpDown.", cell, ".qs")
            cell_geneCS_dat <- qread(file, nthreads = 6)
          }
          
          # replace the chromatin state with state cluster: cell_geneCS_cluster
          if (T) {
            file <- paste0(dir, "/41.", data_type, "_Body_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, 
                           "Body.bs", bin_size, ".", type, ".clusters(", cluster, ")_Body_UpDown.", cell, ".qs")
            
            if (file.exists(file)) {
              cell_geneCS_cluster <- qread(file, nthreads = 6)
            }
            if (! file.exists(file)) {
              cell_geneCS_cluster <- replace_CS_with_cluster(cell_geneCS_dat, cluster_dat)
              
              qsave(cell_geneCS_cluster, file, nthreads = 6)
            }
          }
          
          # each gene continuous cluster number statistics
          if (T) {
            file <- paste0(dir, "/41.", data_type, "_Body_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, 
                           "Body.bs", bin_size, ".", type, ".continuous_max_cluster(", cluster, ").", cell, ".qs")
            
            if (! file.exists(file)) {
              cat(paste0("\t\tNow get each gene continuous cluster number statistics: \n"))
              
              gene_continuous_cluster <- consecutive_seq(cell_geneCS_cluster)
              qsave(gene_continuous_cluster, file = file, nthreads = 6)
              
              rm(gene_continuous_cluster)
            }
          }
        }
      }
    }
    
    # for merged data (in paper)
    if (T) {
      for (cluster in c("chromIDEAS", "kmeans", 1:n_ct)) {
        # test data
        if (F) {
          data_type <- "tx"
          cluster <- "kmeans"
        }
        
        # load cell data
        file <- paste0(dir, "/41.", data_type, "_Body_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, 
                       "Body.bs", bin_size, ".merged.continuous_max_cluster(", cluster, ").merged.qs")
        
        if (! file.exists(file)) {
          # load cell specific data
          if (T) {
            file1 <- paste0(dir, "/41.", data_type, "_Body_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, 
                            "Body.bs", bin_size, ".merged.continuous_max_cluster(", cluster, ").", cell1, ".qs")
            file2 <- paste0(dir, "/41.", data_type, "_Body_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, 
                            "Body.bs", bin_size, ".merged.continuous_max_cluster(", cluster, ").", cell2, ".qs")
            dat1 <- qread(file1, nthreads = 6)
            dat2 <- qread(file2, nthreads = 6)
          }
          
          # make the order identical
          if (! identical(dat1$gene_id, dat2$gene_id)) {
            dat2 <- dat2[match(dat1$gene_id, dat2$gene_id), ]
          }
          
          # get max as representative
          if (T) {
            dat <- data.frame(
              most_state_name = ifelse(dat1$most_state_name == dat2$most_state_name, 
                                       dat1$most_state_name, 
                                       paste0(dat1$most_state_name, "@", dat2$most_state_name)), 
              most_state_number = sapply(1:nrow(dat1), function(x) {
                max(c(as.numeric(dat1$most_state_number[x]), as.numeric(dat2$most_state_number[x])))
              }), 
              gene_id = dat1$gene_id
            )
          }
          
          qsave(dat, file = file, nthreads = 6)
        }
      }
    }
  }
}

# plot the length statistics figure
if (T) {
  # load the environment
  if (T) {
    rm(list = ls())
    options(stringAsFactors = F)
    library(ggplot2)
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
    n_ct <- 5
  }
  
  # common setting
  if (T) {
    data_type <- "tx"
    dir <- paste0("data/saved_data/41.gene_continuous_state/all_txs/", data_type)
  }
  
  for (type in c("merged", "single")) {
    cat(paste0(type, ": \n"))
    
    # cell specific data
    for (cell in c(cell1, cell2)) {
      cat(paste0("\t", cell, ": \n"))
      
      # test data
      if (F) {
        type <- "merged"
        cell <- cell1
      }
      
      # load state continuous length
      if (T) {
        file <- paste0(dir, "/41.", data_type, "_level_profile_Body_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, 
                       "Body.bs", bin_size, ".continuous_max_state.", cell, ".qs")
        gene_continuous_state <- qread(file, nthreads = 6)
      }
      
      # load cluster continuous length
      if (T) {
        for (cluster in c("chromIDEAS", "kmeans", 1:n_ct)) {
          # cluster <- "chromIDEAS"
          file <- paste0(dir, "/41.", data_type, "_Body_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, 
                         "Body.bs", bin_size, ".", type, ".continuous_max_cluster(", cluster, ").", cell, ".qs")
          dat <- qread(file, nthreads = 6)
          
          # assign data
          assign(paste0("gene_continuous_", cluster), dat)
          
          rm(dat, file)
        }
      }
      
      # check the data
      if (T) {
        # gene id
        torf1 <- identical(gene_continuous_state$gene_id, gene_continuous_chromIDEAS$gene_id)
        torf2 <- identical(gene_continuous_state$gene_id, gene_continuous_kmeans$gene_id)
        
        torf3 <- identical(gene_continuous_state$gene_id, gene_continuous_1$gene_id)
        torf4 <- identical(gene_continuous_state$gene_id, gene_continuous_2$gene_id)
        torf5 <- identical(gene_continuous_state$gene_id, gene_continuous_3$gene_id)
        torf6 <- identical(gene_continuous_state$gene_id, gene_continuous_4$gene_id)
        torf7 <- identical(gene_continuous_state$gene_id, gene_continuous_5$gene_id)
        
        if (torf1+torf2+torf3+torf4+torf5+torf6+torf7 < 7) {
          stop("error1")
        }
        
        rm(list = paste0("torf", 1:7))
      }
      
      # summary the length results
      if (T) {
        dat <- data.frame(id = gene_continuous_state$gene_id, 
                          raw_state = as.numeric(gene_continuous_state$most_state_number), 
                          chromideas = as.numeric(gene_continuous_chromIDEAS$most_state_number), 
                          kmeans = as.numeric(gene_continuous_kmeans$most_state_number), 
                          
                          ct1 = as.numeric(gene_continuous_1$most_state_number), 
                          ct2 = as.numeric(gene_continuous_2$most_state_number), 
                          ct3 = as.numeric(gene_continuous_3$most_state_number), 
                          ct4 = as.numeric(gene_continuous_4$most_state_number), 
                          ct5 = as.numeric(gene_continuous_5$most_state_number))
        
        rm(gene_continuous_state, gene_continuous_chromIDEAS, gene_continuous_kmeans, 
           gene_continuous_1, gene_continuous_2, gene_continuous_3, gene_continuous_4, gene_continuous_5)
      }
      
      # format the data
      if (T) {
        ggdat <- melt(dat, id.vars = "id", variable.name = "cluster_method", value.name = "Len")
        ggdat$cluster_method <- factor(ggdat$cluster_method, 
                                       levels = c("kmeans", "chromideas", "raw_state", paste0("ct", 1:5)), 
                                       labels = c("Kmeans", "chromIDEAS", "RawState", paste0("CT", 1:5)))
        
        print(t(do.call(rbind, tapply(ggdat$Len, ggdat$cluster_method, summary))))
      }
      
      # ggplot
      if (T) {
        # statistis group
        if (T) {
          ref <- "chromIDEAS"
          remain <- setdiff(levels(ggdat$cluster_method), ref)
          my_comparisons <- lapply(remain, function(x) {
            c(ref, x)
          })
        }
        
        p <- ggplot(ggdat, aes(x=cluster_method, y=log2(Len))) +
          geom_violin(aes(fill=cluster_method, group=cluster_method), scale = "width") +
          geom_boxplot(aes(group=cluster_method), width=0.1, fill=NA) +
          ggpubr::stat_compare_means(method="wilcox.test", paired = T, label="p.signif", comparisons = my_comparisons) +
          ggpubr::stat_compare_means(method="anova", label.y=20) +
          xlab(NULL) +
          ylab("log2(Length of longest consecutive sequence)") +
          ggtitle(paste0(cell, "+", data_type)) +
          cowplot::theme_cowplot() +
          theme(legend.position = "none",
                strip.background = element_rect(fill = NA))
      }
      
      file <- paste0("results/2.pic/41.longest_consecutive_seq_in_", data_type, "s_profile_Body.", cell, ".", type, ".pdf")
      ggsave(filename = file, plot = p, width = 12, height = 8)
    }
    
    # merged data
    if (type == "merged") {
      cat(paste0("\tmerged: \n"))
      
      # define function
      if (T) {
        load_cell_merged_dat <- function(prefix, cell) {
          file <- paste0(prefix, cell ,".qs")
          dat <- qread(file, nthreads = 6)
          return(dat)
        }
        merge_cell_dat <- function(dat1, dat2) {
          # make the order identical
          if (! identical(dat1$gene_id, dat2$gene_id)) {
            dat2 <- dat2[match(dat1$gene_id, dat2$gene_id), ]
          }
          
          # get max as representative
          if (T) {
            gene_continuous_state <- data.frame(
              most_state_name = ifelse(dat1$most_state_name == dat2$most_state_name, 
                                       dat1$most_state_name, 
                                       paste0(dat1$most_state_name, "@", dat2$most_state_name)), 
              most_state_number = sapply(1:nrow(dat1), function(x) {
                max(c(as.numeric(dat1$most_state_number[x]), as.numeric(dat2$most_state_number[x])))
              }), 
              gene_id = dat1$gene_id
            )
          }
          
          return(gene_continuous_state)
        }
      }
      
      # load state continuous length
      if (T) {
        # load cell specific data
        if (T) {
          prefix <- paste0(dir, "/41.", data_type, "_level_profile_Body_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, 
                          "Body.bs", bin_size, ".continuous_max_state.")
          dat1 <- load_cell_merged_dat(prefix, cell1)
          dat2 <- load_cell_merged_dat(prefix, cell2)
        }
        
        # get max as representative
        gene_continuous_state <- merge_cell_dat(dat1, dat2)
      }
      
      # load cluster continuous length
      if (T) {
        for (cluster in c("chromIDEAS", "kmeans", 1:n_ct)) {
          # cluster <- "chromIDEAS"
          
          # load cell specific data
          if (T) {
            prefix <- paste0(dir, "/41.", data_type, "_Body_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, 
                             "Body.bs", bin_size, ".merged.continuous_max_cluster(", cluster, ").")
            dat1 <- load_cell_merged_dat(prefix, cell1)
            dat2 <- load_cell_merged_dat(prefix, cell2)
          }
          
          # get max as representative
          dat <- merge_cell_dat(dat1, dat2)
          
          # assign data
          assign(paste0("gene_continuous_", cluster), dat)
          
          rm(dat1, dat2, dat, prefix)
        }
      }
      
      # check the data
      if (T) {
        # gene id
        torf1 <- identical(gene_continuous_state$gene_id, gene_continuous_chromIDEAS$gene_id)
        torf2 <- identical(gene_continuous_state$gene_id, gene_continuous_kmeans$gene_id)
        
        torf3 <- identical(gene_continuous_state$gene_id, gene_continuous_1$gene_id)
        torf4 <- identical(gene_continuous_state$gene_id, gene_continuous_2$gene_id)
        torf5 <- identical(gene_continuous_state$gene_id, gene_continuous_3$gene_id)
        torf6 <- identical(gene_continuous_state$gene_id, gene_continuous_4$gene_id)
        torf7 <- identical(gene_continuous_state$gene_id, gene_continuous_5$gene_id)
        
        if (torf1+torf2+torf3+torf4+torf5+torf6+torf7 < 7) {
          stop("error1")
        }
        
        rm(list = paste0("torf", 1:7))
      }
      
      # summary the length results
      if (T) {
        dat <- data.frame(id = gene_continuous_state$gene_id, 
                          raw_state = as.numeric(gene_continuous_state$most_state_number), 
                          chromideas = as.numeric(gene_continuous_chromIDEAS$most_state_number), 
                          kmeans = as.numeric(gene_continuous_kmeans$most_state_number), 
                          
                          ct1 = as.numeric(gene_continuous_1$most_state_number), 
                          ct2 = as.numeric(gene_continuous_2$most_state_number), 
                          ct3 = as.numeric(gene_continuous_3$most_state_number), 
                          ct4 = as.numeric(gene_continuous_4$most_state_number), 
                          ct5 = as.numeric(gene_continuous_5$most_state_number))
        
        rm(gene_continuous_state, gene_continuous_chromIDEAS, gene_continuous_kmeans, 
           gene_continuous_1, gene_continuous_2, gene_continuous_3, gene_continuous_4, gene_continuous_5)
      }
      
      # format the data
      if (T) {
        ggdat <- melt(dat, id.vars = "id", variable.name = "cluster_method", value.name = "Len")
        ggdat$cluster_method <- factor(ggdat$cluster_method, 
                                       levels = c("kmeans", "chromideas", "raw_state", paste0("ct", 1:5)), 
                                       labels = c("Kmeans", "chromIDEAS", "RawState", paste0("CT", 1:5)))
        
        print(t(do.call(rbind, tapply(ggdat$Len, ggdat$cluster_method, summary))))
      }
      
      # ggplot
      if (T) {
        # statistis group
        if (T) {
          ref <- "chromIDEAS"
          remain <- setdiff(levels(ggdat$cluster_method), ref)
          my_comparisons <- lapply(remain, function(x) {
            c(ref, x)
          })
        }
        
        p <- ggplot(ggdat, aes(x=cluster_method, y=log2(Len))) +
          geom_violin(aes(fill=cluster_method, group=cluster_method), scale = "width") +
          geom_boxplot(aes(group=cluster_method), width=0.1, fill=NA) +
          ggpubr::stat_compare_means(method="wilcox.test", paired = T, label="p.signif", comparisons = my_comparisons) +
          ggpubr::stat_compare_means(method="anova", label.y=20) +
          xlab(NULL) +
          ylab("log2(Length of longest consecutive sequence)") +
          ggtitle(paste0("Merged+", data_type)) +
          cowplot::theme_cowplot() +
          theme(legend.position = "none",
                strip.background = element_rect(fill = NA))
      }
      
      file <- paste0("results/2.pic/41.longest_consecutive_seq_in_", data_type, "s_profile_Body.merged.merged.pdf")
      ggsave(filename = file, plot = p, width = 12, height = 8)
    }
  }
}

## merged: 
##     thp1: 
##              Kmeans chromIDEAS  RawState       CT1       CT2       CT3        CT4       CT5
## Min.        2.00000    2.00000   2.00000   2.00000   2.00000   2.00000    2.00000   2.00000
## 1st Qu.    15.00000   14.00000   8.00000   9.00000   9.00000   9.00000   11.00000  10.00000
## Median     29.00000   27.00000  15.00000  17.00000  16.00000  16.00000   21.00000  18.00000
## Mean       92.53389   70.79383  22.84314  26.13375  24.79535  25.05306   47.24888  27.63311
## 3rd Qu.    71.00000   61.00000  29.00000  32.00000  31.00000  31.00000   50.00000  35.00000
## Max.    10350.00000 5180.00000 437.00000 571.00000 571.00000 571.00000 3442.00000 437.00000
##     cd34: 
##             Kmeans chromIDEAS RawState      CT1       CT2      CT3        CT4     CT5
## Min.        2.0000    2.00000   1.0000   2.0000   1.00000   1.0000    2.00000   2.000
## 1st Qu.    15.0000   14.00000   8.0000   9.0000   9.00000   9.0000   10.00000  10.000
## Median     29.0000   27.00000  15.0000  17.0000  16.00000  16.0000   19.00000  19.000
## Mean      101.5603   73.23068  24.6821  28.6264  27.72473  27.9674   46.19773  31.571
## 3rd Qu.    70.0000   61.00000  30.0000  36.0000  34.00000  35.0000   47.00000  38.000
## Max.    10350.0000 3442.00000 844.0000 844.0000 844.00000 844.0000 3442.00000 917.000
##             Kmeans chromIDEAS  RawState       CT1       CT2       CT3        CT4       CT5
## Min.        2.0000    2.00000   2.00000   2.00000   2.00000   2.00000    2.00000   2.00000
## 1st Qu.    17.0000   17.00000  10.00000  12.00000  11.00000  11.00000   13.00000  13.00000
## Median     34.0000   32.00000  19.00000  21.00000  20.00000  20.00000   25.00000  23.00000
## Mean      114.6199   85.88537  29.02316  32.88236  31.59622  31.86476   56.13664  36.03341
## 3rd Qu.    84.0000   73.00000  36.00000  41.00000  39.00000  40.00000   58.00000  44.00000
## Max.    10350.0000 5180.00000 844.00000 844.00000 844.00000 844.00000 3442.00000 917.00000
## single: 
##     thp1: 
##              Kmeans chromIDEAS  RawState       CT1       CT2       CT3        CT4      CT5
## Min.        2.00000    2.00000   2.00000   2.00000   2.00000   2.00000    2.00000   2.0000
## 1st Qu.    15.00000   14.00000   8.00000   9.00000   9.00000   9.00000   11.00000  10.0000
## Median     29.00000   26.00000  15.00000  17.00000  16.00000  16.00000   21.00000  18.0000
## Mean       92.53389   58.64925  22.84314  26.10677  24.78677  25.05235   47.24907  27.6261
## 3rd Qu.    71.00000   59.00000  29.00000  32.00000  31.00000  31.00000   50.00000  35.0000
## Max.    10350.00000 4210.00000 437.00000 571.00000 571.00000 571.00000 3442.00000 437.0000
##     cd34: 
##             Kmeans chromIDEAS RawState       CT1       CT2       CT3        CT4       CT5
## Min.        2.0000    2.00000   1.0000   2.00000   1.00000   1.00000    1.00000   2.00000
## 1st Qu.    15.0000   15.00000   8.0000   8.00000   8.00000   9.00000   10.00000  10.00000
## Median     29.0000   28.00000  15.0000  15.00000  15.00000  17.00000   19.00000  19.00000
## Mean      101.5603   91.29827  24.6821  25.65922  26.08637  28.36614   46.13997  31.65518
## 3rd Qu.    70.0000   66.00000  30.0000  31.00000  31.00000  35.00000   47.00000  38.00000
## Max.    10350.0000 9605.00000 844.0000 844.00000 844.00000 844.00000 3442.00000 917.00000
