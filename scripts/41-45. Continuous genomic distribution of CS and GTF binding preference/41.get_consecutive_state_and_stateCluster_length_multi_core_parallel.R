library(ggplot2)
library(stringr)
library(reshape2)
library(qs)
library(doParallel)
library(parallel)
library(foreach)
cl <- makeCluster(20, outfile="log41.txt")
registerDoParallel(cl)

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
up_bin_num <- 5
down_bin_num <- 5
body_bin_num <- 10
rna_levels <- 10
n_ct <- 5
state_file <- "data/raw_data/2.states/chromIDEAS.state"

# get GSpat
GSpat <- function(cell) {
  # load the environment
  if (T) {
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
    state_file <- "data/raw_data/2.states/chromIDEAS.state"
  }
  
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
  
  data_type <- "tx"
  
  # loop body
  if (T) {
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
}
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    dir <- paste0("data/saved_data/41.gene_continuous_state/all_txs/", data_type)
    
    # for cell specific statistics
    foreach(cell = c(cell1, cell2), .export =ls()) %dopar% GSpat(cell)
    
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
GTpat <- function(cluster) {
  # load the environment
  if (T) {
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
    state_file <- "data/raw_data/2.states/chromIDEAS.state"
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
  
  data_type <- "tx"
  
  # loop body
  if (T) {
    # mkdir
    if (T) {
      dir <- paste0("data/saved_data/41.gene_continuous_state/all_txs/", data_type)
      if (! dir.exists(dir)) {
        dir.create(dir, showWarnings = F, recursive = T)
      }
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
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    for (type in c("merged", "single")) {
      cat(paste0("\t", type, ": \n"))
      
      for (cell in c(cell1, cell2)) {
        cat(paste0("\t", cell, ": \n"))
        
        foreach(cluster = c("chromIDEAS", "kmeans", 1:n_ct), .export =ls()) %dopar% GTpat(cluster)
      }
    }
  }
}

stopCluster(cl)