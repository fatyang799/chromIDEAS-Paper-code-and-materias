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
up_bin_num <- 3
down_bin_num <- 3
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
  # get gene CS (modify based on script41)
  cell_geneCS <- function(genebody_stat, state, up_bin_num, down_bin_num) {
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
        
        # ids <- (tss:tes)
        ids <- (tss_up:tes_down)
        states <- state[ids, 2]
        states <- as.numeric(states)
        
        return(states)
      })
      names(gene_state) <- genebody_stat[, 1]
    }
    
    return(gene_state)
  }
  
  # replace CS with CScluster (script41)
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

# prepare the data
for (data_type in c("gene", "tx")) {
  cat(paste0(data_type, ": \n"))
  # test data
  if (F) {
    data_type <- "gene"
  }
  
  # get gene body location and filter by length>=3
  if (T) {
    file <- paste0("data/saved_data/11.", data_type, "_level_Body_region_length_ge3_bins.qs")
    genebody_stat <- qread(file, nthreads = 6)
    
    genebody_stat <- genebody_stat[, c(1, 4, 2, 3)]
    head(genebody_stat)
  }
  
  for (cell in c(cell1, cell2)) {
    cat(paste0("\t", cell, ": \n"))
    
    # test data
    if (F) {
      data_type <- "gene"
      cell <- cell1
    }
    
    # get cell specific each gene state pattern: cell_geneCS_dat
    if (T) {
      file <- paste0("data/saved_data/56.", data_type, "_level.bs", bin_size, ".states_across_geneBody_up", up_bin_num, "_down", down_bin_num, ".", cell, ".qs")
      
      if (file.exists(file)) {
        cell_geneCS_dat <- qread(file, nthreads = 6)
      }
      if (! file.exists(file)) {
        cell_geneCS_dat <- cell_geneCS(genebody_stat, state[, c("ID", cell)], up_bin_num, down_bin_num)
        
        qsave(cell_geneCS_dat, file, nthreads = 6)
      }
    }
    
    for (method in c("chromIDEAS", "kmeans")) {
      cat(paste0("\t\t", method, ": \n"))
      
      # test data
      if (F) {
        data_type <- "gene"
        method <- "chromIDEAS"
        cell <- cell1
      }
      
      # read the cluster: cluster_dat
      if (T) {
        if (method == "chromIDEAS") {
          cluster_dat <- qread("data/saved_data/29.common_cluster_name_for_merged_2cells_CS_clusters.qs", nthreads = 6)
          
          # only retain the chromatin state and state cluster 2 columns data
          cluster_dat <- cluster_dat[, c("orig.ident", "seurat_clusters")]
        }
        if (method == "kmeans") {
          cluster_dat <- read.table("results/1.tab/40.emission_kmeans_cluster.csv", header = T, sep = ",")
          cluster_dat <- cluster_dat[, c("state", "cluster5")]
        }
        
        # format the data
        if (T) {
          colnames(cluster_dat) <- c("state", "cluster")
          cluster_dat$state <- gsub("S", "", cluster_dat$state)
          cluster_dat$state <- as.numeric(as.vector(cluster_dat$state))
          cluster_dat$cluster <- as.numeric(as.vector(cluster_dat$cluster))
        }
      }
      
      # replace the chromatin state with state cluster: cell_geneCS_cluster
      if (T) {
        file <- paste0("data/saved_data/56.", data_type, "_level.bs", bin_size, ".", method, ".clusters_across_geneBody_up", up_bin_num, "_down", down_bin_num, ".", cell, ".qs")
        
        if (! file.exists(file)) {
          cell_geneCS_cluster <- replace_CS_with_cluster(cell_geneCS_dat, cluster_dat)
          
          qsave(cell_geneCS_cluster, file, nthreads = 6)
        }
      }
    }
  }
}

