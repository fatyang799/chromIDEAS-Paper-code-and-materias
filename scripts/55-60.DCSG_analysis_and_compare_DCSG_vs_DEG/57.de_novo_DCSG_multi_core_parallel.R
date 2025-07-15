library(stringr)
library(qs)
library(doParallel)
library(parallel)
library(foreach)
cl <- makeCluster(20, outfile="log57.txt")
registerDoParallel(cl)

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
up_bin_num <- 3
down_bin_num <- 3

# define functions
if (T) {
  # get gene CS
  get_gene_cs_dat <- function(data_type, cell) {
    # test data
    if (F) {
      cell <- cell1
    }
    
    file <- paste0("data/saved_data/56.", data_type, "_level.bs", bin_size, ".states_across_geneBody_up", up_bin_num, "_down", down_bin_num, ".", cell, ".qs")
    
    dat <- qread(file, nthreads = 6)
    
    return(dat)
  }
  get_gene_cluster_dat <- function(data_type, cell, method) {
    # test data
    if (F) {
      cell <- cell1
    }
    
    file <- paste0("data/saved_data/56.", data_type, "_level.bs", bin_size, ".", method, ".clusters_across_geneBody_up", up_bin_num, "_down", down_bin_num, ".", cell, ".qs")
    
    dat <- qread(file, nthreads = 6)
    
    return(dat)
  }
  
  # summary for dcsg
  dcsg_analysis <- function(cell1_dat, cell2_dat, dists) {
    # prepare hello info
    if (T) {
      start_mess <- paste0("|", paste(rep("-", 100), collapse = ""), "|\n")
      cat(start_mess)
      
      breaks <- round(seq(1, length(cell1_dat), length.out=100))
      breaks <- names(cell1_dat)[breaks]
    }
    
    # dcsg analysis
    if (T) {
      gene_dcsc <- lapply(names(cell1_dat), function(id) {
        # id <- names(cell1_dat)[1]
        # id <- "ENSG00000289544.1"
        
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
        
        # get cell specific state seq
        if (T) {
          dat1 <- cell1_dat[[id]]
          dat2 <- cell2_dat[[id]]
        }
        
        # summary the res
        if (T) {
          res <- data.frame(
            id = id, 
            location_id = 1:length(dat1), 
            location_p = (1:length(dat1)) / length(dat1), 
            cell1 = dat1, 
            cell2 = dat2
          )
          colnames(res) <- c("id", "location_id", "location_p", cell1, cell2)
        }
        
        # distance
        if (T) {
          res$distance <- sapply(1:nrow(res), function(x) {
            # x <- 1
            dists[paste0("S", res[, cell1][x]), paste0("S", res[, cell2][x])]
          })
        }
        
        return(res)
      })
      gene_dcsc <- do.call(rbind, gene_dcsc)
    }
    
    return(gene_dcsc)
  }
  
  # calculate the distance for dcsg cs cluster
  dcsg_cs_cluster_distance <- function(cs_dcsg, cluster_dcsg) {
    # get consistent gene id
    if (T) {
      ids <- unique(cluster_dcsg$id)
      cs_dcsg <- cs_dcsg[cs_dcsg$id %in% ids, ]
    }
    
    # calculate the distance
    if (T) {
      cs_dcsg_dist <- data.frame(
        id = paste0(cs_dcsg$id, cs_dcsg$dcs_loc_absolute), 
        dist = cs_dcsg$distance
      )
      cluster_dcsg_dist <- paste0(cluster_dcsg$id, cluster_dcsg$dcs_loc_absolute)
      
      cluster_dcsg$distance <- cs_dcsg_dist$dist[match(cluster_dcsg_dist, cs_dcsg_dist$id)]
    }
    
    return(cluster_dcsg)
  }
}

# data prepare
loop <- function(id) {
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
    up_bin_num <- 3
    down_bin_num <- 3
  }
  
  # define functions
  if (T) {
    # get gene CS
    get_gene_cs_dat <- function(data_type, cell) {
      # test data
      if (F) {
        cell <- cell1
      }
      
      file <- paste0("data/saved_data/56.", data_type, "_level.bs", bin_size, ".states_across_geneBody_up", up_bin_num, "_down", down_bin_num, ".", cell, ".qs")
      
      dat <- qread(file, nthreads = 6)
      
      return(dat)
    }
    get_gene_cluster_dat <- function(data_type, cell, method) {
      # test data
      if (F) {
        cell <- cell1
      }
      
      file <- paste0("data/saved_data/56.", data_type, "_level.bs", bin_size, ".", method, ".clusters_across_geneBody_up", up_bin_num, "_down", down_bin_num, ".", cell, ".qs")
      
      dat <- qread(file, nthreads = 6)
      
      return(dat)
    }
    
    # summary for dcsg
    dcsg_analysis <- function(cell1_dat, cell2_dat, dists) {
      # prepare hello info
      if (T) {
        start_mess <- paste0("|", paste(rep("-", 100), collapse = ""), "|\n")
        cat(start_mess)
        
        breaks <- round(seq(1, length(cell1_dat), length.out=100))
        breaks <- names(cell1_dat)[breaks]
      }
      
      # dcsg analysis
      if (T) {
        gene_dcsc <- lapply(names(cell1_dat), function(id) {
          # id <- names(cell1_dat)[1]
          # id <- "ENSG00000289544.1"
          
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
          
          # get cell specific state seq
          if (T) {
            dat1 <- cell1_dat[[id]]
            dat2 <- cell2_dat[[id]]
          }
          
          # summary the res
          if (T) {
            res <- data.frame(
              id = id, 
              location_id = 1:length(dat1), 
              location_p = (1:length(dat1)) / length(dat1), 
              cell1 = dat1, 
              cell2 = dat2
            )
            colnames(res) <- c("id", "location_id", "location_p", cell1, cell2)
          }
          
          # distance
          if (T) {
            res$distance <- sapply(1:nrow(res), function(x) {
              # x <- 1
              dists[paste0("S", res[, cell1][x]), paste0("S", res[, cell2][x])]
            })
          }
          
          return(res)
        })
        gene_dcsc <- do.call(rbind, gene_dcsc)
      }
      
      return(gene_dcsc)
    }
    
    # calculate the distance for dcsg cs cluster
    dcsg_cs_cluster_distance <- function(cs_dcsg, cluster_dcsg) {
      # get consistent gene id
      if (T) {
        ids <- unique(cluster_dcsg$id)
        cs_dcsg <- cs_dcsg[cs_dcsg$id %in% ids, ]
      }
      
      # calculate the distance
      if (T) {
        cs_dcsg_dist <- data.frame(
          id = paste0(cs_dcsg$id, cs_dcsg$dcs_loc_absolute), 
          dist = cs_dcsg$distance
        )
        cluster_dcsg_dist <- paste0(cluster_dcsg$id, cluster_dcsg$dcs_loc_absolute)
        
        cluster_dcsg$distance <- cs_dcsg_dist$dist[match(cluster_dcsg_dist, cs_dcsg_dist$id)]
      }
      
      return(cluster_dcsg)
    }
  }
  
  # get parameters
  if (T) {
    data_type <- str_split(id, "@", simplify = T)[, 1]
    method <- str_split(id, "@", simplify = T)[, 2]
  }
  
  # loop body
  if (T) {
    # ------------- CS ------------- #
    
    # load distance value
    if (T) {
      if (method == "chromIDEAS") {
        dists <- qread(file = "data/saved_data/55.cs_distance.qs", nthreads = 6)
      } else {
        file <- "data/saved_data/57.emission_cs_distance.qs"
        
        if (! file.exists(file)) {
          # read emission table
          if (T) {
            emission <- read.table("data/raw_data/1.emission/chromIDEAS.emission.txt", header = T, sep = "\t", fill = T, comment.char = "")
            rownames(emission) <- emission[, 1]
            emission <- emission[, -c(1:2)]
            emission <- scale(emission)
          }
          
          # calculate the distance (modify based on script28)
          if (T) {
            dists <- as.matrix(dist(emission, method = "euclidean"))
          }
          
          # normalization: max-min-normalization
          if (T) {
            values <- c(as.matrix(dists))
            max <- max(values)
            min <- min(values)
            
            dists <- (dists-min)/(max-min)
            dists <- as.matrix(dists)
            
            rm(values, max, min)
          }
          
          qsave(dists, file = file, nthreads = 6)
          rm(emission)
        }
        if (file.exists(file)) {
          dists <- qread(file, nthreads = 6)
        }
      }
    }
    
    # calculate differential CS state genes
    if (T) {
      file <- paste0("data/saved_data/57.", data_type, "_level_de_novo_DCSG_CS_(", method, "_dist).qs")
      
      if (! file.exists(file)) {
        # get gene CS data
        if (T) {
          cell1_dat <- get_gene_cs_dat(data_type, cell1)
          cell2_dat <- get_gene_cs_dat(data_type, cell2)
        }
        
        # calculate distance based on bin
        if (identical(names(cell1_dat), names(cell2_dat))) {
          # 3159848
          dcsg_detail_cs <- dcsg_analysis(cell1_dat, cell2_dat, dists)
        }
        
        # save the data
        qsave(dcsg_detail_cs, file = file, nthreads = 6)
      }
    }
    
    # ------------- cluster ------------- #
    
    # load distance value
    if (T) {
      if (method == "chromIDEAS") {
        dists <- qread(file = "data/saved_data/55.cluster_distance_merged_heatmap.qs", nthreads = 6)
      } else {
        file <- "data/saved_data/57.emission_cs_kmeans_cluster_distance.qs"
        
        if (! file.exists(file)) {
          # load state distance
          if (T) {
            dists_cs <- qread("data/saved_data/57.emission_cs_distance.qs")
          }
          
          # load cluster data
          if (T) {
            cluster_emission_kmeans <- read.table("results/1.tab/40.emission_kmeans_cluster.csv", header = T, sep = ",")
          }
          
          # calculate the distance for each cluster
          if (T) {
            all_clusters <- sort(unique(cluster_emission_kmeans$cluster5))
            dists <- lapply(all_clusters, function(c) {
              # c <- 1
              
              dis <- sapply(all_clusters, function(r_c) {
                # r_c <- all_clusters[1]
                r_states <- cluster_emission_kmeans$state[cluster_emission_kmeans$cluster5 == r_c]
                c_states <- cluster_emission_kmeans$state[cluster_emission_kmeans$cluster5 == c]
                
                mean(dists_cs[c_states, r_states], na.rm = T)
              })
              
              dis <- data.frame(target_c = c, 
                                compared_c = all_clusters, 
                                distance = dis)
              rownames(dis) <- 1:nrow(dis)
              
              return(dis)
            })
            dists <- do.call(rbind, dists)
          }
          
          # long2wide
          if (T) {
            head(dists)
            dists <- dcast(dists, target_c~compared_c, value.var = "distance")
            rownames(dists) <- paste0("C", dists[, 1])
            
            dists <- dists[, -1]
            colnames(dists) <- paste0("C", colnames(dists))
            
            dists <- as.matrix(dists)
          }
          
          qsave(dists, file = file, nthreads = 6)
        }
        if (file.exists(file)) {
          dists <- qread(file, nthreads = 6)
        }
      }
      
      rownames(dists) <- gsub("C", "S", rownames(dists))
      colnames(dists) <- gsub("C", "S", colnames(dists))
    }
    
    # calculate differential CS cluster genes
    if (T) {
      file <- paste0("data/saved_data/57.", data_type, "_level_de_novo_DCSG_CS_cluster_dist_(", method, ").qs")
      
      if (! file.exists(file)) {
        # get gene CS data
        if (T) {
          cell1_dat <- get_gene_cluster_dat(data_type, cell1, method)
          cell2_dat <- get_gene_cluster_dat(data_type, cell2, method)
        }
        
        # calculate distance based on bin
        if (identical(names(cell1_dat), names(cell2_dat))) {
          dcsg_detail_cs <- dcsg_analysis(cell1_dat, cell2_dat, dists)
        }
        
        # save the data
        qsave(dcsg_detail_cs, file = file, nthreads = 6)
      }
    }
  }
}
if (T) {
  foreach(id = paste0(rep(c("gene", "tx"), 2), "@", rep(c("chromIDEAS", "kmeans"), each=2)), .export =ls()) %dopar% loop(id)
}

stopCluster(cl)