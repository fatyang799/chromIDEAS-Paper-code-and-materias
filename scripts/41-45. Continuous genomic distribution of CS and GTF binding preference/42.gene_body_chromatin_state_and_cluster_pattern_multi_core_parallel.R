library(ggplot2)
library(stringr)
library(reshape2)
library(qs)
library(doParallel)
library(parallel)
library(foreach)
cl <- makeCluster(20, outfile="log42.txt")
registerDoParallel(cl)

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
up_bin_num <- 5
down_bin_num <- 5
body_bin_num <- 10
rna_levels <- 10
n_ct <- 5

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
  }
  
  get_states_total_n <- function(raw_state_dat, all_states) {
    # get nonbody region
    if (T) {
      nonbody <- raw_state_dat[, grep("^U|TSS|TES|^D", colnames(raw_state_dat))]
      nonbody <- data.frame(sapply(nonbody, table))
      nonbody <- apply(nonbody, 1, sum)
    }
    
    # get body region
    if (T) {
      body <- raw_state_dat[, grep("^B[0-9]{1,2}_S[0-9]{1,2}", colnames(raw_state_dat))]
      body <- sapply(all_states, function(s) {
        # s <- 0
        subdat <- raw_state_dat[, grep(paste0("^B[0-9]{1,2}_S", s, "$"), colnames(raw_state_dat))]
        subdat <- sapply(subdat, sum)
        subdat <- sum(subdat)
        
        return(subdat)
      })
      names(body) <- all_states
    }
    
    # make sure the order identical
    if (T) {
      nonbody <- nonbody[as.character(all_states)]
      body <- body[as.character(all_states)]
      states_n <- nonbody+body
    }
    
    return(states_n)
  }
  
  select_most_prefer_states <- function(raw_state_dat, ratio) {
    # prepare hello info
    if (T) {
      start_mess <- paste0("|", paste(rep("-", 100), collapse = ""), "|\n")
      cat(start_mess)
      
      raw_state_dat$num <- 1:nrow(raw_state_dat)
      breaks <- round(seq(1, nrow(raw_state_dat), length.out=100))
    }
    
    # get proxy state
    if (T) {
      mat <- data.frame(t(
        apply(raw_state_dat, 1, function(x) {
          # x <- unlist(raw_state_dat[1, ])
          # print process
          if (as.numeric(x["num"]) %in% breaks) {
            if (which(as.numeric(x["num"]) == breaks) == 1) {
              cat("|*")
            }
            if (which(as.numeric(x["num"]) == breaks) == 100) {
              cat("*|\n")
            }
            if (! which(as.numeric(x["num"]) == breaks) %in% c(1, 100)) {
              cat("*")
            }
          }
          
          # get nonbody region
          if (T) {
            nonbody <- x[grep("^U|TSS|TES|^D", names(x))]
          }
          
          # get body region
          if (T) {
            body <- sapply(1:body_bin_num, function(b) {
              # b <- 1
              subdat <- x[grep(paste0("B", b, "_S"), names(x))]
              
              # format the data
              if (T) {
                subdat <- data.frame(id = names(subdat), 
                                     number = as.numeric(subdat))
                subdat$percentage <- subdat$number/sum(subdat$number)
                
                subdat$state <- str_split(subdat$id, "S", simplify = T)[, 2]
                
                # make sure the order is identical
                subdat <- subdat[match(names(ratio), subdat$state), ]
                
                # get prefer dat
                subdat$prefer <- subdat$percentage/ratio
              }
              
              # get most prefer state
              if (T) {
                most_level <- max(subdat$prefer)
                target_state <- subdat$state[subdat$prefer == most_level]
                target_state <- paste(target_state, collapse = "+")
              }
              
              return(target_state)
            })
            names(body) <- paste0("B", 1:body_bin_num)
          }
          
          # merge the results
          if (T) {
            res <- c(nonbody, body)
            res <- res[order]
          }
          
          return(res)
        })
      ))
    }
    
    # format the data
    if (T) {
      mat$gene_id <- raw_state_dat$gene_id
      
      mat <- mat[, c(which(grepl("gene_id", colnames(mat))), 
                     which(!grepl("gene_id", colnames(mat))))]
      
      mat <- data.frame(sapply(mat, function(x) {
        gsub(" ", "", x)
      }))
    }
    
    return(mat)
  }
  
  select_most_prefer_clusters <- function(raw_state_dat, cluster_dat, states_n) {
    # prepare hello info
    if (T) {
      start_mess <- paste0("|", paste(rep("-", 100), collapse = ""), "|\n")
      cat(start_mess)
      
      raw_state_dat$num <- 1:nrow(raw_state_dat)
      breaks <- round(seq(1, nrow(raw_state_dat), length.out=100))
    }
    
    # get the cluster ratio
    if (T) {
      cluster_dat$total_n <- states_n[match(cluster_dat$state, as.numeric(names(states_n)))]
      ratio <- tapply(cluster_dat$total_n, cluster_dat$cluster, sum)
      ratio <- ratio/sum(ratio)
    }
    
    # get proxy state
    if (T) {
      mat <- data.frame(t(
        apply(raw_state_dat, 1, function(x) {
          # x <- unlist(raw_state_dat[1, ])
          # print process
          if (as.numeric(x["num"]) %in% breaks) {
            if (which(as.numeric(x["num"]) == breaks) == 1) {
              cat("|*")
            }
            if (which(as.numeric(x["num"]) == breaks) == 100) {
              cat("*|\n")
            }
            if (! which(as.numeric(x["num"]) == breaks) %in% c(1, 100)) {
              cat("*")
            }
          }
          
          # get nonbody region
          if (T) {
            nonbody <- x[grep("^U|TSS|TES|^D", names(x))]
            nonbody <- cluster_dat$cluster[match(as.numeric(nonbody), cluster_dat$state)]
            names(nonbody) <- grep("^U|TSS|TES|^D", names(x), value = T)
          }
          
          # get body region
          if (T) {
            body <- sapply(1:body_bin_num, function(b) {
              # b <- 1
              subdat <- x[grep(paste0("B", b, "_S"), names(x))]
              
              # format the data
              if (T) {
                subdat <- data.frame(id = names(subdat), 
                                     number = as.numeric(subdat))
                subdat$state <- str_split(subdat$id, "S", simplify = T)[, 2]
                subdat$cluster <- cluster_dat$cluster[match(as.numeric(subdat$state), cluster_dat$state)]
                
                subdat <- tapply(subdat$number, subdat$cluster, sum)
                subdat <- data.frame(cluster = names(subdat), 
                                     number = as.numeric(subdat))
                
                subdat$percentage <- subdat$number/sum(subdat$number)
                
                # make sure the order is identical
                subdat <- subdat[match(names(ratio), subdat$cluster), ]
                
                # get prefer dat
                subdat$prefer <- subdat$percentage/ratio
              }
              
              # get most prefer cluster
              if (T) {
                most_level <- max(subdat$prefer)
                target_cluster <- subdat$cluster[subdat$prefer == most_level]
                target_cluster <- paste(target_cluster, collapse = "+")
              }
              
              return(target_cluster)
            })
            names(body) <- paste0("B", 1:body_bin_num)
          }
          
          # merge the results
          if (T) {
            res <- c(nonbody, body)
            res <- res[order]
          }
          
          return(res)
        })
      ))
      
      mat$gene_id <- raw_state_dat$gene_id
      
      mat <- mat[, c(which(grepl("gene_id", colnames(mat))), 
                     which(!grepl("gene_id", colnames(mat))))]
    }
    
    # format the data
    if (T) {
      mat <- data.frame(sapply(mat, function(x) {
        gsub(" ", "", x)
      }))
    }
    
    return(mat)
  }
}

# get representative state
replace_Cluster <- function(cluster) {
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
  }
  
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
    }
    
    get_states_total_n <- function(raw_state_dat, all_states) {
      # get nonbody region
      if (T) {
        nonbody <- raw_state_dat[, grep("^U|TSS|TES|^D", colnames(raw_state_dat))]
        nonbody <- data.frame(sapply(nonbody, table))
        nonbody <- apply(nonbody, 1, sum)
      }
      
      # get body region
      if (T) {
        body <- raw_state_dat[, grep("^B[0-9]{1,2}_S[0-9]{1,2}", colnames(raw_state_dat))]
        body <- sapply(all_states, function(s) {
          # s <- 0
          subdat <- raw_state_dat[, grep(paste0("^B[0-9]{1,2}_S", s, "$"), colnames(raw_state_dat))]
          subdat <- sapply(subdat, sum)
          subdat <- sum(subdat)
          
          return(subdat)
        })
        names(body) <- all_states
      }
      
      # make sure the order identical
      if (T) {
        nonbody <- nonbody[as.character(all_states)]
        body <- body[as.character(all_states)]
        states_n <- nonbody+body
      }
      
      return(states_n)
    }
    
    select_most_prefer_states <- function(raw_state_dat, ratio) {
      # prepare hello info
      if (T) {
        start_mess <- paste0("|", paste(rep("-", 100), collapse = ""), "|\n")
        cat(start_mess)
        
        raw_state_dat$num <- 1:nrow(raw_state_dat)
        breaks <- round(seq(1, nrow(raw_state_dat), length.out=100))
      }
      
      # get proxy state
      if (T) {
        mat <- data.frame(t(
          apply(raw_state_dat, 1, function(x) {
            # x <- unlist(raw_state_dat[1, ])
            # print process
            if (as.numeric(x["num"]) %in% breaks) {
              if (which(as.numeric(x["num"]) == breaks) == 1) {
                cat("|*")
              }
              if (which(as.numeric(x["num"]) == breaks) == 100) {
                cat("*|\n")
              }
              if (! which(as.numeric(x["num"]) == breaks) %in% c(1, 100)) {
                cat("*")
              }
            }
            
            # get nonbody region
            if (T) {
              nonbody <- x[grep("^U|TSS|TES|^D", names(x))]
            }
            
            # get body region
            if (T) {
              body <- sapply(1:body_bin_num, function(b) {
                # b <- 1
                subdat <- x[grep(paste0("B", b, "_S"), names(x))]
                
                # format the data
                if (T) {
                  subdat <- data.frame(id = names(subdat), 
                                       number = as.numeric(subdat))
                  subdat$percentage <- subdat$number/sum(subdat$number)
                  
                  subdat$state <- str_split(subdat$id, "S", simplify = T)[, 2]
                  
                  # make sure the order is identical
                  subdat <- subdat[match(names(ratio), subdat$state), ]
                  
                  # get prefer dat
                  subdat$prefer <- subdat$percentage/ratio
                }
                
                # get most prefer state
                if (T) {
                  most_level <- max(subdat$prefer)
                  target_state <- subdat$state[subdat$prefer == most_level]
                  target_state <- paste(target_state, collapse = "+")
                }
                
                return(target_state)
              })
              names(body) <- paste0("B", 1:body_bin_num)
            }
            
            # merge the results
            if (T) {
              res <- c(nonbody, body)
              res <- res[order]
            }
            
            return(res)
          })
        ))
      }
      
      # format the data
      if (T) {
        mat$gene_id <- raw_state_dat$gene_id
        
        mat <- mat[, c(which(grepl("gene_id", colnames(mat))), 
                       which(!grepl("gene_id", colnames(mat))))]
        
        mat <- data.frame(sapply(mat, function(x) {
          gsub(" ", "", x)
        }))
      }
      
      return(mat)
    }
    
    select_most_prefer_clusters <- function(raw_state_dat, cluster_dat, states_n) {
      # prepare hello info
      if (T) {
        start_mess <- paste0("|", paste(rep("-", 100), collapse = ""), "|\n")
        cat(start_mess)
        
        raw_state_dat$num <- 1:nrow(raw_state_dat)
        breaks <- round(seq(1, nrow(raw_state_dat), length.out=100))
      }
      
      # get the cluster ratio
      if (T) {
        cluster_dat$total_n <- states_n[match(cluster_dat$state, as.numeric(names(states_n)))]
        ratio <- tapply(cluster_dat$total_n, cluster_dat$cluster, sum)
        ratio <- ratio/sum(ratio)
      }
      
      # get proxy state
      if (T) {
        mat <- data.frame(t(
          apply(raw_state_dat, 1, function(x) {
            # x <- unlist(raw_state_dat[1, ])
            # print process
            if (as.numeric(x["num"]) %in% breaks) {
              if (which(as.numeric(x["num"]) == breaks) == 1) {
                cat("|*")
              }
              if (which(as.numeric(x["num"]) == breaks) == 100) {
                cat("*|\n")
              }
              if (! which(as.numeric(x["num"]) == breaks) %in% c(1, 100)) {
                cat("*")
              }
            }
            
            # get nonbody region
            if (T) {
              nonbody <- x[grep("^U|TSS|TES|^D", names(x))]
              nonbody <- cluster_dat$cluster[match(as.numeric(nonbody), cluster_dat$state)]
              names(nonbody) <- grep("^U|TSS|TES|^D", names(x), value = T)
            }
            
            # get body region
            if (T) {
              body <- sapply(1:body_bin_num, function(b) {
                # b <- 1
                subdat <- x[grep(paste0("B", b, "_S"), names(x))]
                
                # format the data
                if (T) {
                  subdat <- data.frame(id = names(subdat), 
                                       number = as.numeric(subdat))
                  subdat$state <- str_split(subdat$id, "S", simplify = T)[, 2]
                  subdat$cluster <- cluster_dat$cluster[match(as.numeric(subdat$state), cluster_dat$state)]
                  
                  subdat <- tapply(subdat$number, subdat$cluster, sum)
                  subdat <- data.frame(cluster = names(subdat), 
                                       number = as.numeric(subdat))
                  
                  subdat$percentage <- subdat$number/sum(subdat$number)
                  
                  # make sure the order is identical
                  subdat <- subdat[match(names(ratio), subdat$cluster), ]
                  
                  # get prefer dat
                  subdat$prefer <- subdat$percentage/ratio
                }
                
                # get most prefer cluster
                if (T) {
                  most_level <- max(subdat$prefer)
                  target_cluster <- subdat$cluster[subdat$prefer == most_level]
                  target_cluster <- paste(target_cluster, collapse = "+")
                }
                
                return(target_cluster)
              })
              names(body) <- paste0("B", 1:body_bin_num)
            }
            
            # merge the results
            if (T) {
              res <- c(nonbody, body)
              res <- res[order]
            }
            
            return(res)
          })
        ))
        
        mat$gene_id <- raw_state_dat$gene_id
        
        mat <- mat[, c(which(grepl("gene_id", colnames(mat))), 
                       which(!grepl("gene_id", colnames(mat))))]
      }
      
      # format the data
      if (T) {
        mat <- data.frame(sapply(mat, function(x) {
          gsub(" ", "", x)
        }))
      }
      
      return(mat)
    }
  }
  
  # loop body
  if (T) {
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
    
    # get best state cluster across gene body: based on the preference of cluster
    if (T) {
      file <- paste0(dir, "/42.", data_type, "_Body_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, 
                     "Body.bs", bin_size, ".", type, ".clusters(", cluster, ").", cell, ".qs")
      
      if (! file.exists(file)) {
        cat(paste0("\t\tNow get best state cluster across gene body: based on the preference of cluster\n"))
        
        mat <- select_most_prefer_clusters(raw_state_dat, cluster_dat, states_n)
        qsave(mat, file = file, nthreads = 6)
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
        cat(paste0("\t\t", cell, ": \n"))
        
        # test data
        if (F) {
          data_type <- "tx"
          type <- "merged"
          cell <- cell1
        }
        
        # mkdir
        if (T) {
          dir <- paste0("data/saved_data/42.gene_state_pattern/all_txs")
          if (! dir.exists(dir)) {
            dir.create(dir, showWarnings = F, recursive = T)
          }
        }
        
        # load the chromatin state across genes: raw_state_dat
        if (T) {
          file <- paste0("data/saved_data/12.", data_type, "_level_profile_Body_", up_bin_num, "UP_", 
                         down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".raw_mat.", cell, ".qs")
          raw_state_dat <- qread(file)
        }
        
        # common value
        if (T) {
          all_states <- sort(unique(as.numeric(raw_state_dat$TSS)))
        }
        
        # get state ratio across gene body
        if (T) {
          file <- paste0(dir, "/42.", data_type, "_level_profile_Body_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, 
                         "Body.bs", bin_size, ".state_genomic_percentage_ratio.", cell, ".qs")
          if (! file.exists(file)) {
            gene_base_states <- qread(paste0("data/saved_data/12.gene_level_profile_Body_", up_bin_num, "UP_", 
                                             down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".raw_mat.", cell, ".qs"))
            states_n <- get_states_total_n(gene_base_states, all_states)
            rm(gene_base_states)
            
            # get the ratio
            ratio <- states_n/sum(states_n)
            
            qsave(ratio, file = file, nthreads = 6)
          }
          if (file.exists(file)) {
            ratio <- qread(file, nthreads = 6)
          }
        }
        
        # get best chromatin states across gene body: based on the preference of states
        if (T) {
          file <- paste0(dir, "/42.", data_type, "_level_profile_Body_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, 
                         "Body.bs", bin_size, ".state.", cell, ".qs")
          
          if (! file.exists(file)) {
            cat(paste0("\t\tNow get best chromatin states across gene body: based on the preference of states\n"))
            mat <- select_most_prefer_states(raw_state_dat, ratio)
            
            qsave(mat, file, nthreads = 6)
            rm(mat)
          }
        }
        
        foreach(cluster = c("chromIDEAS", "kmeans", 1:n_ct), .export =ls()) %dopar% replace_Cluster(cluster)
      }
    }
  }
}

stopCluster(cl)