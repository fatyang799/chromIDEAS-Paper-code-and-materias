library(stringr)
library(qs)
library(doParallel)
library(parallel)
library(foreach)
cl <- makeCluster(20, outfile="log58.txt")
registerDoParallel(cl)


cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
up_bin_num <- 3
down_bin_num <- 3

# define the funtion
if (T) {
  gene_cluster_conversion <- function(dat, tss_id) {
    # prepare hello info
    if (T) {
      start_mess <- paste0("|", paste(rep("-", 100), collapse = ""), "|\n")
      cat(start_mess)
      
      IDs <- unique(dat$id)
      breaks <- round(seq(1, length(IDs), length.out=100))
      breaks <- IDs[breaks]
    }
    
    # format the data
    if (T) {
      dat$cs_trans <- paste0(dat[, cell2], "to", dat[, cell1])
    }
    
    # all types
    if (T) {
      all_conversion <- sort(unique(dat$cs_trans))
    }
    
    # summary the cs conversion situation
    if (T) {
      cs_trans <- lapply(IDs, function(gene_id) {
        # gene_id <- IDs[1]
        
        # get subdat
        if (T) {
          subdat <- dat[dat$id == gene_id, ]
        }
        
        # print process
        if (gene_id %in% breaks) {
          if (which(gene_id == breaks) == 1) {
            cat("|*")
          }
          if (which(gene_id == breaks) == 100) {
            cat("*|\n")
          }
          if (! which(gene_id == breaks) %in% c(1, 100)) {
            cat("*")
          }
        }
        
        # stat distance
        if (T) {
          convers_dist <- tapply(subdat$distance, subdat$cs_trans, sum)
          convers_dist <- ifelse(all_conversion %in% names(convers_dist), convers_dist[all_conversion], 0)
          names(convers_dist) <- all_conversion
        }
        
        # stat loction
        if (T) {
          tss <- unique(subdat$start)
          tes <- unique(subdat$end)
          strand <- ifelse(tss<tes, "+", "-")
          tss_up <- ifelse(strand == "+", tss - up_bin_num,
                           ifelse(strand == "-", tss + up_bin_num, NA))
          tes_down <- ifelse(strand == "+", tes + down_bin_num,
                             ifelse(strand == "-", tes - down_bin_num, NA))
          
          gene_bin_id <- (tss_up:tes_down)
          dcsg_bin_id <- gene_bin_id[subdat$dcs_loc_absolute]
          
          loc_type <- ifelse(dcsg_bin_id %in% tss_id[[gene_id]], "TSS", "Body")
          convers_tss <- tapply(loc_type, subdat$cs_trans, function(x) {
            sum(x == "TSS")
          })
          convers_tss <- ifelse(all_conversion %in% names(convers_tss), convers_tss[all_conversion], 0)
          names(convers_tss) <- all_conversion
        }
        
        # summary the res
        if (T) {
          max_dist <- max(convers_dist)
          max_dist_type <- names(convers_dist)[convers_dist == max_dist]
          max_dist_N <- length(max_dist_type)
          max_dist_conversion = paste(max_dist_type, collapse = "@")
          
          max_tss <- max(convers_tss)
          max_tss_type <- names(convers_tss)[convers_tss == max_tss]
          max_tss_N <- length(max_tss_type)
          max_tss_conversion = paste(max_tss_type, collapse = "@")
          
          stat <- data.frame(
            id = gene_id, 
            max_dist = max_dist, 
            max_dist_N = max_dist_N, 
            max_dist_conversion = max_dist_conversion, 
            max_tss = max_tss, 
            max_tss_N = max_tss_N, 
            max_tss_conversion = max_tss_conversion)
        }
        
        return(stat)
      })
    }
    
    # merge and format the data
    if (T) {
      cs_trans <- do.call(rbind, cs_trans)
    }
    
    return(cs_trans)
  }
}

# data prepare
loop <- function(method) {
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
  
  # define the funtion
  if (T) {
    gene_cluster_conversion <- function(dat, tss_id) {
      # prepare hello info
      if (T) {
        start_mess <- paste0("|", paste(rep("-", 100), collapse = ""), "|\n")
        cat(start_mess)
        
        IDs <- unique(dat$id)
        breaks <- round(seq(1, length(IDs), length.out=100))
        breaks <- IDs[breaks]
      }
      
      # format the data
      if (T) {
        dat$cs_trans <- paste0(dat[, cell2], "to", dat[, cell1])
      }
      
      # all types
      if (T) {
        all_conversion <- sort(unique(dat$cs_trans))
      }
      
      # summary the cs conversion situation
      if (T) {
        cs_trans <- lapply(IDs, function(gene_id) {
          # gene_id <- IDs[1]
          
          # get subdat
          if (T) {
            subdat <- dat[dat$id == gene_id, ]
          }
          
          # print process
          if (gene_id %in% breaks) {
            if (which(gene_id == breaks) == 1) {
              cat("|*")
            }
            if (which(gene_id == breaks) == 100) {
              cat("*|\n")
            }
            if (! which(gene_id == breaks) %in% c(1, 100)) {
              cat("*")
            }
          }
          
          # stat distance
          if (T) {
            convers_dist <- tapply(subdat$distance, subdat$cs_trans, sum)
            convers_dist <- ifelse(all_conversion %in% names(convers_dist), convers_dist[all_conversion], 0)
            names(convers_dist) <- all_conversion
          }
          
          # stat loction
          if (T) {
            tss <- unique(subdat$start)
            tes <- unique(subdat$end)
            strand <- ifelse(tss<tes, "+", "-")
            tss_up <- ifelse(strand == "+", tss - up_bin_num,
                             ifelse(strand == "-", tss + up_bin_num, NA))
            tes_down <- ifelse(strand == "+", tes + down_bin_num,
                               ifelse(strand == "-", tes - down_bin_num, NA))
            
            gene_bin_id <- (tss_up:tes_down)
            dcsg_bin_id <- gene_bin_id[subdat$dcs_loc_absolute]
            
            loc_type <- ifelse(dcsg_bin_id %in% tss_id[[gene_id]], "TSS", "Body")
            convers_tss <- tapply(loc_type, subdat$cs_trans, function(x) {
              sum(x == "TSS")
            })
            convers_tss <- ifelse(all_conversion %in% names(convers_tss), convers_tss[all_conversion], 0)
            names(convers_tss) <- all_conversion
          }
          
          # summary the res
          if (T) {
            max_dist <- max(convers_dist)
            max_dist_type <- names(convers_dist)[convers_dist == max_dist]
            max_dist_N <- length(max_dist_type)
            max_dist_conversion = paste(max_dist_type, collapse = "@")
            
            max_tss <- max(convers_tss)
            max_tss_type <- names(convers_tss)[convers_tss == max_tss]
            max_tss_N <- length(max_tss_type)
            max_tss_conversion = paste(max_tss_type, collapse = "@")
            
            stat <- data.frame(
              id = gene_id, 
              max_dist = max_dist, 
              max_dist_N = max_dist_N, 
              max_dist_conversion = max_dist_conversion, 
              max_tss = max_tss, 
              max_tss_N = max_tss_N, 
              max_tss_conversion = max_tss_conversion)
          }
          
          return(stat)
        })
      }
      
      # merge and format the data
      if (T) {
        cs_trans <- do.call(rbind, cs_trans)
      }
      
      return(cs_trans)
    }
  }
  
  # loop body
  if (T) {
    # get DCSG representative state conversion
    file <- paste0("data/saved_data/58.", data_type, "_level_DCSG_representative_CS_Conversion(", method, ").qs")
    
    if (! file.exists(file)) {
      
      # get TSS bin IDs
      if (T) {
        tss_id <- qread("data/saved_data/11.tx_level_TSS_region_bins.qs", nthreads = 6)
        gtf <- qread("data/saved_data/4.tx_geneid_thp1_cd34.qs", nthreads = 6)
        
        tss_id$gene_id <- gtf$gene_id[match(tss_id$tx_id, gtf$tx_id)]
        
        tss_id <- split(tss_id$Bin_ID, tss_id$gene_id)
        rm(gtf)
      }
      
      # get Gene bin IDs
      if (T) {
        genebody_stat <- qread(paste0("data/saved_data/11.", data_type, "_level_Body_region_length_ge3_bins.qs"), nthreads = 6)
      }
      
      # get distance based on cs cluster
      if (T) {
        dat <- qread(paste0("data/saved_data/57.", data_type, "_level_de_novo_DCSG_CS_cluster_CSdist_(", method, ").qs"), nthreads = 6)
      }
      
      # merge the data
      if (T) {
        dat$start <- genebody_stat$tss_BinID[match(dat$id, genebody_stat$tx_id)]
        dat$end <- genebody_stat$tes_BinID[match(dat$id, genebody_stat$tx_id)]
        
        rm(genebody_stat)
      }
      
      # calculate the representative conversion
      if (T) {
        cs_conversion <- gene_cluster_conversion(dat, tss_id)
      }
      
      qsave(cs_conversion, file = file, nthreads = 6)
    }
  }
}
if (T) {
  for (data_type in c("gene")) {
    cat(paste0(data_type, ": \n"))
    
    # test data
    if (F) {
      data_type <- "gene"
    }
    
    foreach(method = c("chromIDEAS", "kmeans"), .export =ls()) %dopar% loop(method)
  }
}

stopCluster(cl)
