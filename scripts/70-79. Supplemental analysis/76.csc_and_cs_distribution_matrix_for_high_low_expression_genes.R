# copy of script 36 

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
rna_levels <- 2
state_file <- "data/raw_data/2.states/chromIDEAS.state"

# read the states data
if (T) {
  state <- data.table::fread(state_file, sep = " ", header = T, data.table = F)
  state <- state[, c(1, 5, 6)]
  colnames(state)[1] <- "ID"
  head(state)
}

# read mRNA data
if (T) {
  thp1_mrna <- "data/saved_data/7.4_type_tx_expression_level_classification_thp1.qs"
  thp1_mrna <- qread(thp1_mrna, nthreads = 6)
  
  cd34_mrna <- "data/saved_data/7.4_type_tx_expression_level_classification_cd34.qs"
  cd34_mrna <- qread(cd34_mrna, nthreads = 6)
  
  mrna <- list(thp1_mrna, cd34_mrna)
  names(mrna) <- c("thp1", "cd34")
  rm(thp1_mrna, cd34_mrna)
}

# define functions
if (T) {
  # data prepare (modify: q12)
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
      cell_mrna_non0[["Q1"]] <- c(cell_mrna_0, cell_mrna_non0[["Q1"]])
      
      cell_mrna <- cell_mrna_non0
      
      rm(cell_mrna_0, cell_mrna_non0)
    }
    
    return(cell_mrna)
  }
  
  # from script 12
  get_cell_specific_body_matrix_sum <- function(input, IDs, state_order) {
    # test data
    if (F) {
      input <- "data/saved_data/12.gene_level_profile_Body_5UP_5DW_10Body.bs200.raw_mat.cd34.qs"
      mat <- qread(input, nthreads = 6)
      IDs <- mat$gene_id
      state_order <- sort(as.numeric(unique(mat$TSS)))
    }
    
    # read and get specific data
    if (T) {
      mat <- qread(input, nthreads = 6)
      dat <- mat[mat$gene_id %in% IDs, ]
      
      mess <- sum(! IDs %in% mat$gene_id)
      if (mess > 0) {
        mess <- paste0("There are ", mess, "/", length(IDs), " (", round(mess/length(IDs)*100, 2), "%) genes in IDs not exist in matrix\n")
        cat(mess)
      }
    }
    
    # get total number for each state: nonbody
    if (T) {
      colnames(dat)
      
      nonbody <- dat[, grep("^U|^D|TSS|TES", colnames(dat))]
      
      nonbody <- lapply(state_order, function(s){
        # s <- 0
        nonbody_dat <- sapply(nonbody, function(c_s) {
          sum(s==c_s)
        })
        return(nonbody_dat)
      })
      
      nonbody <- data.frame(do.call(rbind, nonbody))
      nonbody$State <- paste0("S", state_order)
    }
    
    # get total number for each state: body
    if (T) {
      body <- dat[, grep("^B", colnames(dat))]
      
      body <- melt(body, id.vars = NULL, value.name = "value", variable.name = "ID")
      body$Loc <- str_split(body$ID, "_", simplify = T)[, 1]
      body$State <- str_split(body$ID, "_", simplify = T)[, 2]
      body <- dcast(body, State~Loc, value.var = "value", fun.aggregate = sum)
      body <- body[, c(paste0("B", 1:body_bin_num), "State")]
    }
    
    # merge the info 
    if (T) {
      # make sure the order of body and nonbody are identical
      if (T) {
        body <- body[match(paste0("S", state_order), body$State), ]
        nonbody <- nonbody[match(paste0("S", state_order), nonbody$State), ]
      }
      
      # merge the state number data
      if (T) {
        order <- c(
          paste0("U", up_bin_num:1), 
          "TSS", 
          paste0("B", 1:body_bin_num), 
          "TES", 
          paste0("D", 1:down_bin_num)
        )
        
        state_profile <- merge(body, nonbody, by="State")
        state_profile <- state_profile[, c(order, "State")]
        state_profile <- state_profile[match(paste0("S", state_order), state_profile$State), ]
        rownames(state_profile) <- 1:nrow(state_profile)
      }
    }
    
    return(state_profile)
  }
  
  # csc distribution matrix
  replace_cs_with_csc <- function(dat, cluster) {
    # test data
    if (F) {
      dat <- mat$Q1
    }
    
    dat$cluster <- cluster$seurat_clusters[match(dat$State, cluster$orig.ident)]
    dat <- lapply(split(dat[, c(1:(up_bin_num+1+body_bin_num+1+down_bin_num))], dat$cluster), function(subdat) {
      sapply(subdat, sum)
    })
    dat <- data.frame(do.call(rbind, dat))
    dat$cluster <- rownames(dat)
    
    return(dat)
  }
}

# data prepare for each cell
if (T) {
  data_type <- "tx"
  state_order <- sort(unique(state[, 2]))
  
  for (type in c("Body")) {
    cat(paste0(type, ":\n"))
    for (cell in c(cell1, cell2)) {
      cat(paste0("\t", cell, ": \n"))
      
      # test data
      if (F) {
        cell <- cell1
        type <- "Body"
      }
      
      # get input
      if (T) {
        input <- ifelse(type == "TSS", 
                        paste0("data/saved_data/12.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", 
                               down_bin_num, "DW.bs", bin_size, ".raw_mat.", cell, ".qs"), 
                        paste0("data/saved_data/12.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", 
                               down_bin_num, "DW_", body_bin_num, type, ".bs", bin_size, ".raw_mat.", cell, ".qs"))
      }
      
      # get IDs
      if (T) {
        cell_mrna <- mrna[[cell]][[1]]
        ids <- qread(input, nthreads = 6)
        cell_mrna <- cell_mrna[names(cell_mrna) %in% ids$gene_id]
        IDs <- quantitle_mRNA(cell_mrna)
        
        rm(cell_mrna, ids)
      }
      
      # get mRNA specific profile
      if (T) {
        # get function
        if (T) {
          fun <- ifelse(type == "TSS", 
                        get_cell_specific_tss_matrix_sum, 
                        get_cell_specific_body_matrix_sum)
        }
        
        # get output filename
        if (T) {
          file <- ifelse(type == "TSS", 
                         paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                "DW.bs", bin_size, ".raw_mRNA_classification_Q", rna_levels, "_mat.", cell, ".qs"), 
                         paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                "DW_", body_bin_num, "Body.bs", bin_size, ".raw_mRNA_classification_Q", rna_levels, "_mat.", cell, ".qs"))
        }
        
        # calculation
        if (T) {
          if (file.exists(file)) {
            mat <- qread(file, nthreads = 6)
          }
          if (! file.exists(file)) {
            mat <- lapply(IDs, function(x) {
              fun(input, x, state_order)
            })
            qsave(mat, file, nthreads = 6)
          }
        }
        
        rm(IDs)
      }
      
      # get cluster info
      if (T) {
        cluster_file <- ifelse(cell == "merged", 
                               "data/saved_data/29.common_cluster_name_for_merged_2cells_CS_clusters.qs", 
                               "data/saved_data/25.common_cluster_name_for_2_cell_CS_clusters.qs")
        cluster <- qread(cluster_file, nthreads = 6)
        cluster <- cluster[, c("orig.ident", "seurat_clusters", "cell")]
        cluster <- cluster[cluster$cell == cell, ]
        
        cluster$orig.ident <- ifelse(grepl("S", cluster$orig.ident, ignore.case = T), 
                                     cluster$orig.ident, 
                                     paste0("S", cluster$orig.ident))
      }
      
      # format the profile
      if (T) {
        # get output filename
        if (T) {
          file <- ifelse(type == "TSS", 
                         paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                "DW.bs", bin_size, ".format_CSC_distribution_matrix_mRNA_Q", 
                                rna_levels, ".", cell, ".qs"), 
                         paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                "DW_", body_bin_num, "Body.bs", bin_size, ".format_CSC_distribution_matrix_mRNA_Q", 
                                rna_levels, ".", cell, ".qs"))
        }
        
        if (! file.exists(file)) {
          profile <- lapply(names(mat), function(x) {
            # x <- names(mat)[1]
            dat <- mat[[x]]
            ggdat <- replace_cs_with_csc(dat, cluster)
            ggdat$mRNA <- x
            
            return(ggdat)
          })
          
          profile <- data.frame(do.call(rbind, profile))
          profile$mRNA <- factor(profile$mRNA, levels = c(paste0("Q", 1:rna_levels)))
          
          qsave(profile, file, nthreads = 6)
        }
      }
    }
  }
}

# data prepare for merged data
if (T) {
  data_type <- "tx"
  
  for (type in c("Body")) {
    cat(paste0(type, ":\n"))
    
    # test data
    if (F) {
      type <- "Body"
    }
    
    # merge mRNA specific profile
    if (T) {
      # get output filename
      if (T) {
        file <- ifelse(type == "TSS", 
                       paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                              "DW.bs", bin_size, ".raw_mRNA_classification_Q", rna_levels, "_mat.merged.qs"), 
                       paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                              "DW_", body_bin_num, "Body.bs", bin_size, ".raw_mRNA_classification_Q", rna_levels, "_mat.merged.qs"))
      }
      
      # calculation
      if (T) {
        if (file.exists(file)) {
          mat <- qread(file, nthreads = 6)
        }
        if (! file.exists(file)) {
          # get input filename
          if (T) {
            file1 <- ifelse(type == "TSS", 
                            paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                   "DW.bs", bin_size, ".raw_mRNA_classification_Q", rna_levels, "_mat.", cell1, ".qs"), 
                            paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                   "DW_", body_bin_num, "Body.bs", bin_size, ".raw_mRNA_classification_Q", rna_levels, "_mat.", cell1, ".qs"))
            file2 <- ifelse(type == "TSS", 
                            paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                   "DW.bs", bin_size, ".raw_mRNA_classification_Q", rna_levels, "_mat.", cell2, ".qs"), 
                            paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                   "DW_", body_bin_num, "Body.bs", bin_size, ".raw_mRNA_classification_Q", rna_levels, "_mat.", cell2, ".qs"))
          }
          
          # load cell data
          if (T) {
            dat1 <- qread(file1, nthreads = 6)
            dat2 <- qread(file2, nthreads = 6)
          }
          
          # merge the data
          if (identical(names(dat1), names(dat2))) {
            mat <- lapply(names(dat1), function(rna) {
              # rna <- names(dat1)[1]
              
              subdat1 <- dat1[[rna]]
              subdat2 <- dat2[[rna]]
              
              if (identical(colnames(subdat1), colnames(subdat2)) & identical(subdat1$State, subdat2$State)) {
                res <- subdat1[, colnames(subdat1) != "State"] + subdat2[, colnames(subdat2) != "State"]
                res$State <- subdat1$State
                
                return(res)
              } else {
                stop("Error1")
              }
            })
            names(mat) <- names(dat1)
          }
          
          rm(file1, file2, dat1, dat2)
          
          qsave(mat, file, nthreads = 6)
        }
      }
    }
    
    # format the profile
    if (T) {
      # get output filename
      if (T) {
        file <- ifelse(type == "TSS", 
                       paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                              "DW.bs", bin_size, ".format_CSC_distribution_matrix_mRNA_Q", 
                              rna_levels, ".merged.qs"), 
                       paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                              "DW_", body_bin_num, "Body.bs", bin_size, ".format_CSC_distribution_matrix_mRNA_Q", 
                              rna_levels, ".merged.qs"))
      }
      
      if (! file.exists(file)) {
        profile <- lapply(names(mat), function(x) {
          # x <- names(mat)[1]
          dat <- mat[[x]]
          ggdat <- replace_cs_with_csc(dat, cluster)
          ggdat$mRNA <- x
          
          return(ggdat)
        })
        
        profile <- data.frame(do.call(rbind, profile))
        profile$mRNA <- factor(profile$mRNA, levels = c(paste0("Q", 1:rna_levels)))
        
        qsave(profile, file, nthreads = 6)
      }
    }
  }
}
