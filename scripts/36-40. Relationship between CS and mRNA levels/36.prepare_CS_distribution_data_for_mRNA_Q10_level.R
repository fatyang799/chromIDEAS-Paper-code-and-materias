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
  # data prepare
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
      cell_mrna_non0[["subdat0"]] <- cell_mrna_0
      
      cell_mrna <- cell_mrna_non0
      
      rm(cell_mrna_0, cell_mrna_non0)
    }
    
    return(cell_mrna)
  }
  
  # from script 12
  get_cell_specific_tss_matrix_sum <- function(input, IDs, state_order) {
    # test data
    if (F) {
      input <- "data/saved_data/12.gene_level_profile_TSS_5UP_5DW.bs200.raw_mat.cd34.qs"
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
    
    # get total number for each state in each position
    if (T) {
      state_profile <- lapply(state_order, function(s){
        # s <- 0
        cell_dat <- sapply(dat[, -1], function(c_s) {
          sum(s==c_s)
        })
        return(cell_dat)
      })
      
      state_profile <- data.frame(do.call(rbind, state_profile))
      state_profile$State <- paste0("S", state_order)
    }
    
    return(state_profile)
  }
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
  
  # from script 12
  format_cell_specific_percentage_matrix <- function(dat, scale_type) {
    # test data
    if (F) {
      dat <- get_cell_specific_body_matrix_sum("data/saved_data/12.tx_level_profile_Body_5UP_5DW_10Body.bs200.raw_mat.cd34.qs", 
                                               mrna$cd34$subdat_hig, 
                                               sort(unique(state[, "cd34"])))
      scale_type <- "genomic"
    }
    
    if (scale_type == "genomic") {
      dat_col <- data.frame(apply(dat[, grep("^U|^B|^D|TSS|TES", colnames(dat))], 2, function(x) {
        x/sum(x)*100
      }))
      dat_col$State <- dat$State
      
      dat_col <- melt(dat_col, id.vars = "State", variable.name = "Loc", value.name = "Genomic_Percentage")
      dat_col$State <- factor(dat_col$State, levels = paste0("S", sort(as.numeric(gsub("S", "", unique(dat_col$State))))))
      
      if ("TES" %in% unique(dat_col$Loc)) {
        dat_col$Loc <- factor(dat_col$Loc, levels = c(
          paste0("U", up_bin_num:1), 
          "TSS", 
          paste0("B", 1:body_bin_num), 
          "TES", 
          paste0("D", 1:down_bin_num)
        ))
      }
      if (! "TES" %in% unique(dat_col$Loc)) {
        dat_col$Loc <- factor(dat_col$Loc, levels = c(
          paste0("U", up_bin_num:1), 
          "TSS", 
          paste0("D", 1:down_bin_num)
        ))
      }
      
      dat_out <- dat_col
    }
    if (scale_type == "state") {
      bin_num <- sapply(dat[, grep("^U|^B|^D|TSS|TES", colnames(dat))], sum)
      
      dat_row <- data.frame(t(
        apply(dat[, grep("^U|^B|^D|TSS|TES", colnames(dat))], 1, function(x) {
          x <- as.numeric(x)
          # norm the bin number
          (x/sum(x)*100)/(bin_num/bin_num["TSS"])
        })
      ))
      dat_row$State <- dat$State
      
      dat_row <- melt(dat_row, id.vars = "State", variable.name = "Loc", value.name = "State_Percentage")
      dat_row$State <- factor(dat_row$State, levels = paste0("S", sort(as.numeric(gsub("S", "", unique(dat_row$State))))))
      
      if ("TES" %in% unique(dat_row$Loc)) {
        dat_row$Loc <- factor(dat_row$Loc, levels = c(
          paste0("U", up_bin_num:1), 
          "TSS", 
          paste0("B", 1:body_bin_num), 
          "TES", 
          paste0("D", 1:down_bin_num)
        ))
      }
      if (! "TES" %in% unique(dat_row$Loc)) {
        dat_row$Loc <- factor(dat_row$Loc, levels = c(
          paste0("U", up_bin_num:1), 
          "TSS", 
          paste0("D", 1:down_bin_num)
        ))
      }
      
      dat_out <- dat_row
    }
    
    return(dat_out)
  }
}

# data prepare for each cell
if (T) {
  data_type <- "tx"
  state_order <- sort(unique(state[, 2]))
  
  for (type in c("TSS", "Body")) {
    cat(paste0(type, ":\n"))
    for (scale_type in c("genomic", "state")) {
      cat(paste0("\t", scale_type, ":\n"))
      for (cell in c(cell1, cell2)) {
        cat(paste0("\t\t", cell, ": \n"))
        
        # test data
        if (F) {
          cell <- cell1
          type <- "Body"
          scale_type <- "genomic"
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
                           paste0("data/saved_data/36.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                  "DW.bs", bin_size, ".raw_mRNA_classification_Q", rna_levels, "_mat.", cell, ".qs"), 
                           paste0("data/saved_data/36.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
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
        
        # format the profile
        if (T) {
          # get output filename
          if (T) {
            file <- ifelse(type == "TSS", 
                           paste0("data/saved_data/36.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                  "DW.bs", bin_size, ".format_", scale_type, 
                                  "_percentage_mRNA_classification_Q", rna_levels, "_mat.", cell, ".qs"), 
                           paste0("data/saved_data/36.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                  "DW_", body_bin_num, "Body.bs", bin_size, ".format_", scale_type, 
                                  "_percentage_mRNA_classification_Q", rna_levels, "_mat.", cell, ".qs"))
          }
          
          if (! file.exists(file)) {
            profile <- lapply(names(mat), function(x) {
              # x <- names(mat)[1]
              dat <- mat[[x]]
              ggdat <- format_cell_specific_percentage_matrix(dat, scale_type)
              ggdat$mRNA <- x
              
              return(ggdat)
            })
            profile <- data.frame(do.call(rbind, profile))
            profile$mRNA <- factor(profile$mRNA, levels = c("subdat0", paste0("Q", 1:rna_levels)))
            if (scale_type == "genomic") {
              profile <- profile[, c("State", "Loc", "Genomic_Percentage", "mRNA")]
            }
            if (scale_type == "state") {
              profile <- profile[, c("State", "Loc", "State_Percentage", "mRNA")]
            }
            
            qsave(profile, file, nthreads = 6)
          }
        }
      }
    }
  }
}

# data prepare for merged data
if (T) {
  data_type <- "tx"
  
  for (type in c("TSS", "Body")) {
    cat(paste0(type, ":\n"))
    for (scale_type in c("genomic", "state")) {
      cat(paste0("\t", scale_type, ":\n"))
      
      # test data
      if (F) {
        type <- "Body"
        scale_type <- "genomic"
      }
      
      # merge mRNA specific profile
      if (T) {
        # get input filename
        if (T) {
          file1 <- ifelse(type == "TSS", 
                          paste0("data/saved_data/36.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                 "DW.bs", bin_size, ".raw_mRNA_classification_Q", rna_levels, "_mat.", cell1, ".qs"), 
                          paste0("data/saved_data/36.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                 "DW_", body_bin_num, "Body.bs", bin_size, ".raw_mRNA_classification_Q", rna_levels, "_mat.", cell1, ".qs"))
          file2 <- ifelse(type == "TSS", 
                          paste0("data/saved_data/36.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                 "DW.bs", bin_size, ".raw_mRNA_classification_Q", rna_levels, "_mat.", cell2, ".qs"), 
                          paste0("data/saved_data/36.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
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
      }
      
      # format the profile
      if (T) {
        # get output filename
        if (T) {
          file <- ifelse(type == "TSS", 
                         paste0("data/saved_data/36.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                "DW.bs", bin_size, ".format_", scale_type, 
                                "_percentage_mRNA_classification_Q", rna_levels, "_mat.merged.qs"), 
                         paste0("data/saved_data/36.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                "DW_", body_bin_num, "Body.bs", bin_size, ".format_", scale_type, 
                                "_percentage_mRNA_classification_Q", rna_levels, "_mat.merged.qs"))
        }
        
        if (! file.exists(file)) {
          profile <- lapply(names(mat), function(x) {
            # x <- names(mat)[1]
            dat <- mat[[x]]
            ggdat <- format_cell_specific_percentage_matrix(dat, scale_type)
            ggdat$mRNA <- x
            
            return(ggdat)
          })
          profile <- data.frame(do.call(rbind, profile))
          profile$mRNA <- factor(profile$mRNA, levels = c("subdat0", paste0("Q", 1:rna_levels)))
          if (scale_type == "genomic") {
            profile <- profile[, c("State", "Loc", "Genomic_Percentage", "mRNA")]
          }
          if (scale_type == "state") {
            profile <- profile[, c("State", "Loc", "State_Percentage", "mRNA")]
          }
          
          qsave(profile, file, nthreads = 6)
        }
      }
    }
  }
}
