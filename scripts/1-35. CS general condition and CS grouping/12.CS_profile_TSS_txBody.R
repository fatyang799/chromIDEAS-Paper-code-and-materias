# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(stringr)
  library(reshape2)
  library(qs)
}

bin_size <- 200
up_bin_num <- 5
down_bin_num <- 5
body_bin_num <- 10
state_file <- "data/raw_data/2.states/chromIDEAS.state"

# read the states data
if (T) {
  state <- data.table::fread(state_file, sep = " ", header = T, data.table = F)
  state <- state[, c(1, 5, 6)]
  colnames(state)[1] <- "ID"
  head(state)
}

# define the funtion
if (T) {
  # data prepare
  get_cell_specific_tss_matrix <- function(cell, state, mat_file) {
    # test data
    if (F) {
      cell <- "thp1"
      mat_file <- paste0("data/saved_data/11.gene_level_profile_TSS_", up_bin_num, "UP_", down_bin_num, "DW.bs", bin_size, ".bin.qs")
    }
    
    # load the tss profile bins
    if (T) {
      mat <- qread(file = mat_file, nthreads = 6)
    }
    
    # get cell specific tss chromatin state matrix
    if (T) {
      dat <- data.frame(
        t(apply(mat, 1, function(x) {
          # x <- unlist(mat[1, ])
          geneid <- x["gene_id"]
          
          x <- x[grep("^U|TSS|^D", names(x))]
          cell_specific_state <- state[as.numeric(x), cell]
          res <- c(geneid, cell_specific_state)
          
          return(res)
        }))
      )
      colnames(dat)[2:ncol(dat)] <- grep("^U|TSS|^D", colnames(mat), value = T)
      
      dat[, -1] <- sapply(dat[, -1] , as.numeric)
    }
    
    return(dat)
  }
  get_cell_specific_body_matrix <- function(cell, state, mat_file) {
    # test data
    if (F) {
      cell <- "thp1"
      mat_file <- paste0("data/saved_data/11.gene_level_profile_Body_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".bin.qs")
    }
    
    # load the body profile bins
    if (T) {
      mat <- qread(file = mat_file, nthreads = 6)
    }
    
    # prepare hello info
    if (T) {
      state_order <- sort(as.numeric(unique(state[, cell])))
      start_mess <- paste0("|", paste(rep("-", 100), collapse = ""), "|\n")
      cat(start_mess)
      
      mat$num <- 1:nrow(mat)
      breaks <- round(seq(min(mat$num), max(mat$num), length.out=100))
    }
    
    # get cell specific body chromatin state matrix: absolute number
    if (T) {
      dat <- data.frame(t(
        apply(mat, 1, function(x) {
          # x <- unlist(mat[1, ])
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
            nonbody <- x[grepl("^U|TSS|TES|^D", colnames(mat))]
            nonbody <- as.numeric(nonbody)
            nonbody <- state[nonbody, cell]
            names(nonbody) <- grep("^U|TSS|TES|^D", colnames(mat), value = T)
          }
          
          # genebody
          if (T) {
            genebody_bin <- x[grepl("^B", colnames(mat))]
            genebody <- data.frame(
              sapply(genebody_bin, function(gb) {
                # gb <- genebody_bin[1]
                start <- as.numeric(strsplit(gb, "-")[[1]][1])
                end <- as.numeric(strsplit(gb, "-")[[1]][2])
                s_dat <- state[start:end, cell]
                
                s_p <- sapply(state_order, function(s) {
                  sum(s_dat==s)
                })
                names(s_p) <- paste0("S", state_order)
                return(s_p)
              })
            )
            genebody <- c(as.matrix(genebody))
            
            names(genebody) <- paste0(
              rep(grep("^B", colnames(mat), value = T), each=length(state_order)), 
              "_", 
              rep(paste0("S", state_order), body_bin_num)
            )
          }
          
          # merge the results
          if (T) {
            order <- c(
              paste0("U", up_bin_num:1), 
              "TSS", 
              names(genebody), 
              "TES", 
              paste0("D", 1:down_bin_num)
            )
            
            res <- c(nonbody, genebody)
            res <- res[order]
          }
          
          return(res)
        })
      ))
      
      dat$gene_id <- mat$gene_id
      
      dat <- dat[, c(which(grepl("gene_id", colnames(dat))), 
                     which(!grepl("gene_id", colnames(dat))))]
    }
    
    return(dat)
  }
  
  # sum matrix
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
      
      body <- melt(body, id.vars = NULL, variable.name = "ID", value.name = "value")
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
  
  # calculate the percentage: genomic percentage and state percentage
  format_cell_specific_percentage_matrix <- function(input, scale_type) {
    # test data
    if (F) {
      input <- "data/saved_data/12.gene_level_profile_TSS_5UP_5DW.bs200.stateSum_mat.cd34.qs"
      input <- "data/saved_data/12.gene_level_profile_Body_5UP_5DW_10Body.bs200.stateSum_mat.cd34.qs"
      scale_type <- "genomic"
    }
    
    # read the data
    if (T) {
      dat <- qread(input, nthreads = 6)
    }
    
    if (scale_type == "genomic") {
      dat_col <- data.frame(
        apply(dat[, grep("^U|^B|^D|TSS|TES", colnames(dat))], 2, function(x) {
          x/sum(x)*100
        })
      )
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

# data prepare: cell specific tss/body matrix for each gene
if (T) {
  for (type in c("TSS", "Body")) {
    cat(paste0(type, ":\n"))
    for (data_type in c("gene", "tx")) {
      cat(paste0("\t", data_type, ":\n"))
      for (cell in c("thp1", "cd34")) {
        cat(paste0("\t\t", cell, ": "))
        
        # test data
        if (F) {
          type <- "TSS"
          data_type <- "gene"
          cell <- "thp1"
        }
        
        # get files
        if (T) {
          mat_file <- ifelse(type == "TSS", 
                             paste0("data/saved_data/11.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, "DW.bs", bin_size, ".bin.qs"), 
                             paste0("data/saved_data/11.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, type, ".bs", bin_size, ".bin.qs"))
          
          file <- ifelse(type == "TSS", 
                         paste0("data/saved_data/12.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, "DW.bs", bin_size, ".raw_mat.", cell, ".qs"), 
                         paste0("data/saved_data/12.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, type, ".bs", bin_size, ".raw_mat.", cell, ".qs"))
        }
        
        # select function
        if (T) {
          fun <- ifelse(type == "TSS", 
                        get_cell_specific_tss_matrix, 
                        get_cell_specific_body_matrix)
        }
        
        if (! file.exists(file)) {
          cat("running...")
          
          dat <- fun(cell, state, mat_file)
          qsave(dat, file, nthreads = 6)
          
          rm(dat)
        }
        if (file.exists(file)) {
          cat("Done\n")
        }
      }
    }
  }
  
  rm(type, data_type, cell, mat_file, file, fun)
}

# sum matrix: cell specific tss/body matrix for each state, sum matrix (all genes)
if (T) {
  state_order <- sort(unique(state[, 2]))
  
  for (type in c("TSS", "Body")) {
    cat(paste0(type, ":\n"))
    for (data_type in c("gene", "tx")) {
      cat(paste0("\t", data_type, ":\n"))
      for (cell in c("thp1", "cd34")) {
        cat(paste0("\t\t", cell, ": "))
        
        # test data
        if (F) {
          type <- "TSS"
          data_type <- "gene"
          cell <- "thp1"
        }
        
        # get files
        if (T) {
          input <- ifelse(type == "TSS", 
                          paste0("data/saved_data/12.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, "DW.bs", bin_size, ".raw_mat.", cell, ".qs"), 
                          paste0("data/saved_data/12.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, type, ".bs", bin_size, ".raw_mat.", cell, ".qs"))
          
          file <- ifelse(type == "TSS", 
                         paste0("data/saved_data/12.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, "DW.bs", bin_size, ".stateSum_mat.", cell, ".qs"), 
                         paste0("data/saved_data/12.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, type, ".bs", bin_size, ".stateSum_mat.", cell, ".qs"))
        }
        
        # get gene ids
        if (T) {
          IDs <- qread(input, nthreads = 6)
          IDs <- IDs$gene_id
        }
        
        # select function
        if (T) {
          fun <- ifelse(type == "TSS", 
                        get_cell_specific_tss_matrix_sum, 
                        get_cell_specific_body_matrix_sum)
        }
        
        if (! file.exists(file)) {
          cat("running... ")
          
          dat <- fun(input, IDs, state_order)
          qsave(dat, file, nthreads = 6)
          
          rm(dat)
        }
        if (file.exists(file)) {
          cat("Done\n")
        }
      }
    }
  }
  
  rm(type, data_type, cell, input, file, IDs, fun)
}

# color setting
if (T) {
  library(RColorBrewer)
  display.brewer.all()
  colors <- brewer.pal(9, "Set1") 
  scales::show_col(colors, labels=T)
  
  col_cell <- colors[1:2]
  names(col_cell) <- c("cd34", "thp1")
}

# plot CS distribution
if (T) {
  for (type in c("TSS", "Body")) {
    cat(paste0(type, ":\n"))
    for (data_type in c("gene", "tx")) {
      cat(paste0("\t", data_type, ":\n"))
      for (scale_type in c("genomic", "state")) {
        cat(paste0("\t\t", scale_type, ": \n"))
        
        # single cell figure
        if (F) {
          for (cell in c("thp1", "cd34")) {
            cat(paste0("\t\t  plot ", cell, "'s distribution figure...  "))
            
            # test data
            if (F) {
              type <- "Body"
              data_type <- "gene"
              cell <- "thp1"
              scale_type <- "genomic"
            }
            
            # calculate the percentage
            if (T) {
              input <- ifelse(type == "TSS", 
                              paste0("data/saved_data/12.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, "DW.bs", bin_size, ".stateSum_mat.", cell, ".qs"), 
                              paste0("data/saved_data/12.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, type, ".bs", bin_size, ".stateSum_mat.", cell, ".qs"))
              
              dat <- format_cell_specific_percentage_matrix(input, scale_type)
            }
            
            # get parameters for plot
            if (T) {
              colnames(dat) <- c("State", "Loc", "Percentage")
              
              ylabel <- ifelse(scale_type == "genomic", "Genomic_Percentage", 
                               ifelse(scale_type == "state", "State_Percentage", NA))
              
              title <- paste0(cell, ": ", data_type, "+", scale_type)
              
              scale <- "free"
              
              # assistent line
              if (type == "TSS") {
                vline <- "TSS"
              } else if (type == "Body") {
                vline <- c("TSS", "TES")
              }
            }
            
            # plot the distribution
            if (T) {
              p <- ggplot(dat) +
                geom_vline(xintercept = vline, linewidth=0.8, linetype=2, alpha=0.5, color="black") +
                geom_point(aes(x=Loc, y=Percentage, color="black"), size=0.8) +
                geom_line(aes(x=Loc, y=Percentage, group=State, color=State), linewidth=1) +
                ylab(ylabel) +
                scale_x_discrete(name = NULL, breaks=c(paste0("U", up_bin_num), "TSS", "TES", paste0("D", down_bin_num))) +
                ggtitle(title) +
                facet_wrap(~State, scales = scale, ncol = 5) +
                cowplot::theme_cowplot() +
                theme(legend.position = "none",
                      strip.background = element_rect(fill=NA, color=NA))
            }
            
            # save the figure
            if (T) {
              file <- paste0("results/2.pic/12.CS_distribution_", type, "_", cell, "_", data_type, "_level_", scale_type, "_percentage.jpeg")
              ggsave(filename = file, plot = p, width = 12, height = 8)
            }
            
            cat("Done\n")
          }
        }
        
        # plot merged2cells figure
        if (T) {
          # get merged2cells data
          if (T) {
            get_dat <- function(type, cell) {
              input <- ifelse(type == "TSS", 
                              paste0("data/saved_data/12.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, "DW.bs", bin_size, ".stateSum_mat.", cell, ".qs"), 
                              paste0("data/saved_data/12.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, type, ".bs", bin_size, ".stateSum_mat.", cell, ".qs"))
              
              dat <- format_cell_specific_percentage_matrix(input, scale_type)
              dat$cell <- cell
              
              return(dat)
            }
            
            thp1 <- get_dat(type, "thp1")
            cd34 <- get_dat(type, "cd34")
            ggdat <- rbind(thp1, cd34)
            ggdat$group <- paste0(ggdat$State, ggdat$cell)
            
            rm(thp1, cd34, get_dat)
          }
          
          # get parameters for plot
          if (T) {
            colnames(ggdat) <- c("State", "Loc", "Percentage", "cell", "group")
            
            ylabel <- ifelse(scale_type == "genomic", "Genomic_Percentage", 
                             ifelse(scale_type == "state", "State_Percentage", NA))
            
            title <- paste0("merged2cells: ", data_type, "+", scale_type)
            
            scale <- "free"
            
            # assistent line
            if (type == "TSS") {
              vline <- "TSS"
            } else if (type == "Body") {
              vline <- c("TSS", "TES")
            }
          }
          
          # ggplot
          if (T) {
            p <- ggplot(ggdat) +
              geom_vline(xintercept = vline, linewidth=0.8, linetype=2, alpha=0.5, color="black") +
              geom_line(aes(x=Loc, y=Percentage, group=group, color=cell), linewidth=1, alpha=0.6) +
              ylab(ylabel) +
              scale_x_discrete(name = NULL, breaks=c(paste0("U", up_bin_num), "TSS", "TES", paste0("D", down_bin_num))) +
              scale_color_manual(values = col_cell) +
              ggtitle(title) +
              facet_wrap(~State, scales = scale, ncol = 5) +
              cowplot::theme_cowplot() +
              theme(legend.position = "none",
                    strip.background = element_rect(fill=NA, color=NA))
          }
          
          # save the figure
          if (T) {
            file <- paste0("results/2.pic/12.CS_distribution_", type, "_merged2cells_", data_type, "_level_", scale_type, "_percentage.pdf")
            ggsave(filename = file, plot = p, width = 10, height = 12)
          }
        }
      }
    }
  }
}

# merge plot
if (T) {
  # color setting
  if (T) {
    library(RColorBrewer)
    colors <- brewer.pal(9, "Set1") 
    col_cell <- colors[1:2]
    names(col_cell) <- c("cd34", "thp1")
  }
  
  # function
  if (T) {
    format_cell_specific_percentage_matrix <- function(input, scale_type) {
      # test data
      if (F) {
        input <- "data/saved_data/12.gene_level_profile_TSS_5UP_5DW.bs200.stateSum_mat.cd34.qs"
        input <- "data/saved_data/12.gene_level_profile_Body_5UP_5DW_10Body.bs200.stateSum_mat.cd34.qs"
        scale_type <- "genomic"
      }
      
      # read the data
      if (T) {
        dat <- qread(input, nthreads = 6)
      }
      
      if (scale_type == "genomic") {
        dat_col <- data.frame(
          apply(dat[, grep("^U|^B|^D|TSS|TES", colnames(dat))], 2, function(x) {
            x/sum(x)*100
          })
        )
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
  
  # get data
  if (T) {
    get_dat <- function(type, cell) {
      input <- ifelse(type == "TSS", 
                      paste0("data/saved_data/12.tx_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, "DW.bs", bin_size, ".stateSum_mat.", cell, ".qs"), 
                      paste0("data/saved_data/12.tx_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, "DW_", body_bin_num, type, ".bs", bin_size, ".stateSum_mat.", cell, ".qs"))
      
      dat <- format_cell_specific_percentage_matrix(input, "genomic")
      dat$cell <- cell
      
      return(dat)
    }
    
    thp1_tss <- get_dat("TSS", "thp1")
    cd34_tss <- get_dat("TSS", "cd34")
    tss <- rbind(thp1_tss, cd34_tss)
    thp1_body <- get_dat("Body", "thp1")
    cd34_body <- get_dat("Body", "cd34")
    body <- rbind(thp1_body, cd34_body)
    
    tss$type <- "TSS"
    body$type <- "Body"
    
    ggdat <- rbind(tss, body)
    
    rm(thp1_tss, cd34_tss, tss, thp1_body, cd34_body, body, get_dat)
  }
  
  # add outliner
  if (T) {
    head(ggdat)
    
    outline <- as.data.frame(tapply(ggdat$Genomic_Percentage, ggdat$State, max))
    outline$State <- rownames(outline)
    colnames(outline)[1] <- "Genomic_Percentage"
    outline$Genomic_Percentage <- outline$Genomic_Percentage*1.2
    outline$Loc <- "D5"
    outline$cell <- "outline"
    outline$type <- "outline"
    outline$group <- "outline"
    outline <- outline[, colnames(ggdat)]
    ggdat <- rbind(ggdat, outline)
  }
  
  # format the data
  if (T) {
    head(ggdat)
    ggdat$group <- paste0(ggdat$cell, "@", ggdat$State, "@", ggdat$type)
    
    ggdat$Loc <- ifelse(ggdat$type == "TSS", paste0("T_", ggdat$Loc), paste0("B_", ggdat$Loc))
    ggdat$Loc <- factor(ggdat$Loc, levels = c(
      paste0("T_", c(
        paste0("U", up_bin_num:1), "TSS", paste0("D", 1:down_bin_num)
      )), 
      paste0("B_", c(
        paste0("U", up_bin_num:1), "TSS", paste0("B", 1:body_bin_num), "TES", paste0("D", 1:down_bin_num)
      ))
    ))
    
    levels(ggdat$Loc)
    sort(unique(ggdat$Loc))
    
    vline <- c("T_TSS", "B_TSS", "B_TES")
  }
  
  head(ggdat)
  
  p <- ggplot(ggdat) +
    geom_line(aes(x=Loc, y=Genomic_Percentage, group=group, color=cell), linewidth=1, alpha=0.6) +
    geom_vline(xintercept = vline, linewidth=0.8, linetype=1, alpha=0.5, color="black") +
    scale_x_discrete(name = NULL, breaks=c(paste0("U", up_bin_num), "TSS", "TES", paste0("D", down_bin_num))) +
    scale_color_manual(values = col_cell) +
    facet_wrap(~State, scales = "free", ncol = 5) +
    theme_bw() +
    theme(legend.position = "none",
          strip.background = element_rect(fill=NA, color=NA))
  
  ggsave(filename = "results/2.pic/12.CS_distribution_merged_TSS_Body.pdf", plot = p, width = 10, height = 12)
}

# manual
if (F) {
  # line1-23: 读入IDEAS得到的state坐标信息
  
  # line24-456: 定义后续用到的所有函数，实际上是该脚本中的核心部分
  #   line26-274: 定义获取TSS以及GeneBody区域上不同染色质状态的总数
  #     line27-83: 定义TSS区域上不同染色质状态的总数
  #       (1) 读入21脚本中获得所有基因在TSS区域的BinID，并根据BinID替换成对应区域上细胞特异性的染色质状态，得到变量dat
  #       (2) 根据变量dat，统计在TSS区域上，每种染色质状态在所有基因每个区域中的所有数量加和，同时在最后一列添加每种染色质状态的总数，得到变量state_profile
  #     line84-273: 定义GeneBody区域上下游不同染色质状态的总数
  #       (1) 读入21脚本中获得所有基因在GeneBody区域上下游的BinID，并根据BinID替换成对应区域上细胞特异性的染色质状态，由于这里同时存在2个区域数据，故分开进行处理：
  #           - 在nonbody区域（上游至TSS以及TES至下游）：每种染色质状态在所有基因每个区域中的所有数量加和
  #           - 在body区域：因为被切分成body_bin_num块，GeneBody区域内以GenePart_State的形式进行记录
  #           - 将上述2个数据以cbind进行合并得到变量dat
  #       (2) 考虑到在GeneBody区域内的数据会收到bin数目的影响，故根据变量dat，进行一系列数据统计：
  #           - 变量state_sum：统计在GeneBody区域上下游，每种染色质状态在所有基因每个区域中的所有数量加和，同时在最后一列添加每种染色质状态的总数
  #           - 变量bin_sum：统计在GeneBody区域上下游，每个区域（从上游到GeneBody区域到下游）中bin数目的总和
  
  #   line275-389: 为后续绘制不同细胞中每种染色质状态在基因组上分布准备数据，分别准备以基因组为总数和以特定state为总数的数据，以state1为例：
  #           - 基因组为总数(scale_type == genomic)：在全部TSS上区域中，有20%的区域是为state1
  #           - 以特定state为总数(scale_type == state)：在全部state1分布的区域中，有20%的state1分布在TSS上
  #     line276-328:处理TSS区域上的数据
  #       (1) 读入get_cell_specific_tss_matrix函数制作的变量state_profile
  #       (2) 对数据进行scale，分为genomic和state类型：
  #           - scale_type == genomic：以特定区域中总bin数作为总体，计算每个state在该位置上的占比
  #           - scale_type == state：以每个state总bin数作为总体，计算在每个位置上该state的占比
  #     line329-388: 处理GeneBody区域上下游上的数据
  #       (1) 读入get_cell_specific_body_matrix函数制作的变量state_sum
  #       (2) 对数据进行scale，分为genomic和state类型：
  #           - scale_type == genomic：以特定区域中总bin数作为总体，计算每个state在该位置上的占比，对于GeneBody区域无需额外处理
  #           - scale_type == state：以每个state总bin数作为总体，计算在每个位置上该state的占比，对于GeneBody区域，需要考虑在该区域中的bin数目比上下游区域的数目多，
  #                                  故需要除以一个系数来消除该差异。系数的计算方法，以上下游区域bin数目为1，计算每个GeneBody区域平均bin数目是上下游区域bin数目的倍数
  
  #   line390-455: 根据格式化后的数据进行绘制TSS区域以及GeneBody区域上下游上染色质状态的分布模式图
  #     line391-409: 绘制染色质状态的分布模式，对于TSS以及GeneBody区域上下游的数据通用
  #     line410-454: 绘制GeneBody区域上下游的分布
}
