# copy of script 37

# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(stringr)
  library(reshape2)
  library(qs)
  library(ggridges)
  suppressPackageStartupMessages(library(circlize))
}

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
up_bin_num <- 5
down_bin_num <- 5
body_bin_num <- 10
rna_levels <- 2

# plot the CS distribution
if (T) {
  data_type <- "tx"
  
  for (type in c("Body")) {
    cat(paste0(type, ":\n"))
    for (cell in c(cell1, cell2, "merged")) {
      cat(paste0("\t", cell, ": \n"))
      for (cs_type in c("CS", "CSC")) {
        cat(paste0("\t\t", cs_type, ": \n"))
        
        # test data
        if (F) {
          type <- "Body"
          cell <- "merged"
          cell <- cell1
          cs_type <- "CS"
        }
        
        # get input
        if (T) {
          input1 <- ifelse(type == "TSS", 
                           paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                  "DW.bs", bin_size, ".raw_mRNA_classification_Q", rna_levels, "_mat.", cell, ".qs"), 
                           paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                  "DW_", body_bin_num, "Body.bs", bin_size, ".raw_mRNA_classification_Q", rna_levels, "_mat.", cell, ".qs"))
          
          input2 <- ifelse(type == "TSS", 
                           paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                  "DW.bs", bin_size, ".format_CSC_distribution_matrix_mRNA_Q", 
                                  rna_levels, ".", cell, ".qs"), 
                           paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                  "DW_", body_bin_num, "Body.bs", bin_size, ".format_CSC_distribution_matrix_mRNA_Q", 
                                  rna_levels, ".", cell, ".qs"))
          
          input <- ifelse(cs_type == "CS", input1, input2)
          
          rm(input1, input2)
        }
        
        # read and transition into percentage
        if (T) {
          dat <- qread(input, nthreads = 6)
          
          if (class(dat) == "list") {
            format_dat <- lapply(dat, function(subdat) {
              subdat[, 1:(up_bin_num + 1 + body_bin_num + 1 + down_bin_num)] <- sapply(subdat[, 1:(up_bin_num + 1 + body_bin_num + 1 + down_bin_num)], function(x) {
                x/sum(x)*100
              })
              
              return(subdat)
            })
            
            format_dat <- do.call(rbind, format_dat)
            format_dat$mRNA <- str_split(rownames(format_dat), "[.]", simplify = T)[, 1]
          }
          if (class(dat) == "data.frame") {
            format_dat <- dat
            for (rna in paste0("Q", 1:rna_levels)) {
              # rna <- "Q1"
              format_dat[format_dat$mRNA == rna, 1:(up_bin_num + 1 + body_bin_num + 1 + down_bin_num)] <- sapply(format_dat[format_dat$mRNA == rna, 1:(up_bin_num + 1 + body_bin_num + 1 + down_bin_num)], function(x) {
                x/sum(x)*100
              })
              
              rm(rna)
            }
            
            colnames(format_dat)[ncol(format_dat)-1] <- "State"
          }
          format_dat$mRNA <- factor(format_dat$mRNA, levels = paste0("Q", 1:rna_levels))
        }
        
        # wide 2 long
        if (T) {
          dat <- melt(format_dat, id.vars = c("State", "mRNA"), variable.name = "Loc", value.name = "Percentage")
          dat$State <- ifelse(grepl("S", dat$State), dat$State, paste0("S", dat$State))
          dat$State <- factor(dat$State, levels = paste0("S", sort(unique(as.numeric(gsub("S", "", dat$State))))))
          
          dat$mRNA <- factor(dat$mRNA, levels = paste0("Q", 1:rna_levels))
          dat$Loc <- factor(dat$Loc, levels = c(paste0("U", up_bin_num:1), 
                                                "TSS", 
                                                paste0("B", 1:body_bin_num), 
                                                "TES", 
                                                paste0("D", 1:down_bin_num)))
        }
        
        # color setting
        if (T) {
          col_fun <- colorRamp2(c(1, rna_levels), c("blue", "red"))
          cols <- col_fun(1:rna_levels)
          names(cols) <- paste0("Q", 1:rna_levels)
        }
        
        # ggplot
        if (T) {
          head(dat)
          p <- ggplot(dat) +
            geom_line(aes(x=Loc, y=Percentage, group=mRNA, color=mRNA), alpha=0.6, linewidth=1) +
            scale_x_discrete(name = NULL, breaks=c(paste0("U", up_bin_num), "TSS", "TES", paste0("D", down_bin_num))) +
            scale_color_manual(values = cols) +
            ylab("Genomic Percentage (%)") +
            ggtitle(cell) +
            facet_wrap(.~State, scales = "free", ncol = 5) +
            cowplot::theme_cowplot() +
            theme(panel.border = element_rect(color="black"), 
                  strip.background = element_rect(fill=NA, color=NA))
        }
        
        # get output
        if (T) {
          file <- paste0("results/2.pic/77.", cs_type, "_distribution_", cell, "_", data_type, "_level_mRNA_Q", rna_levels, ".pdf")
          
          if (cs_type == "CS") {
            ggsave(file, plot = p, width = 9, height = 11)
          }
          if (cs_type == "CSC") {
            ggsave(file, plot = p, width = 9, height = 2.5)
          }
        }
      }
    }
  }
}

# get order from the distribution
if (T) {
  data_type <- "tx"
  
  for (type in c("Body")) {
    cat(paste0(type, ":\n"))
    for (cell in c(cell1, cell2, "merged")) {
      cat(paste0("\t", cell, ": \n"))
      for (cs_type in c("CS", "CSC")) {
        cat(paste0("\t\t", cs_type, ": \n"))
        
        # test data
        if (F) {
          type <- "Body"
          cell <- "merged"
          cell <- cell1
          cs_type <- "CS"
          cs_type <- "CSC"
        }
        
        # get input
        if (T) {
          input1 <- ifelse(type == "TSS", 
                           paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                  "DW.bs", bin_size, ".raw_mRNA_classification_Q", rna_levels, "_mat.", cell, ".qs"), 
                           paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                  "DW_", body_bin_num, "Body.bs", bin_size, ".raw_mRNA_classification_Q", rna_levels, "_mat.", cell, ".qs"))
          
          input2 <- ifelse(type == "TSS", 
                           paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                  "DW.bs", bin_size, ".format_CSC_distribution_matrix_mRNA_Q", 
                                  rna_levels, ".", cell, ".qs"), 
                           paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                  "DW_", body_bin_num, "Body.bs", bin_size, ".format_CSC_distribution_matrix_mRNA_Q", 
                                  rna_levels, ".", cell, ".qs"))
          
          input <- ifelse(cs_type == "CS", input1, input2)
          
          rm(input1, input2)
        }
        
        # read and transition into percentage and scale by row
        if (T) {
          dat <- qread(input, nthreads = 6)
          
          if (class(dat) == "list") {
            format_dat <- lapply(dat, function(subdat) {
              subdat[, 1:(up_bin_num + 1 + body_bin_num + 1 + down_bin_num)] <- sapply(subdat[, 1:(up_bin_num + 1 + body_bin_num + 1 + down_bin_num)], function(x) {
                x/sum(x)*100
              })
              
              subdat[, 1:(up_bin_num + 1 + body_bin_num + 1 + down_bin_num)] <- t(scale(t(subdat[, 1:(up_bin_num + 1 + body_bin_num + 1 + down_bin_num)])))
              
              return(subdat)
            })
            
            format_dat <- do.call(rbind, format_dat)
            format_dat$mRNA <- str_split(rownames(format_dat), "[.]", simplify = T)[, 1]
          }
          if (class(dat) == "data.frame") {
            format_dat <- dat
            for (rna in paste0("Q", 1:rna_levels)) {
              # rna <- "Q1"
              format_dat[format_dat$mRNA == rna, 1:(up_bin_num + 1 + body_bin_num + 1 + down_bin_num)] <- sapply(format_dat[format_dat$mRNA == rna, 1:(up_bin_num + 1 + body_bin_num + 1 + down_bin_num)], function(x) {
                x/sum(x)*100
              })
              
              format_dat[format_dat$mRNA == rna, 1:(up_bin_num + 1 + body_bin_num + 1 + down_bin_num)] <- t(scale(t(format_dat[format_dat$mRNA == rna, 1:(up_bin_num + 1 + body_bin_num + 1 + down_bin_num)])))
              rm(rna)
            }
            
            colnames(format_dat)[ncol(format_dat)-1] <- "State"
          }
          format_dat$mRNA <- factor(format_dat$mRNA, levels = paste0("Q", 1:rna_levels))
        }
        
        # wide 2 long
        if (T) {
          dat <- melt(format_dat, id.vars = c("State", "mRNA"), variable.name = "Loc", value.name = "Percentage")
          dat$State <- ifelse(grepl("S", dat$State), dat$State, paste0("S", dat$State))
          dat$State <- factor(dat$State, levels = paste0("S", sort(unique(as.numeric(gsub("S", "", dat$State))))))
          
          dat$mRNA <- factor(dat$mRNA, levels = paste0("Q", 1:rna_levels))
          dat$Loc <- factor(dat$Loc, levels = c(paste0("U", up_bin_num:1), 
                                                "TSS", 
                                                paste0("B", 1:body_bin_num), 
                                                "TES", 
                                                paste0("D", 1:down_bin_num)))
        }
        
        # color setting
        if (T) {
          col_fun <- colorRamp2(c(1, rna_levels), c("blue", "red"))
          cols <- col_fun(1:rna_levels)
          names(cols) <- paste0("Q", 1:rna_levels)
        }
        
        # ggplot
        if (T) {
          head(dat)
          p <- ggplot(dat) +
            geom_line(aes(x=Loc, y=Percentage, group=mRNA, color=mRNA), alpha=0.6, linewidth=1) +
            scale_x_discrete(name = NULL, breaks=c(paste0("U", up_bin_num), "TSS", "TES", paste0("D", down_bin_num))) +
            scale_color_manual(values = cols) +
            ylab("Genomic Percentage (%)") +
            ggtitle(cell) +
            facet_wrap(.~State, scales = "free", ncol = 5) +
            cowplot::theme_cowplot() +
            theme(panel.border = element_rect(color="black"), 
                  strip.background = element_rect(fill=NA, color=NA))
          print(p)
        }
        
        # get output
        if (T) {
          file <- paste0("results/2.pic/77.", cs_type, "_distribution_", cell, "_", data_type, "_level_mRNA_Q", rna_levels, ".pdf")
          
          if (cs_type == "CS") {
            ggsave(file, plot = p, width = 9, height = 11)
          }
          if (cs_type == "CSC") {
            ggsave(file, plot = p, width = 9, height = 2.5)
          }
        }
      }
    }
  }
}

# sort the distribution line based on the pseudoLoc
if (T) {
  data_type <- "tx"
  
  for (type in c("Body")) {
    cat(paste0(type, ":\n"))
    for (cell in c(cell1, cell2, "merged")) {
      cat(paste0("\t", cell, ": \n"))
      # test data
      if (F) {
        type <- "Body"
        cell <- "merged"
      }
      
      # get input
      if (T) {
        input <- ifelse(type == "TSS", 
                        paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                               "DW.bs", bin_size, ".raw_mRNA_classification_Q", rna_levels, "_mat.", cell, ".qs"), 
                        paste0("data/saved_data/76.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                               "DW_", body_bin_num, "Body.bs", bin_size, ".raw_mRNA_classification_Q", rna_levels, "_mat.", cell, ".qs"))
      }
      
      # read and transition into percentage and scale by row
      if (T) {
        dat <- qread(input, nthreads = 6)
        
        format_dat <- lapply(dat, function(subdat) {
          subdat[, 1:(up_bin_num + 1 + body_bin_num + 1 + down_bin_num)] <- sapply(subdat[, 1:(up_bin_num + 1 + body_bin_num + 1 + down_bin_num)], function(x) {
            x/sum(x)*100
          })
          
          subdat[, 1:(up_bin_num + 1 + body_bin_num + 1 + down_bin_num)] <- t(scale(t(subdat[, 1:(up_bin_num + 1 + body_bin_num + 1 + down_bin_num)])))
          
          return(subdat)
        })
        format_dat <- do.call(rbind, format_dat)
        
        format_dat$mRNA <- str_split(rownames(format_dat), "[.]", simplify = T)[, 1]
        format_dat$mRNA <- factor(format_dat$mRNA, levels = paste0("Q", 1:rna_levels))
      }
      
      # load pseudo_Loc
      if (T) {
        pseudoLoc_active <- read.table(paste0("results/1.tab/75.pseudoLoc_active.", cell, ".txt"), header = T, sep = "\t", fill = T, comment.char = "")
        pseudoLoc_repressive <- read.table(paste0("results/1.tab/75.pseudoLoc_repressive.", cell, ".txt"), header = T, sep = "\t", fill = T, comment.char = "")
        
        identical(pseudoLoc_active$state, pseudoLoc_repressive$state)
        pseudoLoc <- data.frame(
          state = pseudoLoc_active$state, 
          pseudoloc = ifelse(is.infinite(pseudoLoc_active$loc), pseudoLoc_repressive$loc, pseudoLoc_active$loc), 
          type = ifelse(is.infinite(pseudoLoc_active$loc), "repressive", "active")
        )
        
        rm(pseudoLoc_active, pseudoLoc_repressive)
      }
      
      # merge distribution info into pseudoLoc
      if (T) {
        pseudoLoc <- lapply(paste0("Q", 1:rna_levels), function(rna) {
          # rna <- "Q1"
          subdat <- format_dat[format_dat$mRNA == rna, ]
          rna <- ifelse(rna == "Q2", "active", "repressive")
          
          pseudoLoc_subdat <- pseudoLoc[pseudoLoc$type == rna, ]
          subdat <- subdat[match(pseudoLoc_subdat$state, subdat$State), ]
          
          pseudoLoc_subdat <- cbind(pseudoLoc_subdat, subdat)
          return(pseudoLoc_subdat) 
        })
        pseudoLoc <- data.frame(do.call(rbind, pseudoLoc))
      }
      
      # wide 2 long
      if (T) {
        head(pseudoLoc)
        dat <- melt(pseudoLoc, id.vars = c("state", "type", "pseudoloc"), 
                    variable.name = "Loc", value.name = "Percentage", 
                    measure.vars = c(paste0("U", up_bin_num:1), 
                                     "TSS", 
                                     paste0("B", 1:body_bin_num), 
                                     "TES", 
                                     paste0("D", 1:down_bin_num)))
        dat$state <- ifelse(grepl("S", dat$state), dat$state, paste0("S", dat$state))
        dat$Loc <- factor(dat$Loc, levels = c(paste0("U", up_bin_num:1), 
                                              "TSS", 
                                              paste0("B", 1:body_bin_num), 
                                              "TES", 
                                              paste0("D", 1:down_bin_num)))
      }
      
      # all transition into positive value
      if (T) {
        dat$Percentage <- dat$Percentage - min(dat$Percentage)
      }
      
      # subdat selection
      if (T) {
        active_dat <- dat[dat$type == "active", ]
        repressive_dat <- dat[dat$type == "repressive", ]
        
        ord <- pseudoLoc[pseudoLoc$type == "active", ]
        ord <- ord[order(ord$pseudoloc), ]
        active_dat$state <- factor(active_dat$state, levels = c(ord$state))
        
        ord <- pseudoLoc[pseudoLoc$type == "repressive", ]
        ord <- ord[order(ord$pseudoloc), ]
        if (cell != "cd34") {
          ord <- ord[order(ord$pseudoloc, decreasing = T), ]
        }
        repressive_dat$state <- factor(repressive_dat$state, levels = c(ord$state))
      }
      
      # ggplot
      if (T) {
        head(repressive_dat)
        head(active_dat)
        
        p <- ggplot(active_dat) +
          geom_ridgeline(aes(x=Loc, y=state, height=Percentage, group=state, fill=pseudoloc), alpha=0.8, linewidth=1) + 
          scale_x_discrete(name = NULL, breaks=c(paste0("U", up_bin_num), "TSS", "TES", paste0("D", down_bin_num))) +
          scale_fill_gradientn(colours=c("#FFFFBF", "#FEE090", "#FC8D59", "#D73027")) +
          ylab("Relative Enrichment Level") +
          ggtitle(cell) +
          cowplot::theme_cowplot() +
          theme(panel.border = element_rect(color="black"), 
                strip.background = element_rect(fill=NA, color=NA))
        ggsave(paste0("results/2.pic/77.CS_distribution_based_on_pseudoLoc_", data_type, "_level_active.", cell, ".pdf"), 
               plot = p, width = 10, height = 8)
        
        p <- ggplot(repressive_dat) +
          geom_ridgeline(aes(x=Loc, y=state, height=Percentage, group=state, fill=pseudoloc), alpha=0.8, linewidth=1) + 
          scale_x_discrete(name = NULL, breaks=c(paste0("U", up_bin_num), "TSS", "TES", paste0("D", down_bin_num))) +
          scale_fill_gradientn(colours=rev(c("#E0F3F8", "#91BFDB", "#4575B4"))) +
          ylab("Relative Enrichment Level") +
          ggtitle(cell) +
          cowplot::theme_cowplot() +
          theme(panel.border = element_rect(color="black"), 
                strip.background = element_rect(fill=NA, color=NA))
        ggsave(paste0("results/2.pic/77.CS_distribution_based_on_pseudoLoc_", data_type, "_level_repressive.", cell, ".pdf"), 
               plot = p, width = 10, height = 8)
      }
    }
  }
}
