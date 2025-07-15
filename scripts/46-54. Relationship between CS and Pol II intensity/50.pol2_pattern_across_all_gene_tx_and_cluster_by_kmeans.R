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

# define functions
if (T) {
  # common value
  if (T) {
    name <- c(
      paste0("U", up_bin_num:1), 
      "TSS", 
      paste0("B", 1:body_bin_num), 
      "TES", 
      paste0("D", 1:down_bin_num)
    )
    
    labels <- c(paste0("U", up_bin_num), 
                "TSS", 
                "TES", 
                paste0("D", down_bin_num))
  }
  
  # k-means (modify based on script40)
  cluster_kmeans <- function(mat, k, scale) {
    # test dat
    if (F) {
      k <- 4
      scale <- "row"
    }
    
    df <- mat
    
    # scale data
    if (T) {
      if (scale == "row") {
        df <- data.frame(t(scale(t(df))))
        torf <- apply(df, 1, function(x) {
          sum(is.na(x)) == 0
        })
        df <- df[torf, ]
      }
      if (scale == "col") {
        df <- data.frame(scale(df))
        torf <- apply(df, 1, function(x) {
          sum(is.na(x)) == 0
        })
        df <- df[torf, ]
      }
      if (scale == "none") {
        df <- df
      }
    }
    
    # k-means
    if (T) {
      set.seed(799)
      fit.km <- kmeans(df, centers=k, nstart=25, iter.max=50, algorithm="Lloyd")
    }
    
    return(fit.km)
  }
  determineK <- function(dat, title) {
    # test data
    if (F) {
      title <- "title"
    }
    
    # merge the info
    if (T) {
      kdat <- sapply(dat, function(x) {
        sum(x$withinss)
      })
      wss_dat <- data.frame(K = names(kdat), 
                            WSS = as.numeric(kdat))
      wss_dat$Percentage <- wss_dat$WSS/max(wss_dat$WSS)*100
    }
    
    # format the data
    if (T) {
      wss_dat$K <- factor(wss_dat$K, levels = paste0("K", 1:15))
    }
    
    # ggplot
    if (T) {
      head(wss_dat)
      
      p <- ggplot(wss_dat) +
        geom_line(aes(x=K, y=Percentage), group=1, linewidth=1) +
        scale_x_discrete(name = NULL, breaks = paste0("K", seq(1, 15, 2))) +
        ylab("Within groups sum of squares") +
        ggtitle(title) +
        theme(axis.text = element_text(size = rel(1.2)),
              axis.title = element_text(size = rel(1.2)),
              legend.text = element_text(size = rel(1.2)),
              legend.title = element_text(size = rel(1.2)),
              strip.text = element_text(size = rel(1.2)))
    }
    
    return(p)
  }
  
  # rename the group label
  pattern_group_pk_location <- function(submat, cluster, average_type="mean") {
    # merge the signal
    if (T) {
      profile <- submat
      profile$cluster <- cluster$cluster
    }
    
    # get average signal: mean
    if (T) {
      if (average_type == "mean") {
        fun <- mean
      }
      if (average_type == "median") {
        fun <- median
      }
      
      ggdat <- lapply(split(profile, profile$cluster), function(subdat) {
        sapply(subdat[, name], fun)
      })
      
      ggdat <- data.frame(do.call(rbind, ggdat))
    }
    
    # get mRNA and cluster info
    if (T) {
      ggdat$cluster <- paste0("pol2_", rownames(ggdat))
    }
    
    # statistics pol2 pk location
    if (T) {
      pk_stat <- sapply(1:k, function(cluster_k) {
        # cluster_k <- 1
        subdat <- ggdat[ggdat$cluster == paste0("pol2_", cluster_k), name]
        pk_loc <- names(subdat)[which.max(subdat)]
        pk_loc <- which(name %in% pk_loc)
        
        return(mean(pk_loc))
      })
      names(pk_stat) <- paste0("pol2_", 1:k)
    }
    
    return(pk_stat)
  }
  
  # get average signal
  average_sig <- function(profile, average_type="mean") {
    # get average signal
    if (T) {
      if (average_type == "mean") {
        fun <- mean
      }
      if (average_type == "median") {
        fun <- median
      }
      
      ggdat <- lapply(split(profile, profile$cluster), function(subdat) {
        sapply(subdat[, name], fun)
      })
      
      ggdat <- data.frame(do.call(rbind, ggdat))
    }
    
    # get pol2 cluster info
    if (T) {
      ggdat$cluster <- paste0("pol2_", rownames(ggdat))
    }
    
    # format 
    if (T) {
      ggdat <- melt(ggdat, id.vars = c("cluster"), measure.vars = name, variable.name = "Loc", value.name = "Signal")
      ggdat$cluster <- factor(ggdat$cluster, levels = paste0("pol2_", 1:k))
      ggdat$Loc <- factor(ggdat$Loc, levels = name)
      
      ggdat$Assistant <- ifelse(ggdat$Loc %in% labels, ggdat$Signal, NA)
    }
    
    return(ggdat)
  }
  
  move_file <- function(file, dir) {
    stat <- file.copy(from = file, to = paste0(dir, "/"), overwrite = T, copy.date = T)
    if (stat) {
      stat <- file.remove(file)
    }
  }
}

# apply kmeans to pol2 pattern/signal dat (k: 1-15)
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir1 <- paste0("data/saved_data/50.pol2_pattern/", data_type)
      if (! dir.exists(dir1)) {
        dir.create(dir1, showWarnings = F, recursive = T)
      }
      
      dir2 <- paste0("results/2.pic/50.pol2_pattern/", data_type)
      if (! dir.exists(dir2)) {
        dir.create(dir2, showWarnings = F, recursive = T)
      }
    }
    
    for (type in c("merged", "single")) {
      cat(paste0("\t", type, ": \n"))
      
      # cell specific data
      if (type == "single") {
        for (cell in c(cell1, cell2)) {
          cat(paste0("\t\t", cell, ": \n"))
          
          # test data
          if (F) {
            data_type <- "tx"
            type <- "single"
            cell <- cell1
          }
          
          # get pol2 matrix: mat
          if (T) {
            cell_dat <- function(cell, type) {
              input <- paste0("data/saved_data/48.mk_signal_distribution/tx/48.", data_type, "_Body_", up_bin_num, "UP_", 
                              down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".pol2.mean_mat.", cell, ".", type, ".qs")
              
              mat <- qread(input, nthreads = 6)
              rownames(mat) <- paste0(mat$gene_id, "@", cell)
              
              return(mat)
            }
            mat <- cell_dat(cell, type)
          }
          
          # get kmeans clusters: pat dat
          if (T) {
            file <- paste0(dir1, "/50.pol2_pat_bs", bin_size, "_kmeans1-15_", data_type, "_", cell, ".", type, ".qs")
            
            # clustered by pol2 distribution **pattern**, not the signal level
            if (! file.exists(file)) {
              dat <- lapply(1:15, function(k) {
                # k <- 2
                cat(paste0("Now the K of Kmeans: ", k, "\n"))
                cluster_kmeans(mat[, name], k, scale="row")
              })
              names(dat) <- paste0("K", 1:15)
              
              qsave(dat, file, nthreads = 6)
            }
            if (file.exists(file)) {
              dat <- qread(file, nthreads = 6)
            }
          }
          
          # determine the optimal k value
          if (T) {
            p <- determineK(dat, title=paste0(type, ": ", cell))
            
            file <- paste0("50.", data_type, "s_pol2_pat_determine_Kvalue.bs", bin_size, ".", cell, ".", type, ".pdf")
            ggsave(filename = file, plot = p, width = 8, height = 6)
            move_file(file, dir2)
            
            rm(p)
          }
          
          # get kmeans clusters: sig dat
          if (T) {
            file <- paste0(dir1, "/50.pol2_sig_bs", bin_size, "_kmeans1-15_", data_type, "_", cell, ".", type, ".qs")
            
            # clustered by pol2 distribution **pattern**, not the signal level
            if (! file.exists(file)) {
              dat <- lapply(1:15, function(k) {
                # k <- 2
                cat(paste0("Now the K of Kmeans: ", k, "\n"))
                cluster_kmeans(mat[, name], k, scale="col")
              })
              names(dat) <- paste0("K", 1:15)
              
              qsave(dat, file, nthreads = 6)
            }
            if (file.exists(file)) {
              dat <- qread(file, nthreads = 6)
            }
          }
          
          # determine the optimal k value
          if (T) {
            p <- determineK(dat, title=paste0(type, ": ", cell))
            
            file <- paste0("50.", data_type, "s_pol2_sig_determine_Kvalue.bs", bin_size, ".", cell, ".", type, ".pdf")
            ggsave(filename = file, plot = p, width = 8, height = 6)
            move_file(file, dir2)
            
            rm(p)
          }
        }
      }
      
      # merged data
      if (type == "merged") {
        cat(paste0("\t\tmerged: \n"))
        
        # test data
        if (F) {
          data_type <- "tx"
          type <- "merged"
        }
        
        # get pol2 matrix: mat
        if (T) {
          cell_dat <- function(cell, type) {
            input <- paste0("data/saved_data/48.mk_signal_distribution/tx/48.", data_type, "_Body_", up_bin_num, "UP_", 
                            down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".pol2.mean_mat.", cell, ".", type, ".qs")
            
            mat <- qread(input, nthreads = 6)
            rownames(mat) <- paste0(mat$gene_id, "@", cell)
            
            return(mat)
          }
          dat1 <- cell_dat(cell1, type="single")
          dat2 <- cell_dat(cell2, type="single")
          
          mat <- rbind(dat1, dat2)
        }
        
        # get kmeans clusters: pat dat
        if (T) {
          file <- paste0(dir1, "/50.pol2_pat_bs", bin_size, "_kmeans1-15_", data_type, "_merged.", type, ".qs")
          
          # clustered by pol2 distribution **pattern**, not the signal level
          if (! file.exists(file)) {
            dat <- lapply(1:15, function(k) {
              # k <- 2
              cat(paste0("Now the K of Kmeans: ", k, "\n"))
              cluster_kmeans(mat[, name], k, scale="row")
            })
            names(dat) <- paste0("K", 1:15)
            
            qsave(dat, file, nthreads = 6)
          }
          if (file.exists(file)) {
            dat <- qread(file, nthreads = 6)
          }
        }
        
        # determine the optimal k value
        if (T) {
          p <- determineK(dat, title=paste0(type, ": merged"))
          
          file <- paste0("50.", data_type, "s_pol2_pat_determine_Kvalue.bs", bin_size, ".merged.", type, ".pdf")
          ggsave(filename = file, plot = p, width = 8, height = 6)
          move_file(file, dir2)
          
          rm(p)
        }
        
        # get kmeans clusters: sig dat
        if (T) {
          file <- paste0(dir1, "/50.pol2_sig_bs", bin_size, "_kmeans1-15_", data_type, "_merged.", type, ".qs")
          
          # clustered by pol2 distribution **pattern**, not the signal level
          if (! file.exists(file)) {
            dat <- lapply(1:15, function(k) {
              # k <- 2
              cat(paste0("Now the K of Kmeans: ", k, "\n"))
              cluster_kmeans(mat[, name], k, scale="col")
            })
            names(dat) <- paste0("K", 1:15)
            
            qsave(dat, file, nthreads = 6)
          }
          if (file.exists(file)) {
            dat <- qread(file, nthreads = 6)
          }
        }
        
        # determine the optimal k value
        if (T) {
          p <- determineK(dat, title=paste0(type, ": merged"))
          
          file <- paste0("50.", data_type, "s_pol2_sig_determine_Kvalue.bs", bin_size, ".merged.", type, ".pdf")
          ggsave(filename = file, plot = p, width = 8, height = 6)
          move_file(file, dir2)
          
          rm(p)
        }
      }
    }
  }
}

# get cluster when K in 2:6 (pat row)
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir1 <- paste0("data/saved_data/50.pol2_pattern/", data_type)
      if (! dir.exists(dir1)) {
        dir.create(dir1, showWarnings = F, recursive = T)
      }
      
      dir2 <- paste0("results/2.pic/50.pol2_pattern/", data_type)
      if (! dir.exists(dir2)) {
        dir.create(dir2, showWarnings = F, recursive = T)
      }
    }
    
    for (type in c("merged", "single")) {
      cat(paste0("\t", type, ": \n"))
      
      # cell specific data
      if (type == "single") {
        for (cell in c(cell1, cell2)) {
          cat(paste0("\t\t", cell, ": \n"))
          
          # test data
          if (F) {
            data_type <- "tx"
            type <- "single"
            cell <- cell1
          }
          
          # read kmeans result
          if (T) {
            dat <- qread(paste0(dir1, "/50.pol2_pat_bs", bin_size, "_kmeans1-15_", data_type, "_", cell, ".", type, ".qs"), nthreads = 6)
          }
          
          # get pol2 matrix: mat
          if (T) {
            cell_dat <- function(cell, type) {
              input <- paste0("data/saved_data/48.mk_signal_distribution/tx/48.", data_type, "_Body_", up_bin_num, "UP_", 
                              down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".pol2.mean_mat.", cell, ".", type, ".qs")
              
              mat <- qread(input, nthreads = 6)
              rownames(mat) <- paste0(mat$gene_id, "@", cell)
              
              return(mat)
            }
            mat <- cell_dat(cell, type)
          }
          
          # kmeans: k2-6
          for (k in 2:6) {
            # k <- 5
            cat(paste0("\t\tk=", k, "\n"))
            
            # get cluster group
            if (T) {
              # get kmeans cluster data: cluster
              if (T) {
                cluster <- dat[[paste0("K", k)]]
                cluster <- cluster$cluster
                cluster <- data.frame(ID = names(cluster), 
                                      cluster = as.numeric(cluster))
              }
              
              # get consistent mat: submat
              if (T) {
                submat <- mat[cluster$ID, ]
              }
              
              file <- paste0(dir1, "/50.pol2_pat_bs", bin_size, "_kmeans", k, "_", data_type, "_", cell, ".", type, ".qs")
              if (! file.exists(file)) {
                # get pol2 pat pk location
                if (T) {
                  submat_row <- data.frame(t(scale(t(submat[, name]))))
                  pk_stat <- pattern_group_pk_location(submat_row, cluster, average_type="mean")
                }
                
                # rename the cluster label: assign the cluster label based on the pol2 peak location
                if (T) {
                  # change the name of ggdat
                  for (cluster_k in 1:k) {
                    # cluster_k <- 3
                    new_label <- paste0("cluster", 
                                        rank(pk_stat)[cluster_k])
                    
                    cluster$cluster[cluster$cluster == cluster_k] <- new_label
                  }
                  
                  # format the cluster name
                  cluster$cluster <- as.numeric(gsub("cluster", "", cluster$cluster))
                  
                  rm(new_label, cluster_k)
                }
                
                # rename pol2 pat pk location
                if (T) {
                  names(pk_stat) <- paste0("pol2_", rank(pk_stat))
                  pk_stat <- sort(pk_stat)
                }
                
                res <- list(pk_stat, cluster)
                
                qsave(res, file = file)
              }
              if (file.exists(file)) {
                res <- qread(file, nthreads = 6)
                cluster <- res[[2]]
              }
            }
            
            # plot pol2 pattern distribution
            if (T) {
              # merge and scale the signal
              if (T) {
                profile <- data.frame(t(scale(t(submat[, name]))))
                profile$cluster <- cluster$cluster
              }
              
              # get average signal
              if (T) {
                ggdat <- average_sig(profile, average_type="mean")
              }
              
              # plot overlap pol2 pat
              if (T) {
                head(ggdat)
                
                p_base <- ggplot(ggdat) +
                  geom_segment(aes(x=Loc, y=min(Signal)*0.5, xend=Loc, yend=Assistant), linetype=2, linewidth=0.8, color="grey", na.rm=T) +
                  geom_line(aes(x=Loc, y=Signal, group=cluster, color=cluster), linewidth=1) +
                  scale_x_discrete(name = NULL, breaks = labels) +
                  ylab("Pol2 Enrichment Level") +
                  ggtitle(paste0(cell, ": K", k, " ", type, " + pattern")) +
                  cowplot::theme_cowplot() +
                  theme(panel.border = element_rect(color="black"), 
                        strip.background = element_rect(fill=NA, color=NA), 
                        legend.position = "none")
                
                file <- paste0("50.pol2_pat_kmeans", k, "_overlap_", data_type, "s_profile_pattern.", cell, ".", type, ".pdf")
                ggsave(filename = file, plot = p_base, width = 8, height = 6)
                move_file(file, dir2)
              }
              
              # plot sep pol2 pat
              if (T) {
                p <- p_base + facet_grid(.~cluster)
                
                file <- paste0("50.pol2_pat_kmeans", k, "_sep_", data_type, "s_profile_pattern.", cell, ".", type, ".pdf")
                ggsave(filename = file, plot = p, width = 12, height = 6)
                move_file(file, dir2)
              }
              
              # plot all gene pol2 pat distribution
              if (T) {
                # get all gene pol2 distribution formatted data
                if (T) {
                  # long2wide
                  profile$gene_id <- rownames(profile)
                  ggdat <- melt(profile, id.vars = c("gene_id", "cluster"), measure.vars = name, variable.name = "Loc", value.name = "Signal")
                  
                  # format the data
                  ggdat$cluster <- factor(paste0("pol2_", ggdat$cluster), levels = paste0("pol2_", 1:k))
                  ggdat$Loc <- factor(ggdat$Loc, levels = name)
                }
                
                head(ggdat)
                
                # ggplot2
                if (T) {
                  p <- ggplot(ggdat) +
                    geom_line(aes(x=Loc, y=Signal, group=gene_id, color=cluster), linewidth=1, alpha=0.05) +
                    scale_x_discrete(name = NULL, breaks = labels) +
                    ylab("Pol2 Enrichment Level") +
                    ggtitle(paste0(cell, ": ", type, " + pattern")) +
                    facet_grid(.~cluster) +
                    cowplot::theme_cowplot() +
                    theme(panel.border = element_rect(color="black"), 
                          strip.background = element_rect(fill=NA, color=NA), 
                          legend.position = "none")
                  
                  file <- paste0("50.pol2_pat_kmeans", k, "_all_", data_type, "s_profile_pattern.", cell, ".", type, ".jpeg")
                  ggsave(filename = file, plot = p, width = 12, height = 6)
                  
                  move_file(file, dir2)
                }
              }
            }
            
            # statistics the number of genes for each group
            if (T) {
              # statistics the number
              if (T) {
                group <- as.data.frame.table(table(paste0("pol2_", profile$cluster)))
                head(group)
                colnames(group) <- c("cluster", "Freq")
              }
              
              # calculate the percentage of each cluster
              if (T) {
                group$Percentage <- group$Freq/sum(group$Freq)*100
              }
              
              # format the info
              if (T) {
                group$label <- paste0(round(group$Percentage, 2), "%")
                group$sci <- format(as.numeric(group$Freq), scientific = T, digits = 3)
                
                group$cluster <- factor(group$cluster, levels = sort(unique(group$cluster)))
              }
              
              # ggplot for the number of each group
              if (T) {
                # percentage ggplot
                if (T) {
                  head(group)
                  
                  p <- ggplot(group) +
                    geom_hline(yintercept = seq(0, ceiling(max(group$Percentage)/10)*10, 5), alpha=0.4, colour="darkgrey") +
                    geom_bar(aes(x=cluster, y=Percentage, fill=cluster), stat = "identity") +
                    geom_text(aes(x=cluster, y=Percentage+0.5, label=label)) +
                    scale_x_discrete(name = NULL) +
                    scale_y_continuous(limits = c(0, ceiling(max(group$Percentage)/10)*10), breaks = seq(0, ceiling(max(group$Percentage)/10)*10, 5)) +
                    ylab("Gene Percentage") +
                    ggtitle(paste0(cell, ": ", type, " + pattern")) +
                    cowplot::theme_cowplot() +
                    theme(panel.border = element_rect(color="black"), 
                          strip.background = element_rect(fill=NA, color=NA), 
                          legend.position = "none")
                  
                  file <- paste0("50.pol2_pat_kmeans", k, "_", data_type, "s_profile_pat_percentage.", cell, ".", type, ".pdf")
                  ggsave(filename = file, plot = p, width = 8, height = 6)
                  
                  move_file(file, dir2)
                }
                
                # absolute number ggplot
                if (T) {
                  head(group)
                  
                  p <- ggplot(group) +
                    geom_bar(aes(x=cluster, y=Freq, fill=cluster), stat = "identity") +
                    geom_text(aes(x=cluster, y=Freq+300, label=Freq)) +
                    scale_x_discrete(name = NULL) +
                    ylab("Gene Absolute Number") +
                    ggtitle(paste0(cell, ": ", type, " + pattern")) +
                    cowplot::theme_cowplot() +
                    theme(panel.border = element_rect(color="black"), 
                          strip.background = element_rect(fill=NA, color=NA), 
                          legend.position = "none")
                  
                  file <- paste0("50.pol2_pat_kmeans", k, "_", data_type, "s_profile_pat_num.", cell, ".", type, ".pdf")
                  ggsave(filename = file, plot = p, width = 8, height = 6)
                  
                  move_file(file, dir2)
                }
              }
            }
          }
        }
      }
      
      # merged data
      if (type == "merged") {
        cat(paste0("\t\tmerged: \n"))
        
        # test data
        if (F) {
          data_type <- "tx"
          type <- "merged"
        }
        
        # read kmeans result
        if (T) {
          dat <- qread(paste0(dir1, "/50.pol2_pat_bs", bin_size, "_kmeans1-15_", data_type, "_merged.", type, ".qs"), nthreads = 6)
        }
        
        # get pol2 matrix: mat
        if (T) {
          cell_dat <- function(cell, type) {
            input <- paste0("data/saved_data/48.mk_signal_distribution/tx/48.", data_type, "_Body_", up_bin_num, "UP_", 
                            down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".pol2.mean_mat.", cell, ".", type, ".qs")
            
            mat <- qread(input, nthreads = 6)
            rownames(mat) <- paste0(mat$gene_id, "@", cell)
            
            return(mat)
          }
          dat1 <- cell_dat(cell1, "single")
          dat2 <- cell_dat(cell2, "single")
          
          mat <- rbind(dat1, dat2)
          
          rm(dat1, dat2)
        }
        
        # kmeans: k2-6
        for (k in 2:6) {
          # k <- 5
          cat(paste0("\t\tk=", k, "\n"))
          
          # get cluster group
          if (T) {
            # get kmeans cluster data: cluster
            if (T) {
              cluster <- dat[[paste0("K", k)]]
              cluster <- cluster$cluster
              cluster <- data.frame(ID = names(cluster), 
                                    cluster = as.numeric(cluster))
            }
            
            # get consistent mat: submat
            if (T) {
              submat <- mat[cluster$ID, ]
            }
            
            file <- paste0(dir1, "/50.pol2_pat_bs", bin_size, "_kmeans", k, "_", data_type, "_merged.", type, ".qs")
            if (! file.exists(file)) {
              # get pol2 pat pk location
              if (T) {
                submat_row <- data.frame(t(scale(t(submat[, name]))))
                pk_stat <- pattern_group_pk_location(submat_row, cluster, average_type="mean")
              }
              
              # rename the cluster label: assign the cluster label based on the pol2 peak location
              if (T) {
                # change the name of ggdat
                for (cluster_k in 1:k) {
                  # cluster_k <- 3
                  new_label <- paste0("cluster", 
                                      rank(pk_stat)[cluster_k])
                  
                  cluster$cluster[cluster$cluster == cluster_k] <- new_label
                }
                
                # format the cluster name
                cluster$cluster <- as.numeric(gsub("cluster", "", cluster$cluster))
                
                rm(new_label, cluster_k)
              }
              
              # rename pol2 pat pk location
              if (T) {
                names(pk_stat) <- paste0("pol2_", rank(pk_stat))
                pk_stat <- sort(pk_stat)
              }
              
              res <- list(pk_stat, cluster)
              
              qsave(res, file = file)
            }
            if (file.exists(file)) {
              res <- qread(file, nthreads = 6)
              cluster <- res[[2]]
            }
          }
          
          # plot pol2 pattern distribution
          if (T) {
            # merge and scale the signal
            if (T) {
              profile <- data.frame(t(scale(t(submat[, name]))))
              profile$cluster <- cluster$cluster
            }
            
            # get average signal
            if (T) {
              ggdat <- average_sig(profile, average_type="mean")
            }
            
            # plot overlap pol2 pat
            if (T) {
              head(ggdat)
              
              p_base <- ggplot(ggdat) +
                geom_segment(aes(x=Loc, y=min(Signal)*0.5, xend=Loc, yend=Assistant), linetype=2, linewidth=0.8, color="grey", na.rm=T) +
                geom_line(aes(x=Loc, y=Signal, group=cluster, color=cluster), linewidth=1) +
                scale_x_discrete(name = NULL, breaks = labels) +
                ylab("Pol2 Enrichment Level") +
                ggtitle(paste0("merged: K", k, " ", type, " + pattern")) +
                cowplot::theme_cowplot() +
                theme(panel.border = element_rect(color="black"), 
                      strip.background = element_rect(fill=NA, color=NA), 
                      legend.position = "none")
              
              file <- paste0("50.pol2_pat_kmeans", k, "_overlap_", data_type, "s_profile_pattern.merged.", type, ".pdf")
              ggsave(filename = file, plot = p_base, width = 8, height = 6)
              move_file(file, dir2)
            }
            
            # plot sep pol2 pat
            if (T) {
              p <- p_base + facet_grid(.~cluster)
              
              file <- paste0("50.pol2_pat_kmeans", k, "_sep_", data_type, "s_profile_pattern.merged.", type, ".pdf")
              ggsave(filename = file, plot = p, width = 12, height = 6)
              move_file(file, dir2)
            }
            
            # plot all gene pol2 pat distribution
            if (T) {
              # get all gene pol2 distribution formatted data
              if (T) {
                # long2wide
                profile$gene_id <- rownames(profile)
                ggdat <- melt(profile, id.vars = c("gene_id", "cluster"), measure.vars = name, variable.name = "Loc", value.name = "Signal")
                
                # format the data
                ggdat$cluster <- factor(paste0("pol2_", ggdat$cluster), levels = paste0("pol2_", 1:k))
                ggdat$Loc <- factor(ggdat$Loc, levels = name)
              }
              
              head(ggdat)
              
              # ggplot2
              if (T) {
                p <- ggplot(ggdat) +
                  geom_line(aes(x=Loc, y=Signal, group=gene_id, color=cluster), linewidth=1, alpha=0.05) +
                  scale_x_discrete(name = NULL, breaks = labels) +
                  ylab("Pol2 Enrichment Level") +
                  ggtitle(paste0("merged: ", type, " + pattern")) +
                  facet_grid(.~cluster) +
                  cowplot::theme_cowplot() +
                  theme(panel.border = element_rect(color="black"), 
                        strip.background = element_rect(fill=NA, color=NA), 
                        legend.position = "none")
                
                file <- paste0("50.pol2_pat_kmeans", k, "_all_", data_type, "s_profile_pattern.merged.", type, ".jpeg")
                ggsave(filename = file, plot = p, width = 12, height = 6)
                
                move_file(file, dir2)
              }
            }
          }
          
          # statistics the number of genes for each group
          if (T) {
            # statistics the number
            if (T) {
              group <- as.data.frame.table(table(paste0("pol2_", profile$cluster)))
              head(group)
              colnames(group) <- c("cluster", "Freq")
            }
            
            # calculate the percentage of each cluster
            if (T) {
              group$Percentage <- group$Freq/sum(group$Freq)*100
            }
            
            # format the info
            if (T) {
              group$label <- paste0(round(group$Percentage, 2), "%")
              group$sci <- format(as.numeric(group$Freq), scientific = T, digits = 3)
              
              group$cluster <- factor(group$cluster, levels = sort(unique(group$cluster)))
            }
            
            # ggplot for the number of each group
            if (T) {
              # percentage ggplot
              if (T) {
                head(group)
                
                p <- ggplot(group) +
                  geom_hline(yintercept = seq(0, ceiling(max(group$Percentage)/10)*10, 5), alpha=0.4, colour="darkgrey") +
                  geom_bar(aes(x=cluster, y=Percentage, fill=cluster), stat = "identity") +
                  geom_text(aes(x=cluster, y=Percentage+0.5, label=label)) +
                  scale_x_discrete(name = NULL) +
                  scale_y_continuous(limits = c(0, ceiling(max(group$Percentage)/10)*10), breaks = seq(0, ceiling(max(group$Percentage)/10)*10, 5)) +
                  ylab("Gene Percentage") +
                  ggtitle(paste0("merged: ", type, " + pattern")) +
                  cowplot::theme_cowplot() +
                  theme(panel.border = element_rect(color="black"), 
                        strip.background = element_rect(fill=NA, color=NA), 
                        legend.position = "none")
                
                file <- paste0("50.pol2_pat_kmeans", k, "_", data_type, "s_profile_pat_percentage.merged.", type, ".pdf")
                ggsave(filename = file, plot = p, width = 8, height = 6)
                
                move_file(file, dir2)
              }
              
              # absolute number ggplot
              if (T) {
                head(group)
                
                p <- ggplot(group) +
                  geom_bar(aes(x=cluster, y=Freq, fill=cluster), stat = "identity") +
                  geom_text(aes(x=cluster, y=Freq+300, label=Freq)) +
                  scale_x_discrete(name = NULL) +
                  ylab("Gene Absolute Number") +
                  ggtitle(paste0("merged: ", type, " + pattern")) +
                  cowplot::theme_cowplot() +
                  theme(panel.border = element_rect(color="black"), 
                        strip.background = element_rect(fill=NA, color=NA), 
                        legend.position = "none")
                
                file <- paste0("50.pol2_pat_kmeans", k, "_", data_type, "s_profile_pat_num.merged.", type, ".pdf")
                ggsave(filename = file, plot = p, width = 8, height = 6)
                
                move_file(file, dir2)
              }
            }
          }
        }
      }
    }
  }
}

# get cluster when K in 2:6 (pat col)
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir1 <- paste0("data/saved_data/50.pol2_pattern/", data_type)
      if (! dir.exists(dir1)) {
        dir.create(dir1, showWarnings = F, recursive = T)
      }
      
      dir2 <- paste0("results/2.pic/50.pol2_pattern/", data_type)
      if (! dir.exists(dir2)) {
        dir.create(dir2, showWarnings = F, recursive = T)
      }
    }
    
    for (type in c("merged", "single")) {
      cat(paste0("\t", type, ": \n"))
      
      # cell specific data
      if (type == "single") {
        for (cell in c(cell1, cell2)) {
          cat(paste0("\t\t", cell, ": \n"))
          
          # test data
          if (F) {
            data_type <- "tx"
            type <- "single"
            cell <- cell1
          }
          
          # read kmeans result
          if (T) {
            dat <- qread(paste0(dir1, "/50.pol2_sig_bs", bin_size, "_kmeans1-15_", data_type, "_", cell, ".", type, ".qs"), nthreads = 6)
          }
          
          # get pol2 matrix: mat
          if (T) {
            cell_dat <- function(cell, type) {
              input <- paste0("data/saved_data/48.mk_signal_distribution/tx/48.", data_type, "_Body_", up_bin_num, "UP_", 
                              down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".pol2.mean_mat.", cell, ".", type, ".qs")
              
              mat <- qread(input, nthreads = 6)
              rownames(mat) <- paste0(mat$gene_id, "@", cell)
              
              return(mat)
            }
            mat <- cell_dat(cell, type)
          }
          
          # kmeans: k2-6
          for (k in 2:6) {
            # k <- 5
            cat(paste0("\t\tk=", k, "\n"))
            
            # get cluster group
            if (T) {
              # get kmeans cluster data: cluster
              if (T) {
                cluster <- dat[[paste0("K", k)]]
                cluster <- cluster$cluster
                cluster <- data.frame(ID = names(cluster), 
                                      cluster = as.numeric(cluster))
              }
              
              # get consistent mat: submat
              if (T) {
                submat <- mat[cluster$ID, ]
              }
              
              file <- paste0(dir1, "/50.pol2_sig_bs", bin_size, "_kmeans", k, "_", data_type, "_", cell, ".", type, ".qs")
              if (! file.exists(file)) {
                # get pol2 pat pk location
                if (T) {
                  submat_col <- data.frame(scale(submat[, name]))
                  pk_stat <- pattern_group_pk_location(submat_col, cluster, average_type="mean")
                }
                
                # rename the cluster label: assign the cluster label based on the pol2 peak location
                if (T) {
                  # change the name of ggdat
                  for (cluster_k in 1:k) {
                    # cluster_k <- 3
                    new_label <- paste0("cluster", 
                                        rank(pk_stat)[cluster_k])
                    
                    cluster$cluster[cluster$cluster == cluster_k] <- new_label
                  }
                  
                  # format the cluster name
                  cluster$cluster <- as.numeric(gsub("cluster", "", cluster$cluster))
                  
                  rm(new_label, cluster_k)
                }
                
                # rename pol2 pat pk location
                if (T) {
                  names(pk_stat) <- paste0("pol2_", rank(pk_stat))
                  pk_stat <- sort(pk_stat)
                }
                
                res <- list(pk_stat, cluster)
                
                qsave(res, file = file)
              }
              if (file.exists(file)) {
                res <- qread(file, nthreads = 6)
                cluster <- res[[2]]
              }
            }
            
            # plot pol2 signal distribution
            if (T) {
              # merge and scale the signal
              if (T) {
                profile <- data.frame(scale(submat[, name]))
                profile$cluster <- cluster$cluster
              }
              
              # get average signal
              if (T) {
                ggdat <- average_sig(profile, average_type="mean")
              }
              
              # plot overlap pol2 sig
              if (T) {
                head(ggdat)
                
                p_base <- ggplot(ggdat) +
                  geom_segment(aes(x=Loc, y=min(Signal)*0.5, xend=Loc, yend=Assistant), linetype=2, linewidth=0.8, color="grey", na.rm=T) +
                  geom_line(aes(x=Loc, y=Signal, group=cluster, color=cluster), linewidth=1) +
                  scale_x_discrete(name = NULL, breaks = labels) +
                  ylab("Pol2 Enrichment Level") +
                  ggtitle(paste0(cell, ": K", k, " ", type, " + pattern")) +
                  cowplot::theme_cowplot() +
                  theme(panel.border = element_rect(color="black"), 
                        strip.background = element_rect(fill=NA, color=NA), 
                        legend.position = "none")
                
                file <- paste0("50.pol2_sig_kmeans", k, "_overlap_", data_type, "s_profile_signal.", cell, ".", type, ".pdf")
                ggsave(filename = file, plot = p_base, width = 8, height = 6)
                move_file(file, dir2)
              }
              
              # plot sep pol2 sig
              if (T) {
                p <- p_base + facet_grid(.~cluster)
                
                file <- paste0("50.pol2_sig_kmeans", k, "_sep_", data_type, "s_profile_signal.", cell, ".", type, ".pdf")
                ggsave(filename = file, plot = p, width = 12, height = 6)
                move_file(file, dir2)
              }
              
              # plot all gene pol2 sig distribution
              if (T) {
                # get all gene pol2 distribution formatted data
                if (T) {
                  # long2wide
                  profile$gene_id <- rownames(profile)
                  ggdat <- melt(profile, id.vars = c("gene_id", "cluster"), measure.vars = name, variable.name = "Loc", value.name = "Signal")
                  
                  # format the data
                  ggdat$cluster <- factor(paste0("pol2_", ggdat$cluster), levels = paste0("pol2_", 1:k))
                  ggdat$Loc <- factor(ggdat$Loc, levels = name)
                }
                
                head(ggdat)
                
                # ggplot2
                if (T) {
                  p <- ggplot(ggdat) +
                    geom_line(aes(x=Loc, y=Signal, group=gene_id, color=cluster), linewidth=1, alpha=0.05) +
                    scale_x_discrete(name = NULL, breaks = labels) +
                    ylab("Pol2 Enrichment Level") +
                    ggtitle(paste0(cell, ": ", type, " + pattern")) +
                    facet_grid(.~cluster) +
                    cowplot::theme_cowplot() +
                    theme(panel.border = element_rect(color="black"), 
                          strip.background = element_rect(fill=NA, color=NA), 
                          legend.position = "none")
                  
                  file <- paste0("50.pol2_sig_kmeans", k, "_all_", data_type, "s_profile_signal.", cell, ".", type, ".jpeg")
                  ggsave(filename = file, plot = p, width = 12, height = 6)
                  
                  move_file(file, dir2)
                }
              }
            }
            
            # statistics the number of genes for each group
            if (T) {
              # statistics the number
              if (T) {
                group <- as.data.frame.table(table(paste0("pol2_", profile$cluster)))
                head(group)
                colnames(group) <- c("cluster", "Freq")
              }
              
              # calculate the percentage of each cluster
              if (T) {
                group$Percentage <- group$Freq/sum(group$Freq)*100
              }
              
              # format the info
              if (T) {
                group$label <- paste0(round(group$Percentage, 2), "%")
                group$sci <- format(as.numeric(group$Freq), scientific = T, digits = 3)
                
                group$cluster <- factor(group$cluster, levels = sort(unique(group$cluster)))
              }
              
              # ggplot for the number of each group
              if (T) {
                # percentage ggplot
                if (T) {
                  head(group)
                  
                  p <- ggplot(group) +
                    geom_hline(yintercept = seq(0, ceiling(max(group$Percentage)/10)*10, 5), alpha=0.4, colour="darkgrey") +
                    geom_bar(aes(x=cluster, y=Percentage, fill=cluster), stat = "identity") +
                    geom_text(aes(x=cluster, y=Percentage+0.5, label=label)) +
                    scale_x_discrete(name = NULL) +
                    scale_y_continuous(limits = c(0, ceiling(max(group$Percentage)/10)*10), breaks = seq(0, ceiling(max(group$Percentage)/10)*10, 5)) +
                    ylab("Gene Percentage") +
                    ggtitle(paste0(cell, ": ", type, " + pattern")) +
                    cowplot::theme_cowplot() +
                    theme(panel.border = element_rect(color="black"), 
                          strip.background = element_rect(fill=NA, color=NA), 
                          legend.position = "none")
                  
                  file <- paste0("50.pol2_sig_kmeans", k, "_", data_type, "s_profile_sig_percentage.", cell, ".", type, ".pdf")
                  ggsave(filename = file, plot = p, width = 8, height = 6)
                  
                  move_file(file, dir2)
                }
                
                # absolute number ggplot
                if (T) {
                  head(group)
                  
                  p <- ggplot(group) +
                    geom_bar(aes(x=cluster, y=Freq, fill=cluster), stat = "identity") +
                    geom_text(aes(x=cluster, y=Freq+300, label=Freq)) +
                    scale_x_discrete(name = NULL) +
                    ylab("Gene Absolute Number") +
                    ggtitle(paste0(cell, ": ", type, " + pattern")) +
                    cowplot::theme_cowplot() +
                    theme(panel.border = element_rect(color="black"), 
                          strip.background = element_rect(fill=NA, color=NA), 
                          legend.position = "none")
                  
                  file <- paste0("50.pol2_sig_kmeans", k, "_", data_type, "s_profile_sig_num.", cell, ".", type, ".pdf")
                  ggsave(filename = file, plot = p, width = 8, height = 6)
                  
                  move_file(file, dir2)
                }
              }
            }
          }
        }
      }
      
      # merged data
      if (type == "merged") {
        cat(paste0("\t\tmerged: \n"))
        
        # test data
        if (F) {
          data_type <- "tx"
          type <- "merged"
        }
        
        # read kmeans result
        if (T) {
          dat <- qread(paste0(dir1, "/50.pol2_sig_bs", bin_size, "_kmeans1-15_", data_type, "_merged.", type, ".qs"), nthreads = 6)
        }
        
        # get pol2 matrix: mat
        if (T) {
          cell_dat <- function(cell, type) {
            input <- paste0("data/saved_data/48.mk_signal_distribution/tx/48.", data_type, "_Body_", up_bin_num, "UP_", 
                            down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".pol2.mean_mat.", cell, ".", type, ".qs")
            
            mat <- qread(input, nthreads = 6)
            rownames(mat) <- paste0(mat$gene_id, "@", cell)
            
            return(mat)
          }
          dat1 <- cell_dat(cell1, "single")
          dat2 <- cell_dat(cell2, "single")
          
          mat <- rbind(dat1, dat2)
          
          rm(dat1, dat2)
        }
        
        # kmeans: k2-6
        for (k in 2:6) {
          # k <- 5
          cat(paste0("\t\tk=", k, "\n"))
          
          # get cluster group
          if (T) {
            # get kmeans cluster data: cluster
            if (T) {
              cluster <- dat[[paste0("K", k)]]
              cluster <- cluster$cluster
              cluster <- data.frame(ID = names(cluster), 
                                    cluster = as.numeric(cluster))
            }
            
            # get consistent mat: submat
            if (T) {
              submat <- mat[cluster$ID, ]
            }
            
            file <- paste0(dir1, "/50.pol2_sig_bs", bin_size, "_kmeans", k, "_", data_type, "_merged.", type, ".qs")
            if (! file.exists(file)) {
              # get pol2 pat pk location
              if (T) {
                submat_col <- data.frame(scale(submat[, name]))
                pk_stat <- pattern_group_pk_location(submat_col, cluster, average_type="mean")
              }
              
              # rename the cluster label: assign the cluster label based on the pol2 peak location
              if (T) {
                # change the name of ggdat
                for (cluster_k in 1:k) {
                  # cluster_k <- 3
                  new_label <- paste0("cluster", 
                                      rank(pk_stat)[cluster_k])
                  
                  cluster$cluster[cluster$cluster == cluster_k] <- new_label
                }
                
                # format the cluster name
                cluster$cluster <- as.numeric(gsub("cluster", "", cluster$cluster))
                
                rm(new_label, cluster_k)
              }
              
              # rename pol2 pat pk location
              if (T) {
                names(pk_stat) <- paste0("pol2_", rank(pk_stat))
                pk_stat <- sort(pk_stat)
              }
              
              res <- list(pk_stat, cluster)
              
              qsave(res, file = file)
            }
            if (file.exists(file)) {
              res <- qread(file, nthreads = 6)
              cluster <- res[[2]]
            }
          }
          
          # plot pol2 signal distribution
          if (T) {
            # merge and scale the signal
            if (T) {
              profile <- data.frame(scale(submat[, name]))
              profile$cluster <- cluster$cluster
            }
            
            # get average signal
            if (T) {
              ggdat <- average_sig(profile, average_type="mean")
            }
            
            # plot overlap pol2 sig
            if (T) {
              head(ggdat)
              
              p_base <- ggplot(ggdat) +
                geom_segment(aes(x=Loc, y=min(Signal)*0.5, xend=Loc, yend=Assistant), linetype=2, linewidth=0.8, color="grey", na.rm=T) +
                geom_line(aes(x=Loc, y=Signal, group=cluster, color=cluster), linewidth=1) +
                scale_x_discrete(name = NULL, breaks = labels) +
                ylab("Pol2 Enrichment Level") +
                ggtitle(paste0("merged: K", k, " ", type, " + pattern")) +
                cowplot::theme_cowplot() +
                theme(panel.border = element_rect(color="black"), 
                      strip.background = element_rect(fill=NA, color=NA), 
                      legend.position = "none")
              
              file <- paste0("50.pol2_sig_kmeans", k, "_overlap_", data_type, "s_profile_signal.merged.", type, ".pdf")
              ggsave(filename = file, plot = p_base, width = 8, height = 6)
              move_file(file, dir2)
            }
            
            # plot sep pol2 sig
            if (T) {
              p <- p_base + facet_grid(.~cluster)
              
              file <- paste0("50.pol2_sig_kmeans", k, "_sep_", data_type, "s_profile_signal.merged.", type, ".pdf")
              ggsave(filename = file, plot = p, width = 12, height = 6)
              move_file(file, dir2)
            }
            
            # plot all gene pol2 sig distribution
            if (T) {
              # get all gene pol2 distribution formatted data
              if (T) {
                # long2wide
                profile$gene_id <- rownames(profile)
                ggdat <- melt(profile, id.vars = c("gene_id", "cluster"), measure.vars = name, variable.name = "Loc", value.name = "Signal")
                
                # format the data
                ggdat$cluster <- factor(paste0("pol2_", ggdat$cluster), levels = paste0("pol2_", 1:k))
                ggdat$Loc <- factor(ggdat$Loc, levels = name)
              }
              
              head(ggdat)
              
              # ggplot2
              if (T) {
                p <- ggplot(ggdat) +
                  geom_line(aes(x=Loc, y=Signal, group=gene_id, color=cluster), linewidth=1, alpha=0.05) +
                  scale_x_discrete(name = NULL, breaks = labels) +
                  ylab("Pol2 Enrichment Level") +
                  ggtitle(paste0("merged: ", type, " + pattern")) +
                  facet_grid(.~cluster) +
                  cowplot::theme_cowplot() +
                  theme(panel.border = element_rect(color="black"), 
                        strip.background = element_rect(fill=NA, color=NA), 
                        legend.position = "none")
                
                file <- paste0("50.pol2_sig_kmeans", k, "_all_", data_type, "s_profile_signal.merged.", type, ".jpeg")
                ggsave(filename = file, plot = p, width = 12, height = 6)
                
                move_file(file, dir2)
              }
            }
          }
          
          # statistics the number of genes for each group
          if (T) {
            # statistics the number
            if (T) {
              group <- as.data.frame.table(table(paste0("pol2_", profile$cluster)))
              head(group)
              colnames(group) <- c("cluster", "Freq")
            }
            
            # calculate the percentage of each cluster
            if (T) {
              group$Percentage <- group$Freq/sum(group$Freq)*100
            }
            
            # format the info
            if (T) {
              group$label <- paste0(round(group$Percentage, 2), "%")
              group$sci <- format(as.numeric(group$Freq), scientific = T, digits = 3)
              
              group$cluster <- factor(group$cluster, levels = sort(unique(group$cluster)))
            }
            
            # ggplot for the number of each group
            if (T) {
              # percentage ggplot
              if (T) {
                head(group)
                
                p <- ggplot(group) +
                  geom_hline(yintercept = seq(0, ceiling(max(group$Percentage)/10)*10, 5), alpha=0.4, colour="darkgrey") +
                  geom_bar(aes(x=cluster, y=Percentage, fill=cluster), stat = "identity") +
                  geom_text(aes(x=cluster, y=Percentage+0.5, label=label)) +
                  scale_x_discrete(name = NULL) +
                  scale_y_continuous(limits = c(0, ceiling(max(group$Percentage)/10)*10), breaks = seq(0, ceiling(max(group$Percentage)/10)*10, 5)) +
                  ylab("Gene Percentage") +
                  ggtitle(paste0("merged: ", type, " + pattern")) +
                  cowplot::theme_cowplot() +
                  theme(panel.border = element_rect(color="black"), 
                        strip.background = element_rect(fill=NA, color=NA), 
                        legend.position = "none")
                
                file <- paste0("50.pol2_sig_kmeans", k, "_", data_type, "s_profile_sig_percentage.merged.", type, ".pdf")
                ggsave(filename = file, plot = p, width = 8, height = 6)
                
                move_file(file, dir2)
              }
              
              # absolute number ggplot
              if (T) {
                head(group)
                
                p <- ggplot(group) +
                  geom_bar(aes(x=cluster, y=Freq, fill=cluster), stat = "identity") +
                  geom_text(aes(x=cluster, y=Freq+300, label=Freq)) +
                  scale_x_discrete(name = NULL) +
                  ylab("Gene Absolute Number") +
                  ggtitle(paste0("merged: ", type, " + pattern")) +
                  cowplot::theme_cowplot() +
                  theme(panel.border = element_rect(color="black"), 
                        strip.background = element_rect(fill=NA, color=NA), 
                        legend.position = "none")
                
                file <- paste0("50.pol2_sig_kmeans", k, "_", data_type, "s_profile_sig_num.merged.", type, ".pdf")
                ggsave(filename = file, plot = p, width = 8, height = 6)
                
                move_file(file, dir2)
              }
            }
          }
        }
      }
    }
  }
}
