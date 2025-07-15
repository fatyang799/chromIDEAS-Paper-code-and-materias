# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(stringr)
  library(reshape2)
  library(qs)
  suppressPackageStartupMessages(library(circlize))
  suppressPackageStartupMessages(library(ComplexHeatmap))
}

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
up_bin_num <- 5
down_bin_num <- 5
body_bin_num <- 10
pol2_levels <- 10

# define functions
if (T) {
  # correlation (modify based on script39)
  plot_state_auc_cor_heatmap <- function(cor_dat, cluster) {
    # get cell specific cluster
    if (T) {
      group <- cluster[match(rownames(cor_dat), cluster$orig.ident), "seurat_clusters"]
      group <- paste0("C", group)
    }
    
    # color setting
    if (T) {
      library(RColorBrewer)
      colors <- brewer.pal(9, "Set1") 
      
      colors <- colors[1:2]
      col_fun <- colorRamp2(c(-1, 0, 1), c(colors[2], "white", colors[1]))
    }
    
    # heatmap
    if (T) {
      p <- Heatmap(cor_dat, name = "Pearson", col = col_fun, 
                   border_gp = gpar(col = "black"), rect_gp = gpar(col = "black"),
                   show_heatmap_legend = TRUE, heatmap_legend_param = list(border = T),
                   row_names_max_width = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), row_names_rot = 0,
                   column_names_max_height = unit(6, "cm"), column_names_gp = gpar(fontsize = 10), column_names_rot = 90,
                   cluster_rows = T, cluster_columns = T, show_row_dend = T, show_column_dend = F, 
                   clustering_distance_columns = "euclidean", clustering_distance_rows = "euclidean", 
                   clustering_method_columns = "ward.D2", clustering_method_rows = "ward.D2", 
                   right_annotation = rowAnnotation(Group = group, border = T, 
                                                    col = list(Group = c("C1" = "#e76f51", 
                                                                         "C2" = "#8ecae6", 
                                                                         "C3" = "#f4a261", 
                                                                         "C4" = "#219ebc", 
                                                                         "C5" = "#a7c957")), 
                                                    annotation_legend_param = list(border = "black")))
    }
    
    p <- ggplotify::as.ggplot(p)
    return(p)
  }
  
  # k-means (script39)
  kmeans_determine_k <- function(mat, title) {
    # test dat
    if (F) {
      title <- "test"
    }
    
    # scale data
    if (T) {
      df <- mat
    }
    
    # determine the k value: Within groups sum of squares
    if (T) {
      wssplot <- function(df, nc=15, title, seed=799){
        wss <- c()
        for (i in 1:nc){
          set.seed(seed)
          wss[i] <- sum(kmeans(df, centers=i, nstart=25)$withinss)
        }
        ggdat <- data.frame(k=1:nc, y=wss)
        
        p <- ggplot(ggdat) +
          geom_point(aes(x=k, y=wss), size=3) +
          geom_line(aes(x=k, y=wss, group=1)) +
          scale_x_continuous(breaks = 1:nc) +
          xlab("Number of Clusters") +
          ylab("Within groups sum of squares") +
          ggtitle(title) +
          theme_bw() +
          theme(axis.title = element_text(size = rel(1.2)),
                axis.text = element_text(size = rel(1.2)),
                legend.text = element_text(size = rel(1.2)),
                legend.title = element_text(size = rel(1.2)),
                strip.text = element_text(size = rel(1.2)))
        return(p)
      }
      p <- wssplot(df, nc=15, title, seed=799)
    }
    
    return(p)
  }
  
  # k-means (script50)
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
  
  # line plot clustered by chromIDEAS (script38)
  plot_line_add_merged_line <- function(mat, cluster) {
    # wide2long
    if (T) {
      mat1 <- data.frame(t(scale(t(mat))))
      mat1$State <- rownames(mat1)
      mat1 <- melt(mat1, id.vars = "State", variable.name = "pol2_sig")
      
      mat1$State <- factor(mat1$State, levels = paste0("S", all_states))
      mat1$pol2_sig <- factor(mat1$pol2_sig, levels = all_pol2_sig)
    }
    
    # add clustering data
    if (T) {
      mat1$cluster <- cluster$seurat_clusters[match(mat1$State, cluster$orig.ident)]
      mat1$cluster <- factor(mat1$cluster, levels = sort(unique(mat1$cluster)))
    }
    
    # plot the line plot
    if (T) {
      head(mat1)
      
      p <- ggplot(mat1) +
        geom_line(aes(x=pol2_sig, y=value, group=State, color=State), linewidth=1, alpha=0.3) +
        xlab(NULL) +
        ylab("Relative Preference") +
        facet_grid(.~cluster) +
        cowplot::theme_cowplot() +
        theme(legend.position = "none", 
              panel.border = element_rect(color="black"), 
              strip.background = element_rect(fill=NA, color=NA))
    }
    
    # get merged line
    if (T) {
      merged_line <- tapply(mat1$value, paste0(mat1$cluster, "@", mat1$pol2_sig), median)
      merged_line <- data.frame(id = names(merged_line), 
                                value = as.numeric(merged_line))
      merged_line$cluster <- str_split(merged_line$id, "@", simplify = T)[, 1]
      merged_line$pol2_sig <- str_split(merged_line$id, "@", simplify = T)[, 2]
      
      merged_line$cluster <- factor(merged_line$cluster, levels = levels(mat1$cluster))
      merged_line$pol2_sig <- factor(merged_line$pol2_sig, levels = levels(mat1$pol2_sig))
    }
    
    # add merged line into the ggplot
    if (T) {
      head(merged_line)
      merged_p <- p + ggnewscale::new_scale_color() +
        geom_line(data = merged_line, aes(x=pol2_sig, y=value, group=cluster, color=cluster), linewidth=2)
    }
    
    return(merged_p)
  }
  
  # kmeans plot (modify based on script39)
  plot_state_auc_distribution_kmeans <- function(mat, k, title) {
    # test dat
    if (F) {
      title <- "test"
      k <- 3
    }
    
    # scale data
    if (T) {
      df <- data.frame(t(scale(t(mat))))
    }
    
    # k-means
    if (T) {
      cluster_kmeans <- function(df, k) {
        set.seed(799)
        fit.km <- kmeans(df, centers=k, nstart=25)
        cluster <- fit.km$cluster
        cluster <- data.frame(state = names(cluster), cluster = as.numeric(cluster))
        return(cluster)
      }
      
      clusterK <- cluster_kmeans(df, k)
    }
    
    # add cluster info into ggdat
    if (T) {
      df$State <- rownames(df)
      ggdat <- melt(df, id.vars = "State", value.name = "area", variable.name = "pol2_sig")
      ggdat$cluster <- paste0("C", clusterK[match(ggdat$State, clusterK$state), 2])
    }
    
    # format the data
    if (T) {
      ggdat$State <- factor(ggdat$State, levels = paste0("S", sort(as.numeric(gsub("S", "", unique(ggdat$State))))))
      ggdat$cluster <- factor(ggdat$cluster, levels = paste0("C", 1:k))
      ggdat$pol2_sig <- factor(ggdat$pol2_sig, levels = paste0("Q", 1:pol2_levels))
    }
    
    # ggplot
    if (T) {
      head(ggdat)
      
      p <- ggplot(ggdat) +
        geom_line(aes(x=pol2_sig, y=area, group=State, color=State), linewidth=1, alpha=0.3) +
        ggtitle(title) +
        scale_x_discrete(breaks=NULL) +
        xlab(paste("Expression Level: 0", paste(paste0("Q", 1:10), collapse = "  ->  "), sep = "  ->  ")) +
        ylab("AUC (Area Under Distribution Curve)") +
        facet_grid(.~cluster, scale = "fixed") +
        cowplot::theme_cowplot() +
        theme(legend.position = "none", 
              panel.border = element_rect(color="black"), 
              strip.background = element_rect(fill=NA, color=NA))
    }
    
    # get merged line
    if (T) {
      merged_line <- tapply(ggdat$area, paste(ggdat$cluster, ggdat$pol2_sig), median)
      merged_line <- data.frame(id = names(merged_line), 
                                value = as.numeric(merged_line))
      merged_line$cluster <- str_split(merged_line$id, " ", simplify = T)[, 1]
      merged_line$pol2_sig <- str_split(merged_line$id, " ", simplify = T)[, 2]
      
      merged_line$cluster <- factor(merged_line$cluster, levels = levels(ggdat$cluster))
      merged_line$pol2_sig <- factor(merged_line$pol2_sig, levels = levels(ggdat$pol2_sig))
    }
    
    # add merged line into the ggplot
    if (T) {
      head(merged_line)
      merged_p <- p + ggnewscale::new_scale_color() +
        geom_line(data = merged_line, aes(x=pol2_sig, y=value, group=cluster, color=cluster), linewidth=2)
    }
    
    # get state label
    if (T) {
      state_cluster <- tapply(clusterK$state, paste0("C", clusterK$cluster), function(x) {
        # x <- clusterK$state[clusterK$cluster==1]
        x <- gsub("S", "", x)
        if (length(x) > 20) {
          x1 <- paste0(paste(x[1:10], collapse = ", "), "\n")
          x2 <- paste0(paste(x[11:20], collapse = ", "), "\n")
          x_remin <- paste(x[21:length(x)], collapse = ", ")
          x <- paste0(x1, x2, x_remin)
        } else if (length(x) > 10) {
          x1 <- paste0(paste(x[1:10], collapse = ", "), "\n")
          x_remin <- paste(x[11:length(x)], collapse = ", ")
          x <- paste0(x1, x_remin)
        } else {
          x <- paste(x, collapse = ", ")
        }
        paste0("S", x)
      })
      state_cluster <- data.frame(
        cluster = names(state_cluster), 
        label = state_cluster
      )
    }
    
    # add merged line into the ggplot
    if (T) {
      head(state_cluster)
      label_p <- merged_p + 
        geom_text(data = state_cluster, aes(x=1, y=max(ggdat$area), label=label), hjust = 0, vjust = 1)
    }
    
    out <- list(label_p, clusterK)
    return(out)
  }
}

# state auc tendency clustering
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir <- paste0("results/2.pic/53.state_AUC_acorss_pol2_sig_kmeans/", data_type)
      if (! dir.exists(dir)) {
        dir.create(dir, showWarnings = F, recursive = T)
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
          
          # state auc: pol2 sig only
          if (T) {
            # get auc value data
            if (T) {
              file <- paste0("data/saved_data/52.pol2_ppm/", data_type, "/52.cs_ppm_pol2_Q", pol2_levels, "_AUC_in_", data_type, "s.", cell, ".", type, ".qs")
              ggdat <- qread(file, nthreads = 6)
            }
            
            # read the cluster data
            if (T) {
              file <- "data/saved_data/25.common_cluster_name_for_2_cell_CS_clusters.qs"
              cluster <- qread(file, nthreads = 6)
              cluster$orig.ident <- paste0("S", cluster$orig.ident)
              cluster <- cluster[cluster$cell == cell, ]
            }
            
            # common value
            if (T) {
              all_pol2_sig <- paste0("Q", 1:pol2_levels)
              all_states <- sort(as.numeric(gsub("S", "", cluster$orig.ident)))
            }
            
            # format the data
            if (T) {
              mat <- dcast(ggdat, state~Quantile, value.var = "auc")
              rownames(mat) <- paste0("S", mat$state)
              mat <- mat[, -1]
              mat <- mat[, all_pol2_sig]
            }
            
            # cor heatmap
            if (T) {
              cor_dat <- cor(t(mat), method = "pearson")
              
              p <- plot_state_auc_cor_heatmap(cor_dat, cluster)
              
              file <- paste0(dir, "/53.CS_AUC_cor_heatmap_", data_type, "s.", cell, ".", type, ".pdf")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
            
            # plot merged line clustered by chromIDEAS
            if (T) {
              p <- plot_line_add_merged_line(mat, cluster)
              
              file <- paste0(dir, "/53.CS_AUC_clustered_by_chromIDEAS_", data_type, "s.", cell, ".", type, ".pdf")
              ggsave(filename = file, plot = p, width = 10, height = 3)
            }
            
            # plot merged line clustered by kmeans5
            if (T) {
              # read the cluster data
              if (T) {
                cluster <- read.table("results/1.tab/40.emission_kmeans_cluster.csv", header = T, sep = ",", fill = T, comment.char = "")
                cluster <- cluster[, c("state", "cluster5")]
                colnames(cluster) <- c("orig.ident", "seurat_clusters")
              }
              
              p <- plot_line_add_merged_line(mat, cluster)
              
              file <- paste0(dir, "/53.CS_AUC_clustered_by_kmeans5cluster_", data_type, "s.", cell, ".", type, ".pdf")
              ggsave(filename = file, plot = p, width = 10, height = 3)
            }
            
            # kmeans
            if (T) {
              # determine the k value
              if (T) {
                p <- kmeans_determine_k(mat, title=paste0("pol2 AUC kmeans K"))
                
                file <- paste0(dir, "/53.CS_AUC_determine_K_", data_type, "s.", cell, ".", type, ".pdf")
                ggsave(filename = file, plot = p, width = 8, height = 6)
              }
              
              k_auc <- 3:6
              
              # kmeans distribution
              if (T) {
                kmeans_c <- lapply(k_auc, function(k) {
                  # k <- 5
                  dat <- plot_state_auc_distribution_kmeans(mat, k, title=paste0("kmeans CS AUC: k=", k))
                  
                  # get figure
                  file <- paste0(dir, "/53.CS_AUC_clustered_by_kmeans", k, "_", data_type, "s.", cell, ".", type, ".pdf")
                  ggsave(filename = file, plot = dat[[1]], width = 10, height = 3)
                  
                  # get cluster
                  c <- dat[[2]]
                  colnames(c)[2] <- paste0("K", k)
                  
                  return(c)
                })
                kmeans_c <- do.call(cbind, kmeans_c)
                
                apply(kmeans_c[, seq(1, ncol(kmeans_c), 2)], 1, table)
                
                kmeans_c <- kmeans_c[, c(1, seq(2, ncol(kmeans_c), 2))]
                
                write.table(kmeans_c, file = paste0("results/1.tab/53.CS_AUC_clustered_by_kmeans_", data_type, "s.", cell, ".", type, ".csv"),
                            quote = F, sep = ",", col.names = T, row.names = F)
              }
            }
          }
          
          k <- 4
          
          # state auc: pol2 sig + pol2 pat
          if (T) {
            # get auc value data
            if (T) {
              file <- paste0("data/saved_data/52.pol2_ppm/", data_type, "/52.cs_ppm_pol2_Q", pol2_levels, "_pol2_", k, "pat_AUC_in_", data_type, "s.", cell, ".", type, ".qs")
              ggdat <- qread(file, nthreads = 6)
            }
            
            # read the cluster data
            if (T) {
              file <- "data/saved_data/25.common_cluster_name_for_2_cell_CS_clusters.qs"
              cluster <- qread(file, nthreads = 6)
              cluster$orig.ident <- paste0("S", cluster$orig.ident)
              cluster <- cluster[cluster$cell == cell, ]
            }
            
            # common value
            if (T) {
              all_pol2_pat <- paste0("pol2_", 1:k)
              all_pol2_sig <- paste0("Q", 1:pol2_levels)
              all_states <- sort(as.numeric(gsub("S", "", cluster$orig.ident)))
            }
            
            # format the data
            if (T) {
              mat <- dcast(ggdat, ID~state, value.var = "auc")
              colnames(mat) <- paste0("S", colnames(mat))
              mat$pol2_sig <- str_split(mat[, 1], "@", simplify = T)[, 1]
              mat$pol2_pat <- str_split(mat[, 1], "@", simplify = T)[, 2]
              mat$pol2_pat <- paste0("pol2_", mat$pol2_pat)
              mat <- mat[, -1]
            }
            
            # cor heatmap
            if (T) {
              for (pol2_pat in all_pol2_pat) {
                # pol2_pat <- all_pol2_pat[1]
                
                cor_dat <- cor(mat[mat$pol2_pat == pol2_pat, paste0("S", all_states)], method = "pearson")
                
                p <- plot_state_auc_cor_heatmap(cor_dat, cluster)
                
                file <- paste0(dir, "/53.CS_AUC_cor_heatmap_", data_type, "s.", pol2_pat, ".", cell, ".", type, ".pdf")
                ggsave(filename = file, plot = p, width = 15, height = 10)
              }
            }
            
            # plot merged line clustered by chromIDEAS
            if (T) {
              for (pol2_pat in all_pol2_pat) {
                # pol2_pat <- all_pol2_pat[1]
                
                subdat <- mat[mat$pol2_pat == pol2_pat, ]
                rownames(subdat) <- subdat$pol2_sig
                subdat <- data.frame(t(subdat[, paste0("S", all_states)]))
                
                p <- plot_line_add_merged_line(subdat, cluster)
                
                file <- paste0(dir, "/53.CS_AUC_clustered_by_chromIDEAS_", data_type, "s.", pol2_pat, ".", cell, ".", type, ".pdf")
                ggsave(filename = file, plot = p, width = 10, height = 3)
              }
            }
            
            # plot merged line clustered by kmeans5
            if (T) {
              # read the cluster data
              if (T) {
                cluster <- read.table("results/1.tab/40.emission_kmeans_cluster.csv", header = T, sep = ",", fill = T, comment.char = "")
                cluster <- cluster[, c("state", "cluster5")]
                colnames(cluster) <- c("orig.ident", "seurat_clusters")
              }
              
              for (pol2_pat in all_pol2_pat) {
                # pol2_pat <- all_pol2_pat[1]
                
                subdat <- mat[mat$pol2_pat == pol2_pat, ]
                rownames(subdat) <- subdat$pol2_sig
                subdat <- data.frame(t(subdat[, paste0("S", all_states)]))
                
                p <- plot_line_add_merged_line(subdat, cluster)
                
                file <- paste0(dir, "/53.CS_AUC_clustered_by_kmeans5cluster_", data_type, "s.", pol2_pat, ".", cell, ".", type, ".pdf")
                ggsave(filename = file, plot = p, width = 10, height = 3)
              }
            }
            
            # kmeans
            if (T) {
              for (pol2_pat in all_pol2_pat) {
                # pol2_pat <- all_pol2_pat[1]
                
                subdat <- mat[mat$pol2_pat == pol2_pat, ]
                rownames(subdat) <- subdat$pol2_sig
                subdat <- data.frame(t(subdat[, paste0("S", all_states)]))
                
                # determine the k value
                if (T) {
                  p <- kmeans_determine_k(subdat, title=paste0("pol2 AUC kmeans K"))
                  
                  file <- paste0(dir, "/53.CS_AUC_determine_K_", data_type, "s.", pol2_pat, ".", cell, ".", type, ".pdf")
                  ggsave(filename = file, plot = p, width = 8, height = 6)
                }
                
                k_auc <- 3:6
                
                # kmeans distribution
                if (T) {
                  kmeans_c <- lapply(k_auc, function(k) {
                    # k <- 5
                    dat <- plot_state_auc_distribution_kmeans(subdat, k, title=paste0("kmeans CS AUC: k=", k))
                    
                    # get figure
                    file <- paste0(dir, "/53.CS_AUC_clustered_by_kmeans", k, "_", data_type, "s.", pol2_pat, ".", cell, ".", type, ".pdf")
                    ggsave(filename = file, plot = dat[[1]], width = 10, height = 3)
                    
                    # get cluster
                    c <- dat[[2]]
                    colnames(c)[2] <- paste0("K", k)
                    
                    return(c)
                  })
                  kmeans_c <- do.call(cbind, kmeans_c)
                  
                  apply(kmeans_c[, seq(1, ncol(kmeans_c), 2)], 1, table)
                  
                  kmeans_c <- kmeans_c[, c(1, seq(2, ncol(kmeans_c), 2))]
                  
                  write.table(kmeans_c, file = paste0("results/1.tab/53.CS_AUC_clustered_by_kmeans_", data_type, "s.", pol2_pat, ".", cell, ".", type, ".csv"),
                              quote = F, sep = ",", col.names = T, row.names = F)
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
        
        # state auc: pol2 sig only
        if (T) {
          # get auc value data
          if (T) {
            file <- paste0("data/saved_data/52.pol2_ppm/", data_type, "/52.cs_ppm_pol2_Q", pol2_levels, "_AUC_in_", data_type, "s.merged.", type, ".qs")
            ggdat <- qread(file, nthreads = 6)
          }
          
          # read the cluster data
          if (T) {
            file <- "data/saved_data/29.common_cluster_name_for_merged_2cells_CS_clusters.qs"
            cluster <- qread(file, nthreads = 6)
          }
          
          # common value
          if (T) {
            all_pol2_sig <- paste0("Q", 1:pol2_levels)
            all_states <- sort(as.numeric(gsub("S", "", cluster$orig.ident)))
          }
          
          # format the data
          if (T) {
            mat <- dcast(ggdat, state~Quantile, value.var = "auc")
            rownames(mat) <- paste0("S", mat$state)
            mat <- mat[, -1]
            mat <- mat[, all_pol2_sig]
          }
          
          # cor heatmap
          if (T) {
            cor_dat <- cor(t(mat), method = "pearson")
            
            p <- plot_state_auc_cor_heatmap(cor_dat, cluster)
            
            file <- paste0(dir, "/53.CS_AUC_cor_heatmap_", data_type, "s.merged.", type, ".pdf")
            ggsave(filename = file, plot = p, width = 15, height = 10)
          }
          
          # plot merged line clustered by chromIDEAS
          if (T) {
            p <- plot_line_add_merged_line(mat, cluster)
            
            file <- paste0(dir, "/53.CS_AUC_clustered_by_chromIDEAS_", data_type, "s.merged.", type, ".pdf")
            ggsave(filename = file, plot = p, width = 10, height = 3)
          }
          
          # plot merged line clustered by kmeans5
          if (T) {
            # read the cluster data
            if (T) {
              cluster <- read.table("results/1.tab/40.emission_kmeans_cluster.csv", header = T, sep = ",", fill = T, comment.char = "")
              cluster <- cluster[, c("state", "cluster5")]
              colnames(cluster) <- c("orig.ident", "seurat_clusters")
            }
            
            p <- plot_line_add_merged_line(mat, cluster)
            
            file <- paste0(dir, "/53.CS_AUC_clustered_by_kmeans5cluster_", data_type, "s.merged.", type, ".pdf")
            ggsave(filename = file, plot = p, width = 10, height = 3)
          }
          
          # kmeans
          if (T) {
            # determine the k value
            if (T) {
              p <- kmeans_determine_k(mat, title=paste0("pol2 AUC kmeans K"))
              
              file <- paste0(dir, "/53.CS_AUC_determine_K_", data_type, "s.merged.", type, ".pdf")
              ggsave(filename = file, plot = p, width = 8, height = 6)
            }
            
            k_auc <- 3:6
            
            # kmeans distribution
            if (T) {
              kmeans_c <- lapply(k_auc, function(k) {
                # k <- 5
                dat <- plot_state_auc_distribution_kmeans(mat, k, title=paste0("kmeans CS AUC: k=", k))
                
                # get figure
                file <- paste0(dir, "/53.CS_AUC_clustered_by_kmeans", k, "_", data_type, "s.merged.", type, ".pdf")
                ggsave(filename = file, plot = dat[[1]], width = 10, height = 3)
                
                # get cluster
                c <- dat[[2]]
                colnames(c)[2] <- paste0("K", k)
                
                return(c)
              })
              kmeans_c <- do.call(cbind, kmeans_c)
              
              apply(kmeans_c[, seq(1, ncol(kmeans_c), 2)], 1, table)
              
              kmeans_c <- kmeans_c[, c(1, seq(2, ncol(kmeans_c), 2))]
              
              write.table(kmeans_c, file = paste0("results/1.tab/53.CS_AUC_clustered_by_kmeans_", data_type, "s.merged.", type, ".csv"),
                          quote = F, sep = ",", col.names = T, row.names = F)
            }
          }
        }
        
        k <- 4
        
        # state auc: pol2 sig + pol2 pat
        if (T) {
          # get auc value data
          if (T) {
            file <- paste0("data/saved_data/52.pol2_ppm/", data_type, "/52.cs_ppm_pol2_Q", pol2_levels, "_pol2_", k, "pat_AUC_in_", data_type, "s.merged.", type, ".qs")
            ggdat <- qread(file, nthreads = 6)
          }
          
          # read the cluster data
          if (T) {
            file <- "data/saved_data/29.common_cluster_name_for_merged_2cells_CS_clusters.qs"
            cluster <- qread(file, nthreads = 6)
          }
          
          # common value
          if (T) {
            all_pol2_pat <- paste0("pol2_", 1:k)
            all_pol2_sig <- paste0("Q", 1:pol2_levels)
            all_states <- sort(as.numeric(gsub("S", "", cluster$orig.ident)))
          }
          
          # format the data
          if (T) {
            mat <- dcast(ggdat, ID~state, value.var = "auc")
            colnames(mat) <- paste0("S", colnames(mat))
            mat$pol2_sig <- str_split(mat[, 1], "@", simplify = T)[, 1]
            mat$pol2_pat <- str_split(mat[, 1], "@", simplify = T)[, 2]
            mat$pol2_pat <- paste0("pol2_", mat$pol2_pat)
            mat <- mat[, -1]
          }
          
          # cor heatmap
          if (T) {
            for (pol2_pat in all_pol2_pat) {
              # pol2_pat <- all_pol2_pat[1]
              
              cor_dat <- cor(mat[mat$pol2_pat == pol2_pat, paste0("S", all_states)], method = "pearson")
              
              p <- plot_state_auc_cor_heatmap(cor_dat, cluster)
              
              file <- paste0(dir, "/53.CS_AUC_cor_heatmap_", data_type, "s.", pol2_pat, ".merged.", type, ".pdf")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
          }
          
          # plot merged line clustered by chromIDEAS
          if (T) {
            for (pol2_pat in all_pol2_pat) {
              # pol2_pat <- all_pol2_pat[1]
              
              subdat <- mat[mat$pol2_pat == pol2_pat, ]
              rownames(subdat) <- subdat$pol2_sig
              subdat <- data.frame(t(subdat[, paste0("S", all_states)]))
              
              p <- plot_line_add_merged_line(subdat, cluster)
              
              file <- paste0(dir, "/53.CS_AUC_clustered_by_chromIDEAS_", data_type, "s.", pol2_pat, ".merged.", type, ".pdf")
              ggsave(filename = file, plot = p, width = 10, height = 3)
            }
          }
          
          # plot merged line clustered by kmeans5
          if (T) {
            # read the cluster data
            if (T) {
              cluster <- read.table("results/1.tab/40.emission_kmeans_cluster.csv", header = T, sep = ",", fill = T, comment.char = "")
              cluster <- cluster[, c("state", "cluster5")]
              colnames(cluster) <- c("orig.ident", "seurat_clusters")
            }
            
            for (pol2_pat in all_pol2_pat) {
              # pol2_pat <- all_pol2_pat[1]
              
              subdat <- mat[mat$pol2_pat == pol2_pat, ]
              rownames(subdat) <- subdat$pol2_sig
              subdat <- data.frame(t(subdat[, paste0("S", all_states)]))
              
              p <- plot_line_add_merged_line(subdat, cluster)
              
              file <- paste0(dir, "/53.CS_AUC_clustered_by_kmeans5cluster_", data_type, "s.", pol2_pat, ".merged.", type, ".pdf")
              ggsave(filename = file, plot = p, width = 10, height = 3)
            }
          }
          
          # kmeans
          if (T) {
            for (pol2_pat in all_pol2_pat) {
              # pol2_pat <- all_pol2_pat[1]
              
              subdat <- mat[mat$pol2_pat == pol2_pat, ]
              rownames(subdat) <- subdat$pol2_sig
              subdat <- data.frame(t(subdat[, paste0("S", all_states)]))
              
              # determine the k value
              if (T) {
                p <- kmeans_determine_k(subdat, title=paste0("pol2 AUC kmeans K"))
                
                file <- paste0(dir, "/53.CS_AUC_determine_K_", data_type, "s.", pol2_pat, ".merged.", type, ".pdf")
                ggsave(filename = file, plot = p, width = 8, height = 6)
              }
              
              k_auc <- 3:6
              
              # kmeans distribution
              if (T) {
                kmeans_c <- lapply(k_auc, function(k) {
                  # k <- 5
                  dat <- plot_state_auc_distribution_kmeans(subdat, k, title=paste0("kmeans CS AUC: k=", k))
                  
                  # get figure
                  file <- paste0(dir, "/53.CS_AUC_clustered_by_kmeans", k, "_", data_type, "s.", pol2_pat, ".merged.", type, ".pdf")
                  ggsave(filename = file, plot = dat[[1]], width = 10, height = 3)
                  
                  # get cluster
                  c <- dat[[2]]
                  colnames(c)[2] <- paste0("K", k)
                  
                  return(c)
                })
                kmeans_c <- do.call(cbind, kmeans_c)
                
                apply(kmeans_c[, seq(1, ncol(kmeans_c), 2)], 1, table)
                
                kmeans_c <- kmeans_c[, c(1, seq(2, ncol(kmeans_c), 2))]
                
                write.table(kmeans_c, file = paste0("results/1.tab/53.CS_AUC_clustered_by_kmeans_", data_type, "s.", pol2_pat, ".merged.", type, ".csv"),
                            quote = F, sep = ",", col.names = T, row.names = F)
              }
            }
          }
        }
      }
    }
  }
}


