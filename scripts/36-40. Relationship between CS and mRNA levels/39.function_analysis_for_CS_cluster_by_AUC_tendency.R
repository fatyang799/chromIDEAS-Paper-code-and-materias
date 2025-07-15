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
rna_levels <- 10

# define functions
if (T) {
  # from script 38
  calculate_cell_specific_auc <- function(input, rm_up_down) {
    # test dat
    if (F) {
      input <- "data/saved_data/36.tx_level_profile_Body_5UP_5DW_10Body.bs200.format_genomic_percentage_mRNA_classification_Q10_mat.cd34.qs"
      rm_up_down <- F
    }
    
    # get parameters
    if (T) {
      type <- str_extract(basename(input), "TSS|Body")
      cell <- str_extract(basename(input), paste0(cell1, "|", cell2))
      scale_type <- str_split(str_extract(basename(input), "genomic_percentage|state_percentage"), "_", simplify = T)[, 1]
    }
    
    # get distribution data
    if (T) {
      dat <- qread(input, nthreads = 6)
    }
    
    # format the data
    if (T) {
      colnames(dat) <- c("State", "Loc", "Percentage", "mRNA")
      dat$group <- paste0(dat$mRNA, "@", dat$State)
      mat <- dcast(dat, group~Loc, value.var = "Percentage")
      mat$State <- str_split(mat$group, "@", simplify = T)[, 2]
      mat$mRNA <- str_split(mat$group, "@", simplify = T)[, 1]
      
      if (type == "Body") {
        order <- c(
          paste0("U", up_bin_num:1), 
          "TSS", 
          paste0("B", 1:body_bin_num), 
          "TES", 
          paste0("D", 1:down_bin_num), 
          "State", "mRNA"
        )
      }
      if (type == "TSS") {
        order <- c(
          paste0("U", up_bin_num:1), 
          "TSS", 
          paste0("D", 1:down_bin_num), 
          "State", "mRNA"
        )
      }
      mat <- mat[, order]
    }
    
    # common value
    if (T) {
      state_order <- paste0("S", sort(as.numeric(gsub("S", "", unique(mat$State)))))
      mrna_type <- c("subdat0", paste0("Q", 1:rna_levels))
    }
    
    # remove up and down stream for geneBody region
    pattern <- ifelse(rm_up_down & type == "Body", "TSS|^B|TES", "^U|TSS|^B|TES|^D")
    
    # calculation the AUC value
    if (T) {
      auc <- data.frame(t(sapply(state_order, function(s) {
        # s <- state_order[1]
        auc <- sapply(mrna_type, function(r) {
          # r <- mrna_type[1]
          line_dat <- mat[mat$State == s & mat$mRNA == r, grep(pattern, colnames(mat))]
          line_dat <- as.numeric(unlist(line_dat))
          
          # assume the x gap is 1
          area <- sapply(2:length(line_dat), function(n) {
            # n <- 2
            n2 <- line_dat[n]
            n1 <- line_dat[n-1]
            
            # The formula for calculating the area of a trapezoid
            auc <- (n1+n2) * 1 / 2
            return(auc)
          })
          auc <- sum(area)
          
          return(auc)
        })
        
        names(auc) <- mrna_type
        return(auc)
      })))
      
      auc$State <- state_order
    }
    
    # format the data
    if (F) {
      ggdat <- melt(auc, id.vars = "State", variable.name = "mRNA", value.name = "area")
      ggdat$mRNA <- factor(ggdat$mRNA, levels = mrna_type)
      ggdat$State <- factor(ggdat$State, levels = state_order)
    }
    
    return(auc)
  }
  
  # correlation
  plot_state_auc_cor_heatmap <- function(mat, cluster) {
    # correlation for auc
    if (T) {
      cor_dat <- cor(t(mat), method = "pearson")
    }
    
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
      p <- Heatmap(cor_dat, name = "Pearson", col = col_fun, column_title = paste0("The correlation between the AUC of chromatin states: ", toupper(unique(cluster$cell))), 
                   border_gp = gpar(col = "black"), rect_gp = gpar(col = "black"),
                   show_heatmap_legend = TRUE, heatmap_legend_param = list(border = T),
                   row_names_max_width = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), row_names_rot = 0,
                   column_names_max_height = unit(6, "cm"), column_names_gp = gpar(fontsize = 10), column_names_rot = 60,
                   cluster_rows = T, cluster_columns = T, show_row_dend = T, show_column_dend = F, 
                   clustering_distance_columns = "euclidean", clustering_distance_rows = "euclidean", 
                   clustering_method_columns = "ward.D2", clustering_method_rows = "ward.D2", 
                   right_annotation = rowAnnotation(Group = group, border = T, 
                                                    # gp = gpar(col = "black"), 
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
  
  # k-means
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
  plot_state_auc_distribution_kmeans <- function(mat, k, title) {
    # test dat
    if (F) {
      title <- "test"
      k <- 3
    }
    
    # scale data
    if (T) {
      df <- mat
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
      ggdat <- melt(df, id.vars = "State", value.name = "area", variable.name = "mRNA")
      ggdat$cluster <- paste0("C", clusterK[match(ggdat$State, clusterK$state), 2])
    }
    
    # format the data
    if (T) {
      ggdat$State <- factor(ggdat$State, levels = paste0("S", sort(as.numeric(gsub("S", "", unique(ggdat$State))))))
      ggdat$cluster <- factor(ggdat$cluster, levels = paste0("C", 1:k))
      ggdat$mRNA <- factor(ggdat$mRNA, levels = c("subdat0", paste0("Q", 1:rna_levels)))
    }
    
    # ggplot
    if (T) {
      head(ggdat)
      
      p <- ggplot(ggdat) +
        geom_line(aes(x=mRNA, y=area, group=State, color=State), linewidth=1, alpha=0.3) +
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
      merged_line <- tapply(ggdat$area, paste(ggdat$cluster, ggdat$mRNA), median)
      merged_line <- data.frame(id = names(merged_line), 
                                value = as.numeric(merged_line))
      merged_line$cluster <- str_split(merged_line$id, " ", simplify = T)[, 1]
      merged_line$mRNA <- str_split(merged_line$id, " ", simplify = T)[, 2]
      
      merged_line$cluster <- factor(merged_line$cluster, levels = levels(ggdat$cluster))
      merged_line$mRNA <- factor(merged_line$mRNA, levels = levels(ggdat$mRNA))
    }
    
    # add merged line into the ggplot
    if (T) {
      head(merged_line)
      merged_p <- p + ggnewscale::new_scale_color() +
        geom_line(data = merged_line, aes(x=mRNA, y=value, group=cluster, color=cluster), linewidth=2)
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
  
  move_file <- function(file, dir) {
    stat <- file.copy(from = file, to = paste0(dir, "/"), overwrite = T, copy.date = T)
    if (stat) {
      stat <- file.remove(file)
    }
  }
}

# plot the AUC distribution correlation
if (T) {
  data_type <- "tx"
  
  for (type in c("Body")) {
    cat(paste0(type, ":\n"))
    for (scale_type in c("genomic")) {
      cat(paste0("\t", scale_type, ":\n"))
      for (rm_up_down in c(T, F)) {
        # test data
        if (F) {
          type <- "Body"
          scale_type <- "genomic"
          rm_up_down <- F
          name <- ifelse(rm_up_down, "only", "UD")
          cell <- "merged"
          cell <- cell1
        }
        
        name <- ifelse(rm_up_down, "only", "UD")
        cat(paste0("\t\t", name, ":\n"))
        
        # mkdir
        if (T) {
          dir <- paste0("results/2.pic/39.stata_AUC_distribution_kmeans/Body", name)
          if (! dir.exists(dir)) {
            dir.create(dir, showWarnings = F, recursive = T)
          }
        }
        
        for (cell in c(cell1, cell2, "merged")) {
          cat(paste0("\t\t\t", cell, "\n"))
          
          # get input
          if (T) {
            input <- ifelse(type == "TSS", 
                            paste0("data/saved_data/36.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                   "DW.bs", bin_size, ".format_", scale_type, 
                                   "_percentage_mRNA_classification_Q", rna_levels, "_mat.", cell, ".qs"), 
                            paste0("data/saved_data/36.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                   "DW_", body_bin_num, "Body.bs", bin_size, ".format_", scale_type, 
                                   "_percentage_mRNA_classification_Q", rna_levels, "_mat.", cell, ".qs"))
          }
          
          # read the cluster data
          if (T) {
            file <- ifelse(cell == "merged", 
                           "data/saved_data/29.common_cluster_name_for_merged_2cells_CS_clusters.qs", 
                           "data/saved_data/25.common_cluster_name_for_2_cell_CS_clusters.qs")
            cluster <- qread(file, nthreads = 6)
            cluster <- cluster[cluster$cell == cell, ]
            cluster <- cluster[, c("orig.ident", "seurat_clusters", "cell")]
            
            if (cell != "merged") {
              cluster$orig.ident <- paste0("S", cluster$orig.ident)
            }
            
            cluster$orig.ident <- as.character(cluster$orig.ident)
            cluster$seurat_clusters <- as.character(cluster$seurat_clusters)
            cluster$cell <- as.character(cluster$cell)
          }
          
          # get auc value data
          if (T) {
            ggdat <- calculate_cell_specific_auc(input, rm_up_down)
          }
          
          # scale the data
          if (T) {
            mat <- ggdat[, colnames(ggdat) != "State"]
            mat <- data.frame(t(scale(t(mat))))
          }
          
          # common value
          if (T) {
            state_order <- paste0("S", sort(as.numeric(gsub("S", "", unique(ggdat$State)))))
            mrna_type <- c("subdat0", paste0("Q", 1:rna_levels))
          }
          
          # cor heatmap
          if (T) {
            # get output
            if (T) {
              file <- paste0("39.CS_", type, name, "_", cell, "_", data_type, "_", 
                             scale_type, "_AUC_corHeatmap_mRNA_Q", rna_levels, ".pdf")
            }
            
            # ggplot
            if (T) {
              p <- plot_state_auc_cor_heatmap(mat, cluster)
              
              ggsave(filename = file, plot = p, width = 10, height = 9)
              move_file(file, dir)
            }
          }
          
          # kmeans
          if (T) {
            # get output
            if (T) {
              file <- paste0("39.", type, name, "_", cell, "_", data_type, "_", 
                             scale_type, "_kmeans_Kvalue_mRNA_Q", rna_levels, ".pdf")
            }
            
            # determine the k value
            if (T) {
              p <- kmeans_determine_k(mat, paste0(cell, " Kmeans: ", type, " ", name))
              
              ggsave(filename = file, plot = p, width = 6, height = 5)
              move_file(file, dir)
            }
            
            # kmeans distribution
            if (T) {
              for (k in 2:5) {
                # k <- 5
                dat <- plot_state_auc_distribution_kmeans(mat, k, cell)
                
                p <- dat[[1]]
                # get output
                file <- paste0("39.", type, name, "_", cell, "_", data_type, "_", 
                               scale_type, "_AUC_kmeans_K", k, "_mRNA_Q", rna_levels, ".pdf")
                ggsave(filename = file, plot = p, width = 10, height = 3)
                move_file(file, dir)
                
                c <- dat[[2]]
                # get output
                file <- paste0("results/1.tab/39.CS_", type, name, "_", cell, "_", data_type, "_", 
                               scale_type, "_AUC_kmeans_K", k, "_mRNA_Q", rna_levels, ".csv")
                write.table(c, file = file, quote = F, sep = ",", col.names = T, row.names = F)
              }
            }
          }
        }
      }
    }
  }
}
