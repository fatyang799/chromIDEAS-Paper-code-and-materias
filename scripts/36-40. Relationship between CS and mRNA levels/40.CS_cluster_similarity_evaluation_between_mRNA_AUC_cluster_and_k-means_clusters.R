# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(stringr)
  library(reshape2)
  library(fmsb)
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
emission_file <- "data/raw_data/1.emission/chromIDEAS.emission.txt"
n_ct <- 5

# define functions
if (T) {
  # k-means (script39)
  kmeans_determine_k <- function(emission, title) {
    # test dat
    if (F) {
      title <- "test"
    }
    
    # scale data
    if (T) {
      df <- data.frame(scale(emission))
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
  cluster_kmeans <- function(emission, k, state_order) {
    # test dat
    if (F) {
      k <- 4
    }
    
    # scale data
    if (T) {
      df <- data.frame(scale(emission))
    }
    
    # k-means
    if (T) {
      set.seed(799)
      fit.km <- kmeans(df, centers=k, nstart=25)
    }
    
    # summary the results
    if (T) {
      cluster <- fit.km$cluster
      cluster <- data.frame(state = names(cluster), cluster = as.numeric(cluster))
      
      cluster <- cluster[match(state_order, cluster$state), ]
      rownames(cluster) <- 1:nrow(cluster)
    }
    
    return(cluster)
  }
  
  # calculate the similarity between 2 cluster result
  load_auc_distribution_line_cluster <- function(cell, k, type, name, data_type, scale_type) {
    # test dat
    if (F) {
      k <- 2:5
    }
    
    # get all files
    if (T) {
      files <- paste0("results/1.tab/39.CS_", type, name, "_", cell, "_", data_type, "_", 
                      scale_type, "_AUC_kmeans_K", k, "_mRNA_Q", rna_levels, ".csv")
    }
    
    # read all files
    if (T) {
      dat <- lapply(files, function(x) {
        dat <- read.table(x, header = T, sep = ",")
        colnames(dat)[2] <- str_extract(basename(x), "K[0-9]{1,2}")
        
        return(dat)
      })
      
      dat <- do.call(cbind, dat)
    }
    
    # format the data
    if (T) {
      test <- sapply(dat[, seq(1, ncol(dat), 2)], function(x) {
        paste(x, collapse = "+")
      })
      test <- length(unique(test)) == 1
      
      if (test) {
        dat <- dat[, c(1, seq(2, ncol(dat), 2))]
      } else {
        stop("The state in auc distribution line cluster files are not identical, please check.")
      }
    }
    
    return(dat)
  }
  similarity_between_clusters <- function(cluster1, cluster2) {
    # calculate the ari and ri
    similarity <- data.frame(t(sapply(cluster2, function(c) {
      ari <- flexclust::randIndex(table(c, cluster1),correct = T)
      ri <- flexclust::randIndex(table(c, cluster1),correct = F)
      res <- c(ari, ri)
      names(res) <- c("ari", "ri")
      
      return(res)
    })))
    
    return(similarity)
  }
  plot_radar_figure_with_similarity <- function(similarity, title, color, order=NA) {
    # test dat
    if (F) {
      similarity <- rnorm(8)
      names(similarity) <- c(paste0("cluster", 4:10), "chromIDEAS")
      title <- "Consistency of clustering: ARI"
      order <- c(paste0("cluster", 7:10), "chromIDEAS", paste0("cluster", 4:6))
    }
    
    # format the data
    if (T) {
      ggdat <- data.frame(rbind(
        matrix(rep(max(similarity)*1.1, length(similarity)), ncol = length(similarity)), 
        matrix(rep(min(similarity)*0.9, length(similarity)), ncol = length(similarity)), 
        similarity
      ))
      rownames(ggdat) <- c("max", "min", "value")
      colnames(ggdat) <- names(similarity)
      
      if (sum(is.na(order)) == 0 & length(order) == length(similarity)) {
        ggdat <- ggdat[, order]
      }
    }
    
    # radarchart
    if (T) {
      op <- par(mar = c(1, 1, 1, 1))
      
      radarchart(ggdat, 
                 axistype = 1, # 0:5分别表示6类坐标轴方案
                 axislabcol = "black", # 坐标轴标签颜色
                 caxislabels = paste0(seq(0, 100, 25), "%"), # axistype=1时，自定义坐标轴标签，数量=seg+1
                 calcex = 1, # axistype=1时，坐标轴标签字体大小
                 paxislabels = NULL, # axistype=2时，自定义坐标轴标签，数量=ncol(df)
                 palcex = 1, # axistype=2时，坐标轴标签字体大小
                 
                 pty = 16, # 数据点形状
                 pcol = color, # 数据点之间连线颜色
                 plty = 1, # 数据点之间连线类型，1实线
                 plwd = 2, # 数据点之间连线宽度，默认为1
                 pdensity = NULL, # 数据点围成图形中，以斜线绘制密度，null表示不绘制斜线（不能用0，否则pfcol参数失效），1表示绘制1条斜线
                 pangle = 45, # 数据点围成图形中，以斜线绘制密度，设置斜线的角度
                 pfcol = scales::alpha(color, 0.5), # 数据点围成图形中，填充颜色
                 
                 seg = 4, # 网络圈数
                 cglty = 1, # 雷达背景线条类型，2表示虚线
                 cglwd = 0.8, # 雷达背景线条宽度
                 cglcol = "grey", # 雷达背景线条颜色
                 
                 title = title, 
                 na.itp = F, # 对于NA是否自动插值，默认F
                 centerzero = , # 是否scale，默认F
                 
                 vlabels = NULL, # 自定义每个特征的名字，默认为colnames(df)
                 vlcex = 1 # 自定义每个特征名字的字体大小
      )
      par(op)
    }
  }
}

# read and format emission dat
if (T) {
  # read the emission
  if (T) {
    emission <- read.table(emission_file, header = T, sep = "\t", fill = T, comment.char = "")
  }
  
  # load the ord
  if (T) {
    ord <- qread("data/saved_data/2.order_of_emission_table.qs", nthreads = 6)
  }
  
  # format the emission
  if (T) {
    rownames(emission) <- emission$State
    emission <- emission[ord[[1]], ]
    emission <- emission[, c("ATAC", "H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3", "H3K9me3", "H3K79me2")]
  }
  
  rm(ord)
}

# emission cluster by kmeans
if (T) {
  state_order <- paste0("S", (1:nrow(emission))-1)
  
  # k-means
  if (T) {
    # determine the k value: Within groups sum of squares
    if (T) {
      p <- kmeans_determine_k(emission, "K value for emission k-means clustering")
      ggsave(filename = "results/2.pic/40.kmeans_determine_k_value_for_emission_table.pdf", 
             plot = p, width = 6, height = 5)
    }
    
    # kmeams cluster
    if (T) {
      cluster4 <- cluster_kmeans(emission, 4, state_order)
      colnames(cluster4)[2] <- "cluster4"
      
      for (k in 5:10) {
        cluster <- cluster_kmeans(emission, k, state_order)
        cluster4 <- cbind(cluster4, cluster$cluster)
        colnames(cluster4)[ncol(cluster4)] <- paste0("cluster", k)
      }
      cluster <- cluster4
      
      rm(cluster4, k)
    }
    
    # rename cluster5
    if (T) {
      cluster$cluster5 <- ifelse(cluster$cluster5 == 2, 1, 
                                 ifelse(cluster$cluster5 == 4, 2, 
                                        ifelse(cluster$cluster5 == 5, 3, 
                                               ifelse(cluster$cluster5 == 1, 4, 
                                                      ifelse(cluster$cluster5 == 3, 5, NA)))))
    }
    
    write.table(cluster, file = "results/1.tab/40.emission_kmeans_cluster.csv", 
                quote = F, sep = ",", col.names = T, row.names = F)
  }
  
  # plot the clustered emission heatmap based on the k-means
  if (T) {
    ggdat <- as.matrix(emission)
    group <- cluster[match(rownames(ggdat), cluster$state), ]
    group <- sapply(group[, -1], function(x) {
      paste0("C", x)
    })
    
    # color setting
    if (T) {
      col_fun = colorRamp2(c(0, ceiling(quantile(c(ggdat), 0.95))),
                           c("#FFFFFF", "#3953A4"))
      
      color <- c("#2a9d8f", "#e9c46a", "#f4a261", "#e76f51", "#ff006e", "#8338ec", "#3a86ff", "#a3b18a", "#588157", "#3a5a40")
      names(color) <- paste0("C", 1:10)
    }
    
    # sorting
    if (T) {
      ord <- qread("data/saved_data/2.order_of_emission_table.qs", nthreads = 6)
      ggdat <- ggdat[, ord[[2]]]
    }
    
    # merged cluster results
    if (T) {
      # heatmap
      if (T) {
        p <- Heatmap(ggdat, name = "Probability", col = col_fun, 
                     border_gp = gpar(col = "black"), rect_gp = gpar(col = "black"),
                     show_heatmap_legend = TRUE, heatmap_legend_param = list(border = T),
                     row_names_max_width = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), row_names_rot = 0,
                     column_names_max_height = unit(6, "cm"), column_names_gp = gpar(fontsize = 10), column_names_rot = 60,
                     cluster_rows = F, cluster_columns = F, show_row_dend = F, show_column_dend = F, 
                     right_annotation = rowAnnotation(cluster = group, 
                                                      border = T, 
                                                      col = list(cluster = color),
                                                      annotation_legend_param = list(border = "black", 
                                                                                     at = paste0("C", 1:10))))
        p <- ggplotify::as.ggplot(p)
      }
      
      ggsave(plot = p, filename = "results/2.pic/40.emission_heatmap_kmeans_4-10_cluster_anno.pdf", width = 6, height = 10)
    }
    
    # separeated cluster result
    if (T) {
      # get separeated cluster result
      for (c in 4:10) {
        if (c<10) {
          p <- Heatmap(ggdat, name = "Probability", col = col_fun,
                       border_gp = gpar(col = "black"), rect_gp = gpar(col = "black"),
                       show_heatmap_legend = F, heatmap_legend_param = list(border = T),
                       row_names_max_width = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), row_names_rot = 0,
                       column_names_max_height = unit(6, "cm"), column_names_gp = gpar(fontsize = 10), column_names_rot = 60,
                       row_split = group[, paste0("cluster", c)],
                       cluster_rows = T, cluster_columns = F, show_row_dend = F, show_column_dend = F, 
                       clustering_distance_columns = "euclidean", clustering_distance_rows = "euclidean", 
                       clustering_method_columns = "ward.D2", clustering_method_rows = "ward.D2", 
                       right_annotation = rowAnnotation(cluster = group[, paste0("cluster", c)], 
                                                        border = T, 
                                                        col = list(cluster = color),
                                                        show_legend = F))
        }
        if (c==10) {
          p <- Heatmap(ggdat, name = "Probability", col = col_fun, 
                       border_gp = gpar(col = "black"), rect_gp = gpar(col = "black"),
                       show_heatmap_legend = T, heatmap_legend_param = list(border = T),
                       row_names_max_width = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), row_names_rot = 0,
                       column_names_max_height = unit(6, "cm"), column_names_gp = gpar(fontsize = 10), column_names_rot = 60,
                       row_split = group[, paste0("cluster", c)],
                       cluster_rows = T, cluster_columns = F, show_row_dend = F, show_column_dend = F, 
                       clustering_distance_columns = "euclidean", clustering_distance_rows = "euclidean", 
                       clustering_method_columns = "ward.D2", clustering_method_rows = "ward.D2", 
                       right_annotation = rowAnnotation(cluster = group[, paste0("cluster", c)], 
                                                        border = T, 
                                                        col = list(cluster = color),
                                                        annotation_legend_param = list(border = "black", 
                                                                                       at = paste0("C", 1:c))))
        }
        p <- ggplotify::as.ggplot(p)
        
        assign(paste0("p", c), p)
      }
      
      # patch the figure
      if (T) {
        library(patchwork)
        p <- p4+p5+p6+p7+p8+p9+p10+
          plot_layout(nrow=1, ncol=7, widths = c(rep(1, 6), 1.2))
        p
      }
      
      ggsave(plot = p, filename = "results/2.pic/40.emission_heatmap_kmeans_4-10_cluster_patched.pdf", width = 15, height = 10)
    }
  }
  
  rm(list = c("p", 
              paste0("p", 4:10), 
              "ggdat", "group", "c", "color", "col_fun", "cluster", "emission"))
}

# adjust rand index radar to assess the superiority of various cluster results
if (T) {
  # common used value
  if (T) {
    type <- "Body"
    data_type <- "tx"
    scale_type <- "genomic"
    k_auc <- 2:5
  }
  
  for (name in c("only", "UD")) {
    for (cell in c(cell1, cell2, "merged")) {
      # test dat
      if (F) {
        cell <- "merged"
        cell <- cell1
        name <- "UD"
      }
      
      # get auc distribution line cluster (from scripts 39): as the golden standard
      if (T) {
        cluster_aucDistribution <- load_auc_distribution_line_cluster(cell, k_auc, type, name, data_type, scale_type)
      }
      
      # using kmeans to cluster emission table
      if (T) {
        cluster_emission_kmeans <- read.table("results/1.tab/40.emission_kmeans_cluster.csv", 
                                              header = T, sep = ",")
      }
      
      # read the chromIDEAS cluster
      if (T) {
        file <- ifelse(cell == "merged", 
                       "data/saved_data/29.common_cluster_name_for_merged_2cells_CS_clusters.qs", 
                       "data/saved_data/25.common_cluster_name_for_2_cell_CS_clusters.qs")
        cluster_chrom <- qread(file, nthreads = 6)
        cluster_chrom <- cluster_chrom[cluster_chrom$cell == cell, ]
        cluster_chrom <- cluster_chrom[, c("orig.ident", "seurat_clusters", "cell")]
        
        if (cell != "merged") {
          cluster_chrom$orig.ident <- paste0("S", cluster_chrom$orig.ident)
        }
        
        cluster_chrom$orig.ident <- as.character(cluster_chrom$orig.ident)
        cluster_chrom$seurat_clusters <- as.character(cluster_chrom$seurat_clusters)
        cluster_chrom$cell <- as.character(cluster_chrom$cell)
      }
      
      # create random group as control
      if (T) {
        # get state and cluster data
        if (T) {
          all_cluster <- as.numeric(as.character(cluster_chrom$seurat_clusters))
        }
        
        # random ct
        if (T) {
          file <- paste0("data/saved_data/40.random_cluster_control_n", n_ct, ".", cell, ".qs")
          
          # get control
          if (T) {
            set.seed(799)
            
            cluster_random <- data.frame(
              sapply(1:n_ct, function(i) {
                # i <- 1
                sample(all_cluster, length(all_cluster), replace = F)
              })
            )
            colnames(cluster_random) <- paste0("ct", 1:n_ct)
            cluster_random <- data.frame(state = state_order, 
                                         cluster_random)
            cluster_random[, paste0("ct", 1:n_ct)] <- data.frame(
              sapply(cluster_random[, paste0("ct", 1:n_ct)], as.numeric)
            )
          }
          
          qsave(cluster_random, file = file, nthreads = 6)
        }
        
        rm(all_cluster)
      }
      
      # make the order of all clusters are identical
      if (T) {
        cluster_aucDistribution <- cluster_aucDistribution[match(state_order, cluster_aucDistribution$state), ]
        cluster_emission_kmeans <- cluster_emission_kmeans[match(state_order, cluster_emission_kmeans$state), ]
        cluster_chrom <- cluster_chrom[match(state_order, cluster_chrom$orig.ident), ]
        cluster_random <- cluster_random[match(state_order, cluster_random$state), ]
        
        cluster_compare <- cbind(cluster_emission_kmeans, cluster_chrom, cluster_random)
        cluster_compare <- cluster_compare[, c(paste0("cluster", 4:10), paste0("ct", 1:n_ct), "seurat_clusters")]
        colnames(cluster_compare) <- c(paste0("kmeans", 4:10), paste0("CT", 1:n_ct), "chromIDEAS")
        rm(cluster_emission_kmeans, cluster_chrom, cluster_random)
      }
      
      # compare the similarity compared with auc distribution line cluster (kmeans 4:10)
      if (T) {
        # setting
        if (T) {
          file <- paste0("results/2.pic/40.clustering_consistency_ARI_kmeans4-10_", type, name, "_", cell, "_", data_type, "_", scale_type, "_mRNA.pdf")
          pdf(file = file, width = 8, height = 6)
          colors <- c("#00AFBB", "#E7B800", "#FC4E07", "#a7c957")
          par(mfrow = c(2,2))
        }
        
        # statistics ari
        if (T) {
          ari_stat <- lapply(k_auc, function(k) {
            similarity <- similarity_between_clusters(cluster_aucDistribution[, paste0("K", k)], cluster_compare)
            ari <- similarity$ari
            names(ari) <- rownames(similarity)
            
            return(ari)
          })
          names(ari_stat) <- paste0("K", k_auc)
        }
        
        # plot
        if (T) {
          for (k in k_auc) {
            # setting
            title <- paste0("mRNA AUC cluster", k)
            
            plot_radar_figure_with_similarity(ari_stat[[paste0("K", k)]], 
                                              title, 
                                              color = colors[max(k_auc)+1-k], 
                                              order = c(paste0("CT", 1:n_ct), "chromIDEAS", paste0("kmeans", 4:10)))
          }
        }
        
        # save the results
        if (T) {
          ari <- data.frame(t(do.call(rbind, ari_stat)))
          ari$method <- rownames(ari)
          rownames(ari) <- 1:nrow(ari)
          
          file <- paste0("results/1.tab/40.mRNA_AUC_Body_", name, "_ARI_kmeans4-10_", cell, ".csv")
          write.table(ari, file = file, quote = F, sep = ",", col.names = T, row.names = F)
        }
        
        Sys.sleep(0.05)
        dev.off()
      }
      
      # select kmeans5 only to plot radar figure
      if (T) {
        cluster_compare <- cluster_compare[, c(paste0("kmeans", 5), paste0("CT", 1:n_ct), "chromIDEAS")]
      }
      
      # compare the similarity compared with auc distribution line cluster (kmeans 5)
      if (T) {
        # setting
        if (T) {
          file <- paste0("results/2.pic/40.clustering_consistency_ARI_", type, name, "_", cell, "_", data_type, "_", scale_type, "_mRNA.pdf")
          pdf(file = file, width = 8, height = 6)
          colors <- c("#00AFBB", "#E7B800", "#FC4E07", "#a7c957")
          par(mfrow = c(2,2))
        }
        
        # statistics ari
        if (T) {
          ari_stat <- lapply(k_auc, function(k) {
            similarity <- similarity_between_clusters(cluster_aucDistribution[, paste0("K", k)], cluster_compare)
            ari <- similarity$ari
            names(ari) <- rownames(similarity)
            
            return(ari)
          })
          names(ari_stat) <- paste0("K", k_auc)
        }
        
        # plot
        if (T) {
          for (k in k_auc) {
            # setting
            title <- paste0("mRNA AUC cluster", k)
            
            plot_radar_figure_with_similarity(ari_stat[[paste0("K", k)]], 
                                              title, 
                                              color = colors[max(k_auc)+1-k], 
                                              # order = c(paste0("CT", 1:n_ct), "chromIDEAS", paste0("kmeans", 4:10)))
                                              order = c(paste0("CT", c(3:5)), "chromIDEAS", "kmeans5", paste0("CT", c(1:2))))
          }
        }
        
        # save the results
        if (T) {
          ari <- data.frame(t(do.call(rbind, ari_stat)))
          ari$method <- rownames(ari)
          rownames(ari) <- 1:nrow(ari)
          
          file <- paste0("results/1.tab/40.mRNA_AUC_Body_", name, "_ARI_", cell, ".csv")
          write.table(ari, file = file, quote = F, sep = ",", col.names = T, row.names = F)
        }
        
        Sys.sleep(0.05)
        dev.off()
      }
    }
  }
}
