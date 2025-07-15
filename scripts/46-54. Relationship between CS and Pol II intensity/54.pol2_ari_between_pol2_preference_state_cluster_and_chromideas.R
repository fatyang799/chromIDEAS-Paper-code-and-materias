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
pol2_levels <- 10
k_auc <- 3:6
n_ct <- 5

# define functions
if (T) {
  # calculate the similarity between 2 cluster result (script40)
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

# adjust rand index radar to assess the superiority of various cluster results
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir <- paste0("results/2.pic/54.pol2_ari/", data_type)
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
            # get auc distribution line cluster (from scripts 33): as the golden standard
            if (T) {
              file <- paste0("results/1.tab/53.CS_AUC_clustered_by_kmeans_", data_type, "s.", cell, ".", type, ".csv")
              cluster_aucDistribution <- read.table(file, header = T, sep = ",", fill = T, comment.char = "")
            }
            
            # using kmeans to cluster emission table
            if (T) {
              cluster_emission_kmeans <- read.table("results/1.tab/40.emission_kmeans_cluster.csv", header = T, sep = ",")
            }
            
            # read the chromIDEAS cluster
            if (T) {
              cluster_chrom <- qread("data/saved_data/25.common_cluster_name_for_2_cell_CS_clusters.qs", nthreads = 6)
              cluster_chrom$orig.ident <- paste0("S", cluster_chrom$orig.ident)
              cluster_chrom <- cluster_chrom[cluster_chrom$cell == cell, ]
            }
            
            # get random group as control
            if (T) {
              file <- paste0("data/saved_data/40.random_cluster_control_n", n_ct, ".", cell, ".qs")
              cluster_random <- qread(file, nthreads = 6)
            }
            
            # make the order of all clusters are identical
            if (T) {
              cluster_emission_kmeans <- cluster_emission_kmeans[match(cluster_aucDistribution$state, cluster_emission_kmeans$state), ]
              cluster_chrom <- cluster_chrom[match(cluster_aucDistribution$state, cluster_chrom$orig.ident), ]
              cluster_random <- cluster_random[match(cluster_aucDistribution$state, cluster_random$state), ]
              
              cluster_compare <- cbind(cluster_emission_kmeans, cluster_chrom, cluster_random)
              cluster_compare <- cluster_compare[, c("cluster5", paste0("ct", 1:n_ct), "seurat_clusters")]
              colnames(cluster_compare) <- c("kmeans5", paste0("CT", 1:n_ct), "chromIDEAS")
              
              rm(cluster_emission_kmeans, cluster_chrom, cluster_random)
            }
            
            # compare the similarity compared with auc distribution line cluster
            if (T) {
              # setting
              if (T) {
                file <- paste0(dir, "/54.pol2_ARI_in_", data_type, "s.", cell, ".", type, ".pdf")
                pdf(file = file, width = 8, height = 6)
                colors <- c("#00AFBB", "#E7B800", "#FC4E07", "#a7c957")
                par(mfrow = c(2,2))
              }
              
              # statistics ari
              if (T) {
                ari_stat <- lapply(k_auc, function(kk) {
                  similarity <- similarity_between_clusters(cluster_aucDistribution[, paste0("K", kk)], cluster_compare)
                  ari <- similarity$ari
                  names(ari) <- rownames(similarity)
                  
                  return(ari)
                })
                names(ari_stat) <- paste0("K", k_auc)
              }
              
              # plot
              if (T) {
                for (kk in k_auc) {
                  # setting
                  title <- paste0("Pol2 AUC cluster", kk)
                  
                  plot_radar_figure_with_similarity(ari_stat[[paste0("K", kk)]], 
                                                    title, 
                                                    colors[max(k_auc)+1-kk], 
                                                    c(paste0("CT", c(3:5)), "chromIDEAS", "kmeans5", paste0("CT", c(1:2))))
                }
              }
              
              # save the results
              if (T) {
                ari <- data.frame(t(do.call(rbind, ari_stat)))
                ari$method <- rownames(ari)
                rownames(ari) <- 1:nrow(ari)
                
                file <- paste0("results/1.tab/54.pol2_ARI_in_", data_type, "s.", cell, ".", type, ".csv")
                write.table(ari, file = file, quote = F, sep = ",", col.names = T, row.names = F)
              }
              
              Sys.sleep(0.05)
              dev.off()
            }
          }
          
          k <- 4
          
          # state auc: pol2 sig + pol2 pat
          if (T) {
            for (pol2_pat in paste0("pol2_", 1:k)) {
              # pol2_pat <- "pol2_1"
              
              # get auc distribution line cluster (from scripts 53): as the golden standard
              if (T) {
                file <- paste0("results/1.tab/53.CS_AUC_clustered_by_kmeans_", data_type, "s.", pol2_pat, ".", cell, ".", type, ".csv")
                
                cluster_aucDistribution <- read.table(file, header = T, sep = ",", fill = T, comment.char = "")
              }
              
              # using kmeans to cluster emission table
              if (T) {
                cluster_emission_kmeans <- read.table("results/1.tab/40.emission_kmeans_cluster.csv", header = T, sep = ",")
              }
              
              # read the chromIDEAS cluster
              if (T) {
                cluster_chrom <- qread("data/saved_data/25.common_cluster_name_for_2_cell_CS_clusters.qs", nthreads = 6)
                cluster_chrom$orig.ident <- paste0("S", cluster_chrom$orig.ident)
                cluster_chrom <- cluster_chrom[cluster_chrom$cell == cell, ]
              }
              
              # get random group as control
              if (T) {
                file <- paste0("data/saved_data/40.random_cluster_control_n", n_ct, ".", cell, ".qs")
                cluster_random <- qread(file, nthreads = 6)
              }
              
              # make the order of all clusters are identical
              if (T) {
                cluster_emission_kmeans <- cluster_emission_kmeans[match(cluster_aucDistribution$state, cluster_emission_kmeans$state), ]
                cluster_chrom <- cluster_chrom[match(cluster_aucDistribution$state, cluster_chrom$orig.ident), ]
                cluster_random <- cluster_random[match(cluster_aucDistribution$state, cluster_random$state), ]
                
                cluster_compare <- cbind(cluster_emission_kmeans, cluster_chrom, cluster_random)
                cluster_compare <- cluster_compare[, c("cluster5", paste0("ct", 1:n_ct), "seurat_clusters")]
                colnames(cluster_compare) <- c("kmeans5", paste0("CT", 1:n_ct), "chromIDEAS")
                
                rm(cluster_emission_kmeans, cluster_chrom, cluster_random)
              }
              
              # compare the similarity compared with auc distribution line cluster
              if (T) {
                # setting
                if (T) {
                  file <- paste0(dir, "/54.pol2_ARI_", pol2_pat, "_in_", data_type, "s.", cell, ".", type, ".pdf")
                  pdf(file = file, width = 8, height = 6)
                  colors <- c("#00AFBB", "#E7B800", "#FC4E07", "#a7c957")
                  par(mfrow = c(2,2))
                }
                
                # statistics ari
                if (T) {
                  ari_stat <- lapply(k_auc, function(kk) {
                    similarity <- similarity_between_clusters(cluster_aucDistribution[, paste0("K", kk)], cluster_compare)
                    ari <- similarity$ari
                    names(ari) <- rownames(similarity)
                    
                    return(ari)
                  })
                  names(ari_stat) <- paste0("K", k_auc)
                }
                
                # plot
                if (T) {
                  for (kk in k_auc) {
                    # setting
                    title <- paste0("Pol2 AUC cluster", kk)
                    
                    plot_radar_figure_with_similarity(ari_stat[[paste0("K", kk)]], 
                                                      title, 
                                                      colors[max(k_auc)+1-kk], 
                                                      c(paste0("CT", c(3:5)), "chromIDEAS", "kmeans5", paste0("CT", c(1:2))))
                  }
                }
                
                # save the results
                if (T) {
                  ari <- data.frame(t(do.call(rbind, ari_stat)))
                  ari$method <- rownames(ari)
                  rownames(ari) <- 1:nrow(ari)
                  
                  file <- paste0("results/1.tab/54.pol2_ARI_", pol2_pat, "_in_", data_type, "s.", cell, ".", type, ".csv")
                  write.table(ari, file = file, quote = F, sep = ",", col.names = T, row.names = F)
                }
                
                Sys.sleep(0.05)
                dev.off()
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
          # get auc distribution line cluster (from scripts 33): as the golden standard
          if (T) {
            file <- paste0("results/1.tab/53.CS_AUC_clustered_by_kmeans_", data_type, "s.merged.", type, ".csv")
            cluster_aucDistribution <- read.table(file, header = T, sep = ",", fill = T, comment.char = "")
          }
          
          # using kmeans to cluster emission table
          if (T) {
            cluster_emission_kmeans <- read.table("results/1.tab/40.emission_kmeans_cluster.csv", header = T, sep = ",")
          }
          
          # read the chromIDEAS cluster
          if (T) {
            cluster_chrom <- qread("data/saved_data/29.common_cluster_name_for_merged_2cells_CS_clusters.qs", nthreads = 6)
          }
          
          # get random group as control
          if (T) {
            file <- paste0("data/saved_data/40.random_cluster_control_n", n_ct, ".merged.qs")
            cluster_random <- qread(file, nthreads = 6)
          }
          
          # make the order of all clusters are identical
          if (T) {
            cluster_emission_kmeans <- cluster_emission_kmeans[match(cluster_aucDistribution$state, cluster_emission_kmeans$state), ]
            cluster_chrom <- cluster_chrom[match(cluster_aucDistribution$state, cluster_chrom$orig.ident), ]
            cluster_random <- cluster_random[match(cluster_aucDistribution$state, cluster_random$state), ]
            
            cluster_compare <- cbind(cluster_emission_kmeans, cluster_chrom, cluster_random)
            cluster_compare <- cluster_compare[, c("cluster5", paste0("ct", 1:n_ct), "seurat_clusters")]
            colnames(cluster_compare) <- c("kmeans5", paste0("CT", 1:n_ct), "chromIDEAS")
            
            rm(cluster_emission_kmeans, cluster_chrom, cluster_random)
          }
          
          # compare the similarity compared with auc distribution line cluster
          if (T) {
            # setting
            if (T) {
              file <- paste0(dir, "/54.pol2_ARI_in_", data_type, "s.merged.", type, ".pdf")
              pdf(file = file, width = 8, height = 6)
              colors <- c("#00AFBB", "#E7B800", "#FC4E07", "#a7c957")
              par(mfrow = c(2,2))
            }
            
            # statistics ari
            if (T) {
              ari_stat <- lapply(k_auc, function(kk) {
                similarity <- similarity_between_clusters(cluster_aucDistribution[, paste0("K", kk)], cluster_compare)
                ari <- similarity$ari
                names(ari) <- rownames(similarity)
                
                return(ari)
              })
              names(ari_stat) <- paste0("K", k_auc)
            }
            
            # plot
            if (T) {
              for (kk in k_auc) {
                # setting
                title <- paste0("Pol2 AUC cluster", kk)
                
                plot_radar_figure_with_similarity(ari_stat[[paste0("K", kk)]], 
                                                  title, 
                                                  colors[max(k_auc)+1-kk], 
                                                  c(paste0("CT", c(3:5)), "chromIDEAS", "kmeans5", paste0("CT", c(1:2))))
              }
            }
            
            # save the results
            if (T) {
              ari <- data.frame(t(do.call(rbind, ari_stat)))
              ari$method <- rownames(ari)
              rownames(ari) <- 1:nrow(ari)
              
              file <- paste0("results/1.tab/54.pol2_ARI_in_", data_type, "s.merged.", type, ".csv")
              write.table(ari, file = file, quote = F, sep = ",", col.names = T, row.names = F)
            }
            
            Sys.sleep(0.05)
            dev.off()
          }
        }
        
        k <- 4
        
        # state auc: pol2 sig + pol2 pat
        if (T) {
          for (pol2_pat in paste0("pol2_", 1:k)) {
            # pol2_pat <- "pol2_1"
            
            # get auc distribution line cluster (from scripts 53): as the golden standard
            if (T) {
              file <- paste0("results/1.tab/53.CS_AUC_clustered_by_kmeans_", data_type, "s.", pol2_pat, ".merged.", type, ".csv")
              
              cluster_aucDistribution <- read.table(file, header = T, sep = ",", fill = T, comment.char = "")
            }
            
            # using kmeans to cluster emission table
            if (T) {
              cluster_emission_kmeans <- read.table("results/1.tab/40.emission_kmeans_cluster.csv", header = T, sep = ",")
            }
            
            # read the chromIDEAS cluster
            if (T) {
              cluster_chrom <- qread("data/saved_data/29.common_cluster_name_for_merged_2cells_CS_clusters.qs", nthreads = 6)
            }
            
            # get random group as control
            if (T) {
              file <- paste0("data/saved_data/40.random_cluster_control_n", n_ct, ".merged.qs")
              cluster_random <- qread(file, nthreads = 6)
            }
            
            # make the order of all clusters are identical
            if (T) {
              cluster_emission_kmeans <- cluster_emission_kmeans[match(cluster_aucDistribution$state, cluster_emission_kmeans$state), ]
              cluster_chrom <- cluster_chrom[match(cluster_aucDistribution$state, cluster_chrom$orig.ident), ]
              cluster_random <- cluster_random[match(cluster_aucDistribution$state, cluster_random$state), ]
              
              cluster_compare <- cbind(cluster_emission_kmeans, cluster_chrom, cluster_random)
              cluster_compare <- cluster_compare[, c("cluster5", paste0("ct", 1:n_ct), "seurat_clusters")]
              colnames(cluster_compare) <- c("kmeans5", paste0("CT", 1:n_ct), "chromIDEAS")
              
              rm(cluster_emission_kmeans, cluster_chrom, cluster_random)
            }
            
            # compare the similarity compared with auc distribution line cluster
            if (T) {
              # setting
              if (T) {
                file <- paste0(dir, "/54.pol2_ARI_", pol2_pat, "_in_", data_type, "s.merged.", type, ".pdf")
                pdf(file = file, width = 8, height = 6)
                colors <- c("#00AFBB", "#E7B800", "#FC4E07", "#a7c957")
                par(mfrow = c(2,2))
              }
              
              # statistics ari
              if (T) {
                ari_stat <- lapply(k_auc, function(kk) {
                  similarity <- similarity_between_clusters(cluster_aucDistribution[, paste0("K", kk)], cluster_compare)
                  ari <- similarity$ari
                  names(ari) <- rownames(similarity)
                  
                  return(ari)
                })
                names(ari_stat) <- paste0("K", k_auc)
              }
              
              # plot
              if (T) {
                for (kk in k_auc) {
                  # setting
                  title <- paste0("Pol2 AUC cluster", kk)
                  
                  plot_radar_figure_with_similarity(ari_stat[[paste0("K", kk)]], 
                                                    title, 
                                                    colors[max(k_auc)+1-kk], 
                                                    c(paste0("CT", c(3:5)), "chromIDEAS", "kmeans5", paste0("CT", c(1:2))))
                }
              }
              
              # save the results
              if (T) {
                ari <- data.frame(t(do.call(rbind, ari_stat)))
                ari$method <- rownames(ari)
                rownames(ari) <- 1:nrow(ari)
                
                file <- paste0("results/1.tab/54.pol2_ARI_", pol2_pat, "_in_", data_type, "s.merged.", type, ".csv")
                write.table(ari, file = file, quote = F, sep = ",", col.names = T, row.names = F)
              }
              
              Sys.sleep(0.05)
              dev.off()
            }
          }
        }
      }
    }
  }
}


