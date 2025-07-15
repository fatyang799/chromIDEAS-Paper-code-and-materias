# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(stringr)
  library(fmsb)
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
n_ct <- 5

# define functions
if (T) {
  # cell specific state preference for genomic regions (based on script10)
  state_preference_region <- function(value, index, genomics) {
    # test data
    if (F) {
      value <- state$ID
      index <- state$cd34
    }
    
    # get state
    if (T) {
      cell_state <- split(value, index)
    }
    
    # format the genomic regions
    if (T) {
      genomics$types <- paste0(genomics$data_source, "@", genomics$Type)
      all_types <- sort(unique(genomics$types))
    }
    
    # calculate the log2fc to represent the state preference
    if (T) {
      prefer <- lapply(names(cell_state), function(state_n) {
        # test data
        if (F) {
          state_n <- names(cell_state)[10]
        }
        
        # info print
        if (T) {
          cat(paste0("Now process S", state_n, ":\n"))
        }
        
        # calculate the preference to each region for the state
        if (T) {
          sub_state <- cell_state[[state_n]]
          
          region_prefer <- lapply(all_types, function(region) {
            # region <- all_types[1]
            mess <- paste0("\t(", which(all_types == region), "/", length(all_types), ") ", region, "\n")
            cat(mess)
            
            region_bins <- genomics$ID[genomics$types == region]
            
            # All DEGs
            n <- length(sub_state)
            
            # All genes
            N <- length(value)
            
            # M: All genes in TermA
            M <- length(region_bins)
            
            # k: DEGs in TermA
            k <- length(intersect(sub_state, region_bins))
            
            # enrichment P: phyper(k-1,M,N-M,n,lower.tail = F)
            p <- phyper(k-1, M, N-M, n, lower.tail = F)
            
            # fold change
            term_in_DEG <- k / n
            term_in_BG <- M / N
            fc <- term_in_DEG / term_in_BG
            
            res <- c(p, N, M, n, k, term_in_DEG, term_in_BG, fc, region)
            names(res) <- c("P_hyper", "All_genes_N", "All_genes_in_term_M", "All_DEGs_n", 
                            "DEGs_in_Term_k", "term_in_DEG", "term_in_BG", "FC", "region")
            
            return(res)
          })
        }
        
        # merge the data
        if (T) {
          region_prefer <- data.frame(do.call(rbind, region_prefer))
          region_prefer$State <- state_n
        }
        
        return(region_prefer)
      })
    }
    
    # merge and format the data
    if (T) {
      prefer <- data.frame(do.call(rbind, prefer))
      prefer[, 1:8] <- sapply(prefer[, 1:8], as.numeric)
      prefer$data_source <- str_split(prefer$region, "@", simplify = T)[, 1]
      prefer$regions <- str_split(prefer$region, "@", simplify = T)[, 2]
    }
    
    return(prefer)
  }
  state_preference_region2 <- function(value1, index1, value2, index2, genomics) {
    # get state
    if (T) {
      cell_state1 <- split(value1, index1)
      cell_state2 <- split(value2, index2)
    }
    
    # format the genomic regions
    if (T) {
      genomics$types <- paste0(genomics$data_source, "@", genomics$Type)
      all_types <- sort(unique(genomics$types))
    }
    
    # calculate the log2fc to represent the state preference
    if (T) {
      if (identical(names(cell_state1), names(cell_state2))) {
        all_states <- names(cell_state1)
      }
      
      prefer <- lapply(all_states, function(state_n) {
        # test data
        if (F) {
          state_n <- all_states[10]
        }
        
        # info print
        if (T) {
          cat(paste0("Now process S", state_n, ":\n"))
        }
        
        # calculate the preference to each region for the state
        if (T) {
          sub_state1 <- cell_state1[[state_n]]
          sub_state2 <- cell_state2[[state_n]]
          
          region_prefer <- lapply(all_types, function(region) {
            # region <- all_types[1]
            mess <- paste0("\t(", which(all_types == region), "/", length(all_types), ") ", region, "\n")
            cat(mess)
            
            region_bins <- genomics$ID[genomics$types == region]
            
            # All DEGs
            n <- length(sub_state1) + length(sub_state2)
            
            # All genes
            N <- length(value1) * 2
            
            # M: All genes in TermA
            M <- length(region_bins) * 2
            
            # k: DEGs in TermA
            k <- length(intersect(sub_state1, region_bins)) + length(intersect(sub_state2, region_bins))
            
            # enrichment P: phyper(k-1,M,N-M,n,lower.tail = F)
            p <- phyper(k-1, M, N-M, n, lower.tail = F)
            
            # fold change
            term_in_DEG <- k / n
            term_in_BG <- M / N
            fc <- term_in_DEG / term_in_BG
            
            res <- c(p, N, M, n, k, term_in_DEG, term_in_BG, fc, region)
            names(res) <- c("P_hyper", "All_genes_N", "All_genes_in_term_M", "All_DEGs_n", 
                            "DEGs_in_Term_k", "term_in_DEG", "term_in_BG", "FC", "region")
            
            return(res)
          })
        }
        
        # merge the data
        if (T) {
          region_prefer <- data.frame(do.call(rbind, region_prefer))
          region_prefer$State <- state_n
        }
        
        return(region_prefer)
      })
    }
    
    # merge and format the data
    if (T) {
      prefer <- data.frame(do.call(rbind, prefer))
      prefer[, 1:8] <- sapply(prefer[, 1:8], as.numeric)
      prefer$data_source <- str_split(prefer$region, "@", simplify = T)[, 1]
      prefer$regions <- str_split(prefer$region, "@", simplify = T)[, 2]
    }
    
    return(prefer)
  }
  
  # k-means (based on script40)
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
  
  # calculate the similarity between 2 cluster result: based on the script 34
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

# read the states data
if (T) {
  state_file <- "data/raw_data/2.states/chromIDEAS.state"
  
  state <- data.table::fread(state_file, sep = " ", header = T, data.table = F)
  state <- state[, c(1, 5, 6)]
  colnames(state)[1] <- "ID"
  head(state)
}

# load the genomic region bin ID
if (T) {
  file <- "data/saved_data/10.genomic_region_bin_IDs.qs"
  genomics <- qread(file, nthreads = 6)
}

# TS preference in genomic gtf region
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir1 <- "data/saved_data/45.GTF_preference/"
      if (! dir.exists(dir1)) {
        dir.create(dir1, showWarnings = F, recursive = T)
      }
      dir2 <- "results/2.pic/45.GTF_preference/"
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
          
          # get cell specific state prefer gtf regions
          if (T) {
            file <- paste0(dir1, "/45.state_prefer_genomic_region_stat_", type, ".", cell, ".qs")
            
            if (file.exists(file)) {
              cell_state_prefer <- qread(file, nthreads = 6)
            }
            if (! file.exists(file)) {
              cell_state_prefer <- state_preference_region(state$ID, state[, cell], genomics)
              
              qsave(cell_state_prefer, file = file, nthreads = 6)
            }
          }
          
          # format the data
          if (T) {
            head(cell_state_prefer)
            
            # convert into log2fc
            if (T) {
              non0min <- cell_state_prefer$FC
              non0min <- min(non0min[non0min>0])
              
              cell_state_prefer$Log2FC <- log2(cell_state_prefer$FC + non0min)
            }
            
            # long2wide
            if (T) {
              ggdat <- dcast(cell_state_prefer, State~region, value.var = "Log2FC")
              rownames(ggdat) <- paste0("S", ggdat[, 1])
              ggdat <- ggdat[, -1]
            }
            
            # format the colnames
            if (T) {
              colnames(ggdat) <- str_split(colnames(ggdat), "@", simplify = T)[, 2]
              colnames(ggdat) <- ifelse(colnames(ggdat) %in% c("cds", "tes", "tss", "utr3", "utr5"), toupper(colnames(ggdat)), 
                                        ifelse(colnames(ggdat) %in% c("exon", "intron", "intergenic", "repeats"), str_to_title(colnames(ggdat)), "CpG"))
              ggdat <- ggdat[, c("CpG", "TSS", "UTR5", "UTR3", "TES", "Exon", "CDS", "Intron", "Intergenic", "Repeats")]
            }
          }
          
          # kmeans clustering
          if (T) {
            # determine the K value
            if (T) {
              p <- kmeans_determine_k(ggdat, paste0("K value for genomic preference clustering in ", cell))
              ggsave(filename = paste0(dir2, "/45.wss_Kvalue_for_CS_gtf.", type, ".", cell, ".pdf"), 
                     plot = p, width = 6, height = 5)
            }
            
            k_genomic_distribution <- 3:6
            
            # get cluster
            if (T) {
              # k cluster for k_genomic_distribution (3:6)
              if (T) {
                cluster <- lapply(k_genomic_distribution, function(k) {
                  dat <- cluster_kmeans(ggdat, k, state_order = paste0("S", sort(as.numeric(gsub("S", "", rownames(ggdat))))))
                  colnames(dat)[2] <- paste0("K", k)
                  
                  return(dat)
                })
              }
              
              # merge the results
              if (T) {
                cluster <- data.frame(do.call(cbind, cluster))
                
                apply(cluster[, seq(1, ncol(cluster), 2)], 1, table)
                
                cluster <- cluster[, c("state", paste0("K", k_genomic_distribution))]
              }
              
              # save the results
              file <- paste0("results/1.tab/45.kmeans_CS_gtf_preference_for_", type, ".", cell, ".csv")
              write.table(cluster, file = file, quote = F, sep = ",", col.names = T, row.names = F)
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
        
        # get state specific prefer gtf regions
        if (T) {
          file <- paste0(dir1, "/45.state_prefer_genomic_region_stat_", type, ".merged.qs")
          
          if (file.exists(file)) {
            state_prefer_gtf <- qread(file, nthreads = 6)
          }
          if (! file.exists(file)) {
            state_prefer_gtf <- state_preference_region2(state$ID, state$cd34, state$ID, state$thp1, genomics)
            
            qsave(state_prefer_gtf, file = file, nthreads = 6)
          }
        }
        
        # format the data
        if (T) {
          head(state_prefer_gtf)
          
          # convert into log2fc
          if (T) {
            non0min <- state_prefer_gtf$FC
            non0min <- min(non0min[non0min>0])
            
            state_prefer_gtf$Log2FC <- log2(state_prefer_gtf$FC + non0min)
          }
          
          # long2wide
          if (T) {
            ggdat <- dcast(state_prefer_gtf, State~region, value.var = "Log2FC")
            rownames(ggdat) <- paste0("S", ggdat[, 1])
            ggdat <- ggdat[, -1]
          }
          
          # format the colnames
          if (T) {
            colnames(ggdat) <- str_split(colnames(ggdat), "@", simplify = T)[, 2]
            colnames(ggdat) <- ifelse(colnames(ggdat) %in% c("cds", "tes", "tss", "utr3", "utr5"), toupper(colnames(ggdat)), 
                                      ifelse(colnames(ggdat) %in% c("exon", "intron", "intergenic", "repeats"), str_to_title(colnames(ggdat)), "CpG"))
            ggdat <- ggdat[, c("CpG", "TSS", "UTR5", "UTR3", "TES", "Exon", "CDS", "Intron", "Intergenic", "Repeats")]
          }
        }
        
        # kmeans clustering
        if (T) {
          # determine the K value
          if (T) {
            p <- kmeans_determine_k(ggdat, paste0("K value for genomic preference clustering in merged"))
            ggsave(filename = paste0(dir2, "/45.wss_Kvalue_for_CS_gtf.", type, ".merged.pdf"), 
                   plot = p, width = 6, height = 5)
          }
          
          k_genomic_distribution <- 3:6
          
          # get cluster
          if (T) {
            # k cluster for k_genomic_distribution (3:6)
            if (T) {
              cluster <- lapply(k_genomic_distribution, function(k) {
                dat <- cluster_kmeans(ggdat, k, state_order = paste0("S", sort(as.numeric(gsub("S", "", rownames(ggdat))))))
                colnames(dat)[2] <- paste0("K", k)
                
                return(dat)
              })
            }
            
            # merge the results
            if (T) {
              cluster <- data.frame(do.call(cbind, cluster))
              
              apply(cluster[, seq(1, ncol(cluster), 2)], 1, table)
              
              cluster <- cluster[, c("state", paste0("K", k_genomic_distribution))]
            }
            
            # save the results
            file <- paste0("results/1.tab/45.kmeans_CS_gtf_preference_for_", type, ".merged.csv")
            write.table(cluster, file = file, quote = F, sep = ",", col.names = T, row.names = F)
          }
        }
      }
    }
  }
}

# adjust rand index radar to assess the superiority of various cluster results
if (T) {
  # value derive from above kmeans setting
  k_genomic_distribution <- 3:6
  
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir1 <- "data/saved_data/45.GTF_preference/"
      if (! dir.exists(dir1)) {
        dir.create(dir1, showWarnings = F, recursive = T)
      }
      dir2 <- "results/2.pic/45.GTF_preference/"
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
          
          # get kmeans clusters based on genomic preference for each state: as the golden standard
          if (T) {
            file <- paste0("results/1.tab/45.kmeans_CS_gtf_preference_for_", type, ".", cell, ".csv")
            cluster_golden <- read.table(file, header = T, sep = ",", fill = T, comment.char = "")
            
            state_order <- cluster_golden$state
          }
          
          # using kmeans to cluster emission table
          if (T) {
            cluster_emission_kmeans <- read.table("results/1.tab/40.emission_kmeans_cluster.csv", header = T, sep = ",")
            cluster_emission_kmeans <- cluster_emission_kmeans[, c("state", "cluster5")]
          }
          
          # read the chromIDEAS cluster
          if (T) {
            cluster_chrom <- qread("data/saved_data/25.common_cluster_name_for_2_cell_CS_clusters.qs", nthreads = 6)
            cluster_chrom <- cluster_chrom[cluster_chrom$cell == cell, ]
            cluster_chrom$orig.ident <- paste0("S", cluster_chrom$orig.ident)
          }
          
          # get random group control
          if (T) {
            file <- paste0("data/saved_data/40.random_cluster_control_n", n_ct, ".", cell, ".qs")
            cluster_random <- qread(file, nthreads = 6)
          }
          
          # make the order of all clusters are identical
          if (T) {
            cluster_golden <- cluster_golden[match(state_order, cluster_golden$state), ]
            cluster_emission_kmeans <- cluster_emission_kmeans[match(state_order, cluster_emission_kmeans$state), ]
            cluster_chrom <- cluster_chrom[match(state_order, cluster_chrom$orig.ident), ]
            cluster_random <- cluster_random[match(state_order, cluster_random$state), ]
            
            cluster_compare <- cbind(cluster_emission_kmeans, cluster_chrom, cluster_random)
            cluster_compare <- cluster_compare[, c(paste0("cluster", 5), paste0("ct", 1:n_ct), "seurat_clusters")]
            colnames(cluster_compare) <- c(paste0("kmeans", 5), paste0("CT", 1:n_ct), "chromIDEAS")
            
            rm(cluster_emission_kmeans, cluster_chrom, cluster_random)
          }
          
          # compare the similarity compared with auc distribution line cluster
          if (T) {
            # setting
            if (T) {
              file <- paste0(dir2, "/45.GTF_ARI_", data_type, "_", type, ".", cell, ".pdf")
              pdf(file = file, width = 8, height = 6)
              colors <- c("#00AFBB", "#E7B800", "#FC4E07", "#a7c957")
              par(mfrow = c(2,2))
            }
            
            # statistics ari
            if (T) {
              ari_stat <- lapply(k_genomic_distribution, function(k) {
                similarity <- similarity_between_clusters(cluster_golden[, paste0("K", k)], cluster_compare)
                ari <- similarity$ari
                names(ari) <- rownames(similarity)
                
                return(ari)
              })
              names(ari_stat) <- paste0("K", k_genomic_distribution)
            }
            
            # plot
            if (T) {
              for (k in k_genomic_distribution) {
                # setting
                title <- paste0("GTF Preference cluster", k)
                
                plot_radar_figure_with_similarity(ari_stat[[paste0("K", k)]], 
                                                  title, 
                                                  colors[max(k_auc)+1-k], 
                                                  c(paste0("CT", 1:n_ct), "chromIDEAS", paste0("kmeans", 4:10)))
              }
            }
            
            # save the results
            if (T) {
              ari <- data.frame(t(do.call(rbind, ari_stat)))
              ari$method <- rownames(ari)
              rownames(ari) <- 1:nrow(ari)
              
              file <- paste0("results/1.tab/45.gtf_preference_cluster_ARI_", type, ".", cell, ".csv")
              write.table(ari, file = file, quote = F, sep = ",", col.names = T, row.names = F)
            }
            
            Sys.sleep(0.05)
            dev.off()
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
        
        # get kmeans clusters based on genomic preference for each state: as the golden standard
        if (T) {
          file <- paste0("results/1.tab/45.kmeans_CS_gtf_preference_for_", type, ".merged.csv")
          cluster_golden <- read.table(file, header = T, sep = ",", fill = T, comment.char = "")
          
          state_order <- cluster_golden$state
        }
        
        # using kmeans to cluster emission table
        if (T) {
          cluster_emission_kmeans <- read.table("results/1.tab/40.emission_kmeans_cluster.csv", header = T, sep = ",")
          cluster_emission_kmeans <- cluster_emission_kmeans[, c("state", "cluster5")]
        }
        
        # read the chromIDEAS cluster
        if (T) {
          cluster_chrom <- qread("data/saved_data/29.common_cluster_name_for_merged_2cells_CS_clusters.qs", nthreads = 6)
        }
        
        # get random group control
        if (T) {
          file <- paste0("data/saved_data/40.random_cluster_control_n", n_ct, ".merged.qs")
          cluster_random <- qread(file, nthreads = 6)
        }
        
        # make the order of all clusters are identical
        if (T) {
          cluster_golden <- cluster_golden[match(state_order, cluster_golden$state), ]
          cluster_emission_kmeans <- cluster_emission_kmeans[match(state_order, cluster_emission_kmeans$state), ]
          cluster_chrom <- cluster_chrom[match(state_order, cluster_chrom$orig.ident), ]
          cluster_random <- cluster_random[match(state_order, cluster_random$state), ]
          
          cluster_compare <- cbind(cluster_emission_kmeans, cluster_chrom, cluster_random)
          cluster_compare <- cluster_compare[, c(paste0("cluster", 5), paste0("ct", 1:n_ct), "seurat_clusters")]
          colnames(cluster_compare) <- c(paste0("kmeans", 5), paste0("CT", 1:n_ct), "chromIDEAS")
          
          rm(cluster_emission_kmeans, cluster_chrom, cluster_random)
        }
        
        # compare the similarity compared with auc distribution line cluster
        if (T) {
          # setting
          if (T) {
            file <- paste0(dir2, "/45.GTF_ARI_", data_type, "_", type, ".merged.pdf")
            pdf(file = file, width = 8, height = 6)
            colors <- c("#00AFBB", "#E7B800", "#FC4E07", "#a7c957")
            par(mfrow = c(2,2))
          }
          
          # statistics ari
          if (T) {
            ari_stat <- lapply(k_genomic_distribution, function(k) {
              similarity <- similarity_between_clusters(cluster_golden[, paste0("K", k)], cluster_compare)
              ari <- similarity$ari
              names(ari) <- rownames(similarity)
              
              return(ari)
            })
            names(ari_stat) <- paste0("K", k_genomic_distribution)
          }
          
          # plot
          if (T) {
            for (k in k_genomic_distribution) {
              # setting
              title <- paste0("GTF Preference cluster", k)
              
              plot_radar_figure_with_similarity(ari_stat[[paste0("K", k)]], 
                                                title, 
                                                colors[7-k], 
                                                c(paste0("CT", 1:n_ct), "chromIDEAS", paste0("kmeans", 4:10)))
            }
          }
          
          # save the results
          if (T) {
            ari <- data.frame(t(do.call(rbind, ari_stat)))
            ari$method <- rownames(ari)
            rownames(ari) <- 1:nrow(ari)
            
            file <- paste0("results/1.tab/45.gtf_preference_cluster_ARI_", type, ".merged.csv")
            write.table(ari, file = file, quote = F, sep = ",", col.names = T, row.names = F)
          }
          
          Sys.sleep(0.05)
          dev.off()
        }
      }
    }
  }
}
