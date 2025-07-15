# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(reshape2)
  suppressPackageStartupMessages(library(stringr))
  library(qs)
}

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
body_bin_num <- 10
rna_levels <- 10
length_leveles <- 10
cells <- c(cell1, cell2)

# read the chromatin states cluster results
if (T) {
  file <- "data/saved_data/25.common_cluster_name_for_2_cell_CS_clusters.qs"
  cluster <- qread(file, nthreads = 6)
  cluster$seurat_clusters <- as.vector(cluster$seurat_clusters)
  cluster$orig.ident <- paste0("S", cluster$orig.ident)
}

# read the emission table
if (T) {
  emission <- read.table("data/raw_data/1.emission/chromIDEAS.emission.txt", 
                         header = T, sep = "\t")
  rownames(emission) <- emission[, 1]
  emission <- emission[, -c(1:2)]
}

# loop to get specific histone signal distribution
if (T) {
  deg_histone <- lapply(cells, function(cell) {
    # test data
    if (F) {
      cell <- cell1
    }
    
    # sub_cluster
    if (T) {
      cell_cluster <- cluster[cluster$cell == cell, ]
    }
    
    # format the data
    if (T) {
      stat_dat <- emission
      stat_dat <- stat_dat[cell_cluster$orig.ident, ]
      stat_dat$cluster <- cell_cluster$seurat_clusters
    }
    
    # deg analysis for histone signal
    if (T) {
      mks <- colnames(stat_dat)[-ncol(stat_dat)]
      clusters <- sort(as.numeric(unique(stat_dat$cluster)))
      
      deg_histone <- lapply(clusters, function(c) {
        subres <- lapply(mks, function(mk) {
          # test data
          if (F) {
            c <- 1
            mk <- mks[1]
          }
          
          target_dat <- stat_dat[, mk][as.numeric(stat_dat$cluster) == c]
          remain_dat <- stat_dat[, mk][as.numeric(stat_dat$cluster) != c]
          log2fc <- log2(mean(target_dat)/mean(remain_dat))
          type <- ifelse(log2fc>0, "greater", "less")
          
          p <- wilcox.test(x=target_dat, y=remain_dat, alternative=type)$p.value
          
          res <- data.frame(
            cluster = c, 
            mk = mk, 
            log2fc = log2fc, 
            p = p
          )
          
          return(res)
        })
        subres <- do.call(rbind, subres)
        return(subres)
      })
      deg_histone <- do.call(rbind, deg_histone)
    }
    
    deg_histone$cell <- cell
    
    return(deg_histone)
  })
  deg_histone <- do.call(rbind, deg_histone)
  
  file <- "results/1.tab/30.deg_epi.csv"
  write.table(deg_histone, file = file, quote = F, sep = ",", col.names = T, row.names = F)
}

# heatmap to show diff histone signal
if (T) {
  # filter the data
  if (T) {
    ggdat <- deg_histone
    ggdat$log2fc[ggdat$p>0.05] <- 0
  }
  
  # format the data
  if (T) {
    ggdat$mk <- paste0(ggdat$mk, "@", ggdat$cluster)
    head(ggdat)
    ggdat <- dcast(ggdat, cell~mk, value.var = "log2fc")
    rownames(ggdat) <- ggdat[, 1]
    ggdat <- ggdat[, -1]
    ggdat <- as.matrix(ggdat)
  }
  
  # color setting
  if (T) {
    library(RColorBrewer)
    display.brewer.all()
    colors <- brewer.pal(9, "Set1") 
    scales::show_col(colors, labels=T)
    
    range <- range(ggdat, na.rm = T)
    range <- sort(c(range, 0))
    col_fun <- colorRamp2(range, c(colors[2], "white", colors[1]))
  }
  
  # sorting 
  if (T) {
    colnames(ggdat)
    ggdat <- ggdat[, paste0(rep(c("ATAC", "H3K4me1", "H3K4me3", "H3K27ac", "H3K36me3", "H3K79me2", "H3K9me3", "H3K27me3"), length(unique(deg_histone$cluster))), 
                            "@", 
                            rep(sort(unique(deg_histone$cluster)), each=length(unique(deg_histone$mk))))]
  }
  
  # heatmap
  if (T) {
    p <- Heatmap(ggdat, name = "log2FC", col = col_fun, 
                 border_gp = gpar(col = "black"), rect_gp = gpar(col = "black"),
                 show_heatmap_legend = TRUE, heatmap_legend_param = list(border = T),
                 column_split = str_split(colnames(ggdat), "@", simplify = T)[, 2], cluster_row_slices = FALSE, na_col = "white",
                 row_names_max_width = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), row_names_rot = 0,
                 column_names_max_height = unit(6, "cm"), column_names_gp = gpar(fontsize = 10), column_names_rot = 0,
                 cluster_rows = F, cluster_columns = F, show_row_dend = F, show_column_dend = F)
  }
  
  # save the fig
  if (T) {
    p <- ggplotify::as.ggplot(p)
    
    ggsave(filename = "results/2.pic/30.cluster_histone_feature_heatmap.pdf", plot = p, width = 8, height = 8)
  }
  
  qsave(ggdat, file = "data/saved_data/30.histone_deg_heatmap_mat.qs", nthreads = 6)
  
  rm(p, ggdat, range, col_fun, colors)
}

# loop to do deg
if (T) {
  file <- "data/saved_data/30.deg_analysis_for_all_cluster_WNNmatrix_genes.qs"
  if (file.exists(file)) {
    results <- qread(file, nthreads = 6)
  }
  if (! file.exists(file)) {
    results <- lapply(cells, function(cell) {
      # test data
      if (F) {
        cell <- cell1
      }
      
      # get cluster info
      if (T) {
        cell_cluster <- cluster[cluster$cell == cell, ]
      }
      
      # get gene body chromatin states percentage matrix
      if (T) {
        gene_body <- qread(paste0("data/saved_data/23.seurat_gene_CS_percentage_input_mat_", cell, ".qs"), 
                           nthreads = 6)
        gene_body$n <- 1:nrow(gene_body)
      }
      
      # do differential gene analysis
      deg_all_clusters <- lapply(sort(unique(cell_cluster$seurat_clusters)), function(c) {
        # c <- sort(unique(cell_cluster$seurat_clusters))[1]
        # info
        if (T) {
          start_mess <- paste0("Now process differential gene analysis: ", cell, " cluster", c, " vs remain clusters\n")
          cat(start_mess)
          
          total_line_mess <- paste0("|", paste(rep("-", 100), collapse = ""), "|\n")
          cat(total_line_mess)
          
          breaks <- round(seq(1, nrow(gene_body), length.out=100))
        }
        
        # get group info 
        if (T) {
          s <- cell_cluster[cell_cluster$seurat_clusters == c, "orig.ident"]
          group <- ifelse(colnames(gene_body)[grep("^S[0-9]{1,2}", colnames(gene_body))] %in% s, paste0("group", c), "group_other")
        }
        
        # deg analysis Wilcoxon
        if (T) {
          deg_res <- data.frame(t(
            apply(gene_body, 1, function(x) {
              # x <- unlist(gene_body[1, ])
              # split the data
              if (T) {
                nrow <- as.numeric(x[length(x)])
                x <- x[-length(x)]
              }
              
              # print info
              if (nrow %in% breaks) {
                if (which(nrow == breaks) == 1) {
                  cat("|*")
                }
                if (which(nrow == breaks) == 100) {
                  cat("*|\n")
                }
                if (! which(nrow == breaks) %in% c(1, 100)) {
                  cat("*")
                }
              }
              
              # get represent value
              if (T) {
                target_group <- x[group!="group_other"]
                remain_group  <- x[group=="group_other"]
                
                target_group_med <- median(target_group)
                remain_group_med <- median(remain_group)
                
                target_group_mean <- mean(target_group)
                remain_group_mean <- mean(remain_group)
              }
              
              # wilcox p value calculation
              if (T) {
                type <- ifelse(target_group_mean>remain_group_mean, "greater", "less")
                
                result <- wilcox.test(target_group, remain_group, alternative = type, exact = F)
                p <- result$p.value
              }
              
              # get statistical value
              if (T) {
                target_group_max <- max(target_group)
                remain_group_max <- max(remain_group)
                
                target_group_min <- min(target_group)
                remain_group_min <- min(remain_group)
                
                target_group_sd <- sd(target_group)
                remain_group_sd <- sd(remain_group)
                
                target_group_n0 <- sum(target_group == 0)
                remain_group_n0 <- sum(remain_group == 0)
                
                target_group_non0 <- sum(target_group != 0)
                remain_group_non0 <- sum(remain_group != 0)
                
                target_group_sum <- sum(target_group)
                remain_group_sum <- sum(remain_group)
              }
              
              res <- c(p, target_group_med, remain_group_med, target_group_mean, remain_group_mean, 
                       target_group_max, target_group_min, target_group_sd, target_group_n0, target_group_non0, target_group_sum, 
                       remain_group_max, remain_group_min, remain_group_sd, remain_group_n0, remain_group_non0, remain_group_sum)
              names(res) <- c("P_value", "target_group_med", "remain_group_med", "target_group_mean", "remain_group_mean", 
                              "target_group_max", "target_group_min", "target_group_sd", "target_group_n0", "target_group_non0", "target_group_sum", 
                              "remain_group_max", "remain_group_min", "remain_group_sd", "remain_group_n0", "remain_group_non0", "remain_group_sum")
              
              return(res)
            })
          ))
        }
        
        # format the results
        if (T) {
          deg_res$FDR <- p.adjust(deg_res$P_value, "fdr")
          deg_res$Target_Cluster <- paste0(cell, "@cluster", c)
          deg_res$ID <- rownames(deg_res)
        }
        
        return(deg_res)
      })
      
      deg_all_clusters <- data.frame(do.call(rbind, deg_all_clusters))
      
      return(deg_all_clusters)
    })
    qsave(results, file = file, nthreads = 6)
  }
}

# format the results
if (T) {
  file <- "data/saved_data/30.deg_analysis_formated_results_for_all_cluster_WNNmatrix_genes.qs"
  if (! file.exists(file)) {
    # simple format
    if (T) {
      deg <- data.frame(do.call(rbind, results))
      head(deg)
      
      deg$cell <- str_split(deg$Target_Cluster, "@", simplify = T)[, 1]
      deg$cluster <- str_split(deg$Target_Cluster, "@", simplify = T)[, 2]
      
      deg$geneID <- str_split(deg$ID, "@", simplify = T)[, 1]
      deg$geneLoc <- str_split(deg$ID, "@", simplify = T)[, 2]
      
      rm(results)
    }
    
    # calculation the difference between target cluster and remain cluster
    if (T) {
      # statistics for mean value and median value
      if (T) {
        head(deg)
        sapply(deg[, 2:5], summary)
        # 使用mean value进行计算
        
        ggdat <- melt(deg[, grep("med|mean|cell", colnames(deg))], id.vars = "cell")
        head(ggdat)
        ggdat$variable <- factor(ggdat$variable, levels = c("target_group_med", "remain_group_med", "target_group_mean", "remain_group_mean"))
        p <- ggplot(ggdat) +
          geom_violin(aes(x=variable, y=value, group=variable, fill=variable), scale = "width", width=1) +
          xlab("Represent Value") +
          ylab("Average Percentage Within A Gene Part") +
          facet_grid(cell~.) +
          theme_bw() +
          theme(axis.title = element_text(size = rel(1.2)),
                axis.text.x = element_text(size = rel(1.2), angle = 90, hjust = 1, vjust = 1),
                axis.text.y = element_text(size = rel(1.2)),
                legend.position = "none", 
                strip.text = element_text(size = rel(1.2)),
                strip.background = element_rect(fill = NA, color=NA))
        ggsave(filename = "results/2.pic/30.med_mean_value_distribution_thp1_cd34.pdf", plot = p, width = 10, height = 6)
      }
      
      # calculate log2fc_mean
      if (T) {
        # calculation
        if (T) {
          # 由于很多数据中target_group_mean数据为0，导致结果中无法存在下调基因，故分子分母各添加一个小数1e-10
          deg$log2fc_mean <- log2((deg$target_group_mean+1e-10)/(deg$remain_group_mean+1e-10))
          head(deg)
          deg$log2fc_mean_state <- ifelse(deg$log2fc_mean == 0, "Notsig", 
                                          ifelse(deg$log2fc_mean > 0, "Up", "Down"))
          
          deg$log2fc_mean_cell_state <- paste0(deg$cell, "_", deg$log2fc_mean_state)
        }
        
        # statistics
        if (T) {
          tapply(deg$log2fc_mean, deg$log2fc_mean_cell_state, summary)
          p <- ggplot(deg) +
            geom_density(aes(x=log2fc_mean, fill=cluster)) +
            facet_grid(cluster~cell, scales = "free") +
            theme_bw() +
            theme(axis.title = element_text(size = rel(1.2)),
                  axis.text = element_text(size = rel(1.2)),
                  legend.position = "none",
                  strip.text = element_text(size = rel(1.2)),
                  strip.background = element_rect(fill = NA, color=NA))
          ggsave(filename = "results/2.pic/30.mean_value_distribution_within_cluster_thp1_cd34.pdf", plot = p, width = 10, height = 6)
        }
        
        # get 0 parts
        if (T) {
          head(deg)
          deg$type <- ifelse(deg$target_group_mean == 0 | deg$remain_group_mean == 0, "num0", "norm")
          
          p <- ggplot(deg) +
            geom_density(aes(x=log2fc_mean, fill=type, group=type)) +
            scale_fill_manual(values=c("red", "grey"),
                              breaks = c("norm", "num0")) +
            facet_grid(type~cell, scales = "free") +
            theme_bw() +
            theme(axis.title = element_text(size = rel(1.2)),
                  axis.text = element_text(size = rel(1.2)),
                  legend.position = "none", 
                  strip.text = element_text(size = rel(1.2)),
                  strip.background = element_rect(fill = NA, color=NA))
          ggsave(filename = "results/2.pic/30.mean_value_distribution_distinct0_thp1_cd34.pdf", plot = p, width = 10, height = 6)
        }
      }
      
      # calculate log2fc_sum
      if (T) {
        # calculation
        if (T) {
          # 由于很多数据中target_group_mean数据为0，导致结果中无法存在下调基因，故分子分母各添加一个小数1e-10
          deg$log2fc_sum <- log2((deg$target_group_sum+1e-10)/(deg$remain_group_sum+1e-10))
          head(deg)
          deg$log2fc_sum_state <- ifelse(deg$log2fc_sum == 0, "Notsig", 
                                         ifelse(deg$log2fc_sum > 0, "Up", "Down"))
          
          deg$log2fc_sum_cell_state <- paste0(deg$cell, "_", deg$log2fc_sum_state)
        }
        
        # statistics
        if (T) {
          tapply(deg$log2fc_sum, deg$log2fc_sum_cell_state, summary)
          ggplot(deg) +
            geom_density(aes(x=log2fc_sum, fill=cluster)) +
            facet_grid(cluster~cell, scales = "free") +
            theme_bw() +
            theme(axis.title = element_text(size = rel(1.2)),
                  axis.text = element_text(size = rel(1.2)),
                  legend.position = "none",
                  strip.text = element_text(size = rel(1.2)))
        }
        
        # get 0 parts
        if (T) {
          head(deg)
          
          ggplot(deg) +
            geom_density(aes(x=log2fc_sum, fill=type, group=type)) +
            scale_fill_manual(values=c("red", "grey"),
                              breaks = c("norm", "num0")) +
            facet_grid(type~cell, scales = "free") +
            theme_bw() +
            theme(axis.title = element_text(size = rel(1.2)),
                  axis.text = element_text(size = rel(1.2)),
                  legend.position = "none", 
                  strip.text = element_text(size = rel(1.2)))
        }
      }
    }
    
    qsave(deg, file = file, nthreads = 6)
  }
  if (file.exists(file)) {
    deg <- qread(file, nthreads = 6)
  }
}

# calculation the difference cutoff between target cluster and remain cluster
if (T) {
  head(deg)
  cutoff <- lapply(cells, function(cell) {
    sapply(c("mean", "sum"), function(type) {
      # test data
      if (F) {
        cell <- cell1
        type <- "mean"
      }
      
      x <- deg[deg$cell == cell & deg$type == "norm", paste0("log2fc_", type)]
      nontop_x <- x[x<=quantile(x, 0.98) & x>=quantile(x, 0.02)]
      xz = (x - mean(nontop_x)) / sd(nontop_x)
      zp = pnorm(xz)
      zp <- ifelse(zp>0.5, 1-zp, zp)
      
      # P<0.05
      filtered <- x[zp<0.05]
      
      cutoff0 <- tapply(filtered, ifelse(filtered>0, "pos", "neg"), function(x) {
        min(abs(x))
      })
      name <- names(cutoff0)
      cutoff0 <- ifelse(name == "pos", cutoff0, -cutoff0)
      names(cutoff0) <- name
      
      return(cutoff0)
    })
  })
  
  # format the results
  names(cutoff) <- cells
  cutoff <- data.frame(do.call(rbind, cutoff))
  cutoff$cell <- rep(cells, each=2)
  cutoff
}

# filter the results with logfc
if (T) {
  head(deg)
  cutoff
  table(deg$log2fc_mean_cell_state, deg$log2fc_sum_cell_state)
  
  # update up and down state
  if (T) {
    for (cell in cells) {
      for (type in c("mean", "sum")) {
        # cell <- cell1
        # type <- "mean"
        cutoff_value_p <- cutoff[cutoff$cell == cell & cutoff[, type]>0, type]
        cutoff_value_n <- cutoff[cutoff$cell == cell & cutoff[, type]<0, type]
        deg[deg$cell == cell, paste0("log2fc_", type, "_state")] <- ifelse(deg[deg$cell == cell, paste0("log2fc_", type)]>cutoff_value_p, 
                                                                           "Up", 
                                                                           ifelse(deg[deg$cell == cell, paste0("log2fc_", type)]<cutoff_value_n, 
                                                                                  "Down", "NotSig"))
        rm(cutoff_value_p, cutoff_value_n)
      }
    }
    
    deg$log2fc_mean_cell_state <- paste0(deg$cell, "_", deg$log2fc_mean_state)
    deg$log2fc_sum_cell_state <- paste0(deg$cell, "_", deg$log2fc_sum_state)
  }
  
  # results similarity
  table(deg$log2fc_mean_cell_state, deg$log2fc_sum_cell_state)
  
  # filter only save Up gene parts
  table(deg$log2fc_mean_state, deg$FDR<0.05)
  torf <- (deg$log2fc_mean_state == "Up") | (deg$log2fc_mean_state == "NotSig" & deg$FDR<0.05 & deg$target_group_mean>deg$remain_group_mean)
  
  # 27820
  sum(torf)
  
  mess <- paste0("There are ", sum(torf), "/", nrow(deg), " (", round(sum(torf)/nrow(deg)*100, 2), "%) gene parts pass the criterions\n")
  cat(mess)
  
  deg <- deg[torf, ]
}

# save the data
file <- "data/saved_data/30.filtered_deg_analysis_for_all_cluster_WNNmatrix_genes.qs"
qsave(deg, file = file, nthreads = 6)


# 在知道每个细胞中chromatin states的分群结果后，下一步便是功能注释
# 这里主要从转录层面对分群结果进行注释
# 对基因body染色质状态的percentage矩阵进行差异分析，将相同cluster的作为一组，与剩余所有染色质状态进行比较差异
# line1-24: 读入不同细胞中的染色质状态分群结果
# line25-45: 读入不同细胞中所有基因part中染色质状态的占比数据，同时对数据进行过滤：
#            基因所占bin数一定要>=body_bin_num，以保证每个基因part至少独立享受1个bin
# line46-193: 进行差异检查分析：
#                 （1）因为每个基因part中所有的比例之和总为1，故比较的原则是看在基因每个part部分中，不同cluster中染色质状态所占比例之和的大小
#                      类似一个公司中有不同党派，选择最终人数最多的党派时，直接计算总和即可
#                 （2）差异检验本质上是对均值的检验，故在此其实不适合使用。但是也可以进行计算，计算方法使用Wilcoxon
#                 （3）目标差异基因挑选前，需要选择合适的值来当作代表，传统差异分析是以均值进行代表，这里由于不确定何种值适合，故进行一系列探索：
#                      1) 分别计算不同cluster在每个基因part所占比例的均值和中位值，发现中位值多为0，故选择使用mean值进行计算差异
#                      2) 由于很多数据中target_group_mean数据为0，导致log2(deg$target_group_mean/deg$remain_group_mean)几乎无下调基因, 故分子分母各添加一个小数1e-10
#                      3) 一系列计算指标的含义：
#                        - log2fc_sum： 以cluster中所有state占比总和计算
#                        - log2fc_mean： 以cluster中所有state占比均值计算
#                        - log2fc_{sum,mean}： cutoff的计算，根据数据分布而定：
#                      4) 绘制log2fc_{sum,mean}值的分布
#                      5) log2fc_{sum,mean}分布大致符合正态分布
#                      6) 将log2fc_{sum,mean}值进行zscale，然后计算P值，过滤标准选择P<0.05
#                      7) 考虑到实际意义中mean值可以理解为矫正了每个cluster中state数目的结果，故优先使用log2fc_mean的结果进行后续分析，仅保留上调基因part作为差异，具体原因和单细胞仅使用pos结果一致
# line194-333: 数据格式化，同时对log2fc_{sum,mean}数值的分布模式进行探索
# line334-368: 分别将log2fc_{sum,mean}值转为正态分布，找出P<0.05时的cutoff
# line369-425: 根据cutoff值，将基因分为上下调基因，并统计相应结果

