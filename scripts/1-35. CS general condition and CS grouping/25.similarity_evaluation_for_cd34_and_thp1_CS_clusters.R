# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(qs)
  library(reshape2)
  library(stringr)
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(circlize))
}

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
body_bin_num <- 10
rna_levels <- 10
length_leveles <- 10
emission_file <- "data/raw_data/1.emission/chromIDEAS.emission.txt"

# read the emission
if (T) {
  emission <- read.table(emission_file, header = T, sep = "\t")
  rownames(emission) <- emission[, 1]
  emission <- emission[, -c(1:2)]
}

# read the cluster data
if (T) {
  cluster1 <- read.table(paste0("results/1.tab/24.", cell1, "_states_group_info.csv"), header = T, sep = ",")
  cluster2 <- read.table(paste0("results/1.tab/24.", cell2, "_states_group_info.csv"), header = T, sep = ",")
  
  if (identical(cluster1$orig.ident, cluster2$orig.ident)) {
    table(cluster1$seurat_clusters, cluster2$seurat_clusters)
  }
}

# get corresponding relationship
if (T) {
  # cd34:
  #     0->1
  #     1->2
  #     2->3
  #     3->4
  #     4->5
  # thp1:
  #     0->1
  #     1->2
  #     2->5
  #     3->3
  #     4->4
  cluster1$seurat_clusters <- cluster1$seurat_clusters+1
  cluster2$seurat_clusters <- cluster2$seurat_clusters+1
  
  relation <- as.data.frame.matrix(table(cluster1$seurat_clusters, cluster2$seurat_clusters))
  n <- ifelse(nrow(relation)>=ncol(relation), 1, 2)
  
  # number of cell1 cluster is more
  if (n == 1) {
    # calculate the relationship base on the maximum overlap
    if (T) {
      new_relation <- sapply(rownames(relation), function(c1) {
        colnames(relation)[which.max(relation[c1, ])]
      })
      
      new_relation <- data.frame(
        cell1 = paste0(cell1, "_", rownames(relation)),
        cell2 = paste0(cell2, "_", new_relation)
      )
    }
    
    # rename to make corresponding relationship more human friendly
    if (T) {
      new_relation
      dup_c <- new_relation[, 3-n]
      dup_c <- unique(dup_c[duplicated(dup_c)])
      dup_c <- new_relation[new_relation[, 3-n] %in% dup_c, n]
      
      cluster1$seurat_clusters
      for (i in 1:nrow(cluster1)) {
        # i <- 16
        old_c <- cluster1$seurat_clusters[i]
        new_c <- new_relation[new_relation[, 1] == paste0(cell1, "_", old_c), 2]
        new_c <- as.numeric(str_extract(new_c, "[0-9]{1,2}$"))
        
        cluster1$seurat_clusters[i] <- ifelse(paste0(cell1, "_", old_c) %in% dup_c, as.numeric(paste0(new_c, old_c)), new_c)
      }
    }
  }
  
  # number of cell2 cluster is more
  if (n == 2) {
    # calculate the relationship base on the maximum overlap
    if (T) {
      new_relation <- sapply(colnames(relation), function(c2) {
        rownames(relation)[which.max(relation[, c2])]
      })
      
      new_relation <- data.frame(
        cell1 = paste0(cell1, "_", new_relation),
        cell2 = paste0(cell2, "_", colnames(relation))
      )
    }
    
    # rename to make corresponding relationship more human friendly
    if (T) {
      new_relation
      dup_c <- new_relation[, 3-n]
      dup_c <- unique(dup_c[duplicated(dup_c)])
      dup_c <- new_relation[new_relation[, 3-n] %in% dup_c, n]
      
      cluster2$seurat_clusters
      for (i in 1:nrow(cluster2)) {
        # i <- 7
        old_c <- cluster2$seurat_clusters[i]
        new_c <- new_relation[new_relation[, 2] == paste0("cd34_", old_c), 1]
        new_c <- as.numeric(str_extract(new_c, "[0-9]{1,2}$"))
        
        cluster2$seurat_clusters[i] <- ifelse(paste0("cd34_", old_c) %in% dup_c, as.numeric(paste0(new_c, old_c)), new_c)
      }
    }
  }
  
  rm(old_c, new_c, dup_c, i, n, new_relation, relation)
}

# state cluster similarity heatmap
if (T) {
  # format the data
  if (T) {
    cluster2 <- cluster2[match(cluster1$orig.ident, cluster2$orig.ident), ]
  }
  
  # statistics
  if (T) {
    if (identical(cluster1$orig.ident, cluster2$orig.ident)) {
      ggdat <- as.data.frame.matrix(table(cluster1$seurat_clusters, cluster2$seurat_clusters))
      rownames(ggdat) <- paste0(cell1, "_cluster", rownames(ggdat))
      colnames(ggdat) <- paste0(cell2, "_cluster", colnames(ggdat))
      ggdat <- as.matrix(ggdat)
    } else {
      mess <- "The order of data is not identical, please check\n"
      stop(mess)
    }
  }
  
  # color setting
  if (T) {
    library(RColorBrewer)
    colors <- brewer.pal(9, "Set1") 
    
    colors <- colors[1:2]
  }
  
  # heatmap
  if (T) {
    p <- Heatmap(ggdat, col = c("white", colors[1]), 
                 border_gp = gpar(col = "black"), rect_gp = gpar(col = "black"),
                 show_heatmap_legend = F, 
                 row_names_max_width = unit(6, "cm"), row_names_gp = gpar(fontsize = 15), row_names_centered = T,
                 column_names_max_height = unit(6, "cm"), column_names_gp = gpar(fontsize = 15), column_names_rot = 0, column_names_centered = T,
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(ggdat[i, j], x, y, gp = gpar(fontsize = 13))
                 }, 
                 cluster_rows = F, cluster_columns = F, show_row_dend = F, show_column_dend = F)
  }
  
  # save the figure
  if (T) {
    p <- ggplotify::as.ggplot(p)
    
    file <- paste0("results/2.pic/25.similarity_between_", cell2, "_", cell1, "_CS_clusters_heatmap.pdf")
    ggsave(filename = file, plot = p, width = 8, height = 6)
  }
  
  rm(p, ggdat)
}

# state cluster similarity sankey
if (T) {
  library(networkD3)
  
  # link data prepare
  if (T) {
    df <- as.data.frame.table(table(cluster1$seurat_clusters, cluster2$seurat_clusters))
    colnames(df) <- c(cell1, cell2, "value")
    
    df[, cell1] <- paste0(cell1, "_cluster", df[, cell1])
    df[, cell2] <- paste0(cell2, "_cluster", df[, cell2])
    
    df$group1 <- str_split(df[, cell1], "_", simplify = T)[, 2]
    df$group2 <- str_split(df[, cell2], "_", simplify = T)[, 2]
    
    df$group1 <- str_extract(df$group1, "cluster[0-9]")
    df$group2 <- str_extract(df$group2, "cluster[0-9]")
    
    df$group <- ifelse(df$group1 == df$group2, "main", "other")
    df$group <- factor(df$group)
  }
  
  # node data prepare
  if (T) {
    df.nodes <- data.frame(
      name=unique(c(df[, cell1], df[, cell2]))
    )
    df.nodes$group <- factor(df.nodes$name)
    
    df$IDsource <- match(df[, cell1], df.nodes$name)-1
    df$IDtarget <- match(df[, cell2], df.nodes$name)-1  
  }
  
  # color setting
  if (T) {
    my_color <- paste0('d3.scaleOrdinal() .domain([', 
                       paste0('"', paste(levels(df$group), collapse = '", "'), '"'), 
                       ', ', 
                       paste0('"', paste(levels(df.nodes$group), collapse = '", "'), '"'), 
                       ']) .range([', 
                       paste0('"', paste(c("#577590", "#90BE6D"), collapse = '", "'), '"'), 
                       ', ', 
                       paste0('"', paste(c("#F3722C", "#F8961E", "#43AA8B", "#4D908E"), collapse = '", "'), '"'), 
                       '])')
  }
  
  # sankey plot
  if (T) {
    Sankey.p <- sankeyNetwork(Links = df, Nodes = df.nodes,
                              Source = "IDsource", Target = "IDtarget",
                              Value = "value", NodeID = "name", 
                              LinkGroup="group", NodeGroup="group",
                              colourScale=my_color, fontFamily="arial", 
                              nodeWidth=20, nodePadding=15, sinksRight=T,
                              fontSize=20)
    Sankey.p
  }
  
  # save the figure
  if (T) {
    htmlwidgets::saveWidget(Sankey.p, 
                            file=paste0("results/2.pic/25.Sankey_between_", cell2, "_", cell1, "_CS_clusters.html"))
    webshot::webshot(paste0("results/2.pic/25.Sankey_between_", cell2, "_", cell1, "_CS_clusters.html"), 
                     paste0("results/2.pic/25.Sankey_between_", cell2, "_", cell1, "_CS_clusters.pdf"))
  }
  
  rm(df, df.nodes, my_color, Sankey.p)
}

# plot the weight for each CS
if (T) {
  # merge the data
  if (T) {
    cluster1$cell <- cell1
    cluster2$cell <- cell2
    dat <- rbind(cluster1, cluster2)
  }
  
  # format the data
  if (T) {
    ggdat <- melt(dat, id.vars = c("orig.ident", "cell", "seurat_clusters"), measure.vars = c("gene.weight", "emission.weight"))
    
    if (all(ggdat$value>0)) {
      ggdat$value <- ifelse(ggdat$variable == "gene.weight", ggdat$value, 0-ggdat$value)  
    }
    ggdat$orig.ident <- gsub("S", "", ggdat$orig.ident)
    ggdat$orig.ident <- factor(ggdat$orig.ident, levels = seq(min(as.numeric(ggdat$orig.ident)), max(as.numeric(ggdat$orig.ident))))
  }
  
  # annotation
  if (T) {
    anno <- data.frame(x=c(1, 1), y=c(0.5, -0.5), label=c("gene.weight", "emission.weight"))  
  }
  
  # ggplot
  if (T) {
    head(ggdat)
    p <- ggplot(ggdat) +
      geom_hline(yintercept = c(0.5, -0.5), colour="black", linetype=2, alpha=0.4) +
      geom_hline(yintercept = 0, colour="black", linetype=1, alpha=0.4) +
      geom_segment(aes(x=orig.ident, y=value, xend=orig.ident, yend=0, colour = orig.ident)) +
      geom_text(data = anno, aes(x=x, y=y, label=label, angle=90, alpha=0.5)) +
      geom_point(aes(x=orig.ident, y=value, colour = orig.ident), size=4) +
      scale_y_continuous(breaks = seq(-1,1,0.5), limits=c(-1, 1), labels = abs(seq(-1,1,0.5))) +
      theme_bw() +
      xlab("State") +
      ylab("Weight") +
      facet_grid(cell~.) +
      theme(axis.title = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1.2)),
            panel.border = element_rect(color="black"), 
            strip.background = element_rect(fill=NA, color=NA), 
            legend.position = "none",
            strip.text = element_text(size = rel(1.2)))
  }
  ggsave("results/2.pic/25.WNN_emission_gene_weight.pdf", plot = p, width = 12, height = 6)
  
  # ggplot
  if (T) {
    head(ggdat)
    p <- ggplot(ggdat) +
      geom_hline(yintercept = c(0.5, -0.5), colour="black", linetype=2, alpha=0.4) +
      geom_hline(yintercept = 0, colour="black", linetype=1, alpha=0.4) +
      geom_segment(aes(x=orig.ident, y=value, xend=orig.ident, yend=0, colour = orig.ident)) +
      geom_text(data = anno, aes(x=x, y=y, label=label, angle=90, alpha=0.5)) +
      geom_point(aes(x=orig.ident, y=value, colour = orig.ident), size=4) +
      scale_y_continuous(breaks = seq(-1,1,0.5), limits=c(-1, 1), labels = abs(seq(-1,1,0.5))) +
      theme_bw() +
      xlab("State") +
      ylab("Weight") +
      facet_grid(cell~seurat_clusters, scales = "free_x") +
      theme(axis.title = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1.2)),
            panel.border = element_rect(color="black"), 
            strip.background = element_rect(fill=NA, color=NA), 
            legend.position = "none",
            strip.text = element_text(size = rel(1.2)))
  }
  ggsave("results/2.pic/25.WNN_emission_gene_weight_cell_specific_cluster_faceted.pdf", plot = p, width = 12, height = 6)
  
  rm(dat, ggdat, anno)
}

# plot emission figure based on the cluster data
if (T) {
  emission <- emission[cluster1$orig.ident, ]
  
  # cell specific cluster data
  for (cell in c(cell1, cell2)) {
    # cell <- cell1
    # get cell specific cluster data
    if (T) {
      if (cell == cell1) {
        dat <- cluster1
      }
      if (cell == cell2) {
        dat <- cluster2
      }
    }
    
    # color setting
    if (T) {
      cluster <- sort(unique(dat$seurat_clusters))
      col_fun <- colorRamp2(c(1, length(cluster)), c("blue", "red"))
      cluster <- col_fun(1:length(cluster))
      names(cluster) <- sort(unique(dat$seurat_clusters))
      
      col_fun = colorRamp2(c(0, ceiling(quantile(c(as.matrix(emission)), 0.95))),
                           c("#FFFFFF", "#3953A4"))
    }
    
    # ordering
    if (T) {
      emission <- emission[, c("ATAC", "H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3", "H3K9me3", "H3K79me2")]
    }
    
    # heatmap plot
    if (T) {
      dat$seurat_clusters <- factor(dat$seurat_clusters, levels = sort(unique(dat$seurat_clusters)))
      p <- Heatmap(as.matrix(emission), name = "Probablity", col = col_fun, 
                   border_gp = gpar(col = "black"), rect_gp = gpar(col = "black"),
                   show_heatmap_legend = TRUE, heatmap_legend_param = list(border = T),
                   row_split = dat$seurat_clusters, cluster_row_slices = FALSE, 
                   row_names_max_width = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), row_names_rot = 0,
                   column_names_max_height = unit(6, "cm"), column_names_gp = gpar(fontsize = 10), column_names_rot = 60,
                   cluster_rows = T, cluster_columns = F, show_row_dend = F, show_column_dend = F, 
                   clustering_distance_columns = "euclidean", clustering_distance_rows = "euclidean", 
                   clustering_method_columns = "ward.D2", clustering_method_rows = "ward.D2", 
                   right_annotation = rowAnnotation(clusters = factor(dat$seurat_clusters), 
                                                    col = list(clusters = cluster)))
    }
    
    # save the results
    if (T) {
      file <- paste0("results/2.pic/25.emission_heatmap_", cell, "_CS_cluster.pdf")
      
      p <- ggplotify::as.ggplot(p)
      ggsave(filename = file, plot = p, width = 3, height = 10)
    }
  }
  
  rm(cell, dat, cluster, col_fun, p, file)
}

# replot the histone signal based on the clusters
if (T) {
  # get cluster data
  if (T) {
    dat <- rbind(cluster1, cluster2)
    dat <- dat[, c("orig.ident", "seurat_clusters", "cell")]
    dat$orig.ident <- as.numeric(gsub("S", "", dat$orig.ident))
    dat$orig.ident <- factor(dat$orig.ident, levels = min(dat$orig.ident):max(dat$orig.ident))
  }
  
  # prepare the data
  if (T) {
    head(dat)
    
    dat <- lapply(c(cell1, cell2), function(cell) {
      # cell <- cell1
      subdat <- dat[dat$cell == cell, ]
      subdat <- cbind(subdat, emission[paste0("S", as.character(subdat$orig.ident)), ])
      
      return(subdat)
    })
    dat <- do.call(rbind, dat)
    
    dat <- melt(dat, id.vars = colnames(dat)[1:3])
    colnames(dat)[4] <- "mk"
    
    dat$seurat_clusters <- factor(dat$seurat_clusters, levels = sort(unique(dat$seurat_clusters)))
    dat$mk <- factor(dat$mk, levels = c("ATAC", "H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3", "H3K9me3", "H3K79me2"))
  }
  
  # plot
  if (T) {
    head(dat)
    
    # define function
    if (T) {
      plot_violin <- function(dat, cell) {
        # cell <- cell1
        cell_dat <- dat[dat$cell == cell, ]
        
        p <- ggplot(cell_dat) +
          geom_violin(aes(x=seurat_clusters, y=value, group=seurat_clusters, fill=seurat_clusters), scale="width") +
          geom_jitter(aes(x=seurat_clusters, y=value)) +
          ylab("Histone Signal") +
          ggtitle(cell) +
          facet_wrap(~mk, scales = "free") +
          cowplot::theme_cowplot() +
          theme(axis.title = element_text(size = rel(1.2)),
                axis.text = element_text(size = rel(1)),
                legend.position = "none",
                strip.background = element_rect(fill = NA), 
                strip.text = element_text(size = rel(1), face="bold"))
        
        return(p)
      }
    }
    
    # plot cell specific plot
    if (T) {
      plot_violin(dat, cell1)
      ggsave(paste0("results/2.pic/25.states_cluster_", cell1, "_emission_features_vlnplot.pdf"), 
             width = 10, height = 10)
      
      plot_violin(dat, cell2)
      ggsave(paste0("results/2.pic/25.states_cluster_", cell2, "_emission_features_vlnplot.pdf"), 
             width = 10, height = 10)
    }
    
    qsave(dat, "data/saved_data/25.merged_emission_table_violin_mat.qs", nthreads = 6)
  }
}

# save the data
if (T) {
  dat <- dat[, grep("orig.ident|seurat_clusters|cell", colnames(dat))]
  dat <- dat[! duplicated(dat), ]
  
  file <- "data/saved_data/25.common_cluster_name_for_2_cell_CS_clusters.qs"
  qsave(dat, file)
}



# line1-36: 读入emission数据以及染色质状态分群结果
# line37-113: 由于前面cd34和thp1细胞中染色质状态分群是独立进行的，因此他们的分群名称不一致。
#             为了找到2个细胞中分群名称的对应关系，这里使用最大重叠原则来判断2个细胞之间分群的对应关系
#             并将对应关系进行修改以保证后续分析更符合人类阅读习惯，即分群数字从1开始，而非0
# line114-160: cd34和thp1细胞间分群的保守性，chromatin states转换热图
# line161-229: cd34和thp1细胞间分群的保守性，sankey plot
# line230-296: 对WNN算法中权重的含义进行可视化：棒棒糖图
# line297-353: 按照聚类结果可视化emission热图：热图
# line354-422: 按照聚类结果可视化每个cluster的组蛋白信号强度值：小提琴图
