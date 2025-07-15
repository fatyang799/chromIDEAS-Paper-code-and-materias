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
  # separated data
  if (T) {
    file <- "data/saved_data/25.common_cluster_name_for_2_cell_CS_clusters.qs"
    cluster <- qread(file, nthreads = 6)
    cluster$orig.ident <- paste0("S", cluster$orig.ident)
    
    cluster1 <- cluster[cluster$cell == cell1, ]
    cluster2 <- cluster[cluster$cell == cell2, ]
  }
  
  # merged data
  if (T) {
    file <- "results/1.tab/28.states_group_info.csv"
    merged <- read.table(file, header = T, sep = ",", fill = T, comment.char = "")
    colnames(merged)
    merged <- merged[, c("orig.ident", "seurat_clusters")]
    merged$seurat_clusters <- merged$seurat_clusters+1
    merged$cell <- "merged"
  }
  
  rm(file, cluster)
}

# get corresponding relationship
if (T) {
  # 0->1
  # 1->2
  # 2->5
  # 3->4
  # 4->3
  
  # simple statistics
  if (T) {
    table(cluster1$seurat_clusters, merged$seurat_clusters)
    table(cluster2$seurat_clusters, merged$seurat_clusters)
  }
  
  # cell1/cell2 as the reference: cell1 vs merged
  if (T) {
    # calculate the relationship base on the maximum overlap
    if (identical(cluster1$orig.ident, merged$orig.ident)) {
      relation1 <- table(cluster1$seurat_clusters, merged$seurat_clusters)
      new_relation1 <- sapply(rownames(relation1), function(c1) {
        colnames(relation1)[which.max(relation1[c1, ])]
      })
      
      new_relation1 <- data.frame(
        cell1 = paste0(cell1, "_", rownames(relation1)),
        merged = paste0("merged_", new_relation1)
      )
      
      new_relation1
    } else {
      stop(paste0("The state of between merged and ", cell1, " is not identical"))
    }
    
    # rename to make corresponding relationship more human friendly
    if (T) {
      new_relation1
      
      merged$seurat_clusters
      for (i in 1:nrow(merged)) {
        # i <- 1
        old_c <- merged$seurat_clusters[i]
        new_c <- new_relation1[new_relation1[, 2] == paste0("merged_", old_c), 1]
        new_c <- as.numeric(str_extract(new_c, "[0-9]{1,2}$"))
        
        merged$seurat_clusters_cell1[i] <- new_c
      }
      rm(old_c, new_c, i)
    }
  }
  
  # cell1/cell2 as the reference: cell2 vs merged
  if (T) {
    # calculate the relationship base on the maximum overlap
    if (identical(cluster2$orig.ident, merged$orig.ident)) {
      relation2 <- table(cluster2$seurat_clusters, merged$seurat_clusters)
      new_relation2 <- sapply(rownames(relation2), function(c1) {
        colnames(relation2)[which.max(relation2[c1, ])]
      })
      
      new_relation2 <- data.frame(
        cell2 = paste0(cell2, "_", rownames(relation2)),
        merged = paste0("merged_", new_relation2)
      )
      
      new_relation2
    } else {
      stop("The state of between merged and thp1 is not identical")
    }
    
    # rename to make corresponding relationship more human friendly
    if (T) {
      new_relation2
      
      merged$seurat_clusters
      for (i in 1:nrow(merged)) {
        # i <- 1
        old_c <- merged$seurat_clusters[i]
        new_c <- new_relation2[new_relation2[, 2] == paste0("merged_", old_c), 1]
        new_c <- as.numeric(str_extract(new_c, "[0-9]{1,2}$"))
        
        merged$seurat_clusters_cell2[i] <- new_c
      }
      rm(old_c, new_c, i)
    }
  }
  
  # rename the cluster
  if (identical(merged$seurat_clusters_cell1, merged$seurat_clusters_cell2)) {
    merged$seurat_clusters <- merged$seurat_clusters_cell1
  } else {
    stop("The relationship is not identical between thp1 and cd34 compared with merged data")
  }
  
  rm(relation1, relation2, new_relation1, new_relation2)
}

# state cluster similarity sankey
if (T) {
  library(networkD3)
  
  # link data prepare
  if (T) {
    # statistics
    if (T) {
      cluster1 <- cluster1[match(merged$orig.ident, cluster1$orig.ident), ]
      cluster2 <- cluster2[match(merged$orig.ident, cluster2$orig.ident), ]
      dat <- as.data.frame.table(table(paste(cluster1$seurat_clusters, 
                                             cluster2$seurat_clusters, 
                                             merged$seurat_clusters, sep = "_")))
      dat[, cell1] <- paste0(cell1, "_", str_split(dat$Var1, "_", simplify = T)[, 1])
      dat[, cell2] <- paste0(cell2, "_", str_split(dat$Var1, "_", simplify = T)[, 2])
      dat$merged <- paste0("merged_", str_split(dat$Var1, "_", simplify = T)[, 3])
    }
    
    # format the results
    if (T) {
      link1 <- aggregate(x = dat$Freq, by=list(dat[, cell1], dat$merged), sum)
      link2 <- aggregate(x = dat$Freq, by=list(dat$merged, dat[, cell2]), sum)
      df <- rbind(link1, link2)
      
      colnames(df) <- c('source',"target","value")
    }
    
    # group setting
    if (T) {
      df$group <- ifelse(str_split(df$source, "_", simplify = T)[, 2] == str_split(df$target, "_", simplify = T)[, 2], 
                         "main", "other")
      
      df$group <- factor(df$group)
    }
  }
  
  # node data prepare
  if (T) {
    df.nodes <- data.frame(
      name=unique(c(df$source, df$target))
    )
    df.nodes$group <- factor(df.nodes$name)
    
    df$IDsource <- match(df$source, df.nodes$name)-1
    df$IDtarget <- match(df$target, df.nodes$name)-1
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
                       paste0('"', paste(rep(c("#F3722C", "#F8961E", "#43AA8B", "#4D908E"), 3), collapse = '", "'), '"'), 
                       '])')
  }
  
  # sankey plot
  if (T) {
    Sankey.p <- sankeyNetwork(Links = df, Nodes = df.nodes,
                              Source = "IDsource", Target = "IDtarget",
                              Value = "value", NodeID = "name", 
                              LinkGroup="group", NodeGroup="group",
                              colourScale=my_color, fontFamily="arial", 
                              nodeWidth=30, nodePadding=15, sinksRight=T,
                              fontSize=20)
    Sankey.p
  }
  
  # save the figure
  if (T) {
    htmlwidgets::saveWidget(Sankey.p, 
                            file=paste0("results/2.pic/29.Sankey_between_merged_and_separated_CS_clusters.html"))
    webshot::webshot(paste0("results/2.pic/29.Sankey_between_merged_and_separated_CS_clusters.html"), 
                     paste0("results/2.pic/29.Sankey_between_merged_and_separated_CS_clusters.pdf"))
  }
  
  rm(dat, df, df.nodes, link1, link2, Sankey.p, my_color)
}

# plot emission figure based on the cluster data
if (T) {
  emission <- emission[merged$orig.ident, ]
  
  # color setting
  if (T) {
    cluster <- sort(unique(merged$seurat_clusters))
    col_fun <- colorRamp2(c(1, length(cluster)), c("blue", "red"))
    cluster <- col_fun(1:length(cluster))
    names(cluster) <- sort(unique(merged$seurat_clusters))
    
    col_fun = colorRamp2(c(0, ceiling(quantile(c(as.matrix(emission)), 0.95))),
                         c("#FFFFFF", "#3953A4"))
  }
  
  # ordering
  if (T) {
    emission <- emission[, c("ATAC", "H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3", "H3K9me3", "H3K79me2")]
  }
  
  # heatmap plot
  if (T) {
    merged$seurat_clusters <- factor(merged$seurat_clusters, levels = sort(unique(merged$seurat_clusters)))
    p <- Heatmap(as.matrix(emission), name = "Probablity", col = col_fun, 
                 border_gp = gpar(col = "black"), rect_gp = gpar(col = "black"),
                 show_heatmap_legend = TRUE, heatmap_legend_param = list(border = T),
                 row_split = merged$seurat_clusters, cluster_row_slices = FALSE, 
                 row_names_max_width = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), row_names_rot = 0,
                 column_names_max_height = unit(6, "cm"), column_names_gp = gpar(fontsize = 10), column_names_rot = 60,
                 cluster_rows = T, cluster_columns = F, show_row_dend = F, show_column_dend = F, 
                 clustering_distance_columns = "euclidean", clustering_distance_rows = "euclidean", 
                 clustering_method_columns = "ward.D2", clustering_method_rows = "ward.D2", 
                 right_annotation = rowAnnotation(clusters = factor(merged$seurat_clusters), 
                                                  col = list(clusters = cluster)))
    draw(p, show_heatmap_legend = T)
  }
  
  # save the results
  if (T) {
    file <- "results/2.pic/29.emission_heatmap_merged_2cells_CS_cluster.pdf"
    
    p <- ggplotify::as.ggplot(p)
    ggsave(filename = file, plot = p, width = 3, height = 10)
  }
  
  rm(cluster, col_fun, p, file)
}

# replot the histone signal based on the clusters
if (T) {
  # prepare the data
  if (T) {
    head(merged)
    
    ggdat <- cbind(merged, emission[merged$orig.ident, ])
    
    ggdat <- melt(ggdat, id.vars = colnames(ggdat)[1:5], variable.name = "mk")
    ggdat$seurat_clusters <- factor(ggdat$seurat_clusters, levels = sort(unique(ggdat$seurat_clusters)))
    ggdat$mk <- factor(ggdat$mk, levels = c("ATAC", "H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3", "H3K9me3", "H3K79me2"))
  }
  
  # plot
  if (T) {
    head(ggdat)
    
    p <- ggplot(ggdat) +
      geom_violin(aes(x=seurat_clusters, y=value, group=seurat_clusters, fill=seurat_clusters), scale="width") +
      geom_jitter(aes(x=seurat_clusters, y=value)) +
      ylab("Histone Signal") +
      ggtitle("merged_cells") +
      facet_wrap(~mk, scales = "free") +
      cowplot::theme_cowplot() +
      theme(axis.title = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1)),
            legend.position = "none",
            strip.background = element_rect(fill = NA), 
            strip.text = element_text(size = rel(1), face="bold"))
    
    ggsave("results/2.pic/29.states_cluster_merged_2cells_emission_features_vlnplot.pdf", plot = p, width = 10, height = 10)
  }
  
  # merge all 3 situation plot
  if (T) {
    # get cell specific data
    if (T) {
      dat <- qread("data/saved_data/25.merged_emission_table_violin_mat.qs", nthreads = 6)
    }
    
    # merge the data
    if (T) {
      ggdat <- ggdat[, colnames(dat)]
      ggdat <- rbind(ggdat, dat)
    }
    
    head(ggdat)
    
    p <- ggplot(ggdat) +
      geom_violin(aes(x=mk, y=value, group=mk, fill=mk), scale="width") +
      geom_jitter(aes(x=mk, y=value)) +
      ylab("Histone Signal") +
      facet_grid(cell~seurat_clusters, scales = "free_x") +
      cowplot::theme_cowplot() +
      theme(axis.title = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1)),
            legend.position = "none",
            panel.border = element_rect(color="black"), 
            strip.background = element_rect(fill = NA), 
            strip.text = element_text(size = rel(1), face="bold"))
    
    ggsave(filename = "results/2.pic/29.states_cluster_merged_3_situation_emission_features_vlnplot.pdf", plot = p, 
           width = 15, height = 8)
  }
}

# save the data
if (T) {
  file <- "data/saved_data/29.common_cluster_name_for_merged_2cells_CS_clusters.qs"
  qsave(merged, file, nthreads = 6)
}





