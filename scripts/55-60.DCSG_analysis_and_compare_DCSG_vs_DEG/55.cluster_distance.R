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

# read distance dat
if (T) {
  file <- "results/1.tab/28.distance_for_states.csv"
  dists <- read.table(file, header = T, sep = ",", fill = T, comment.char = "")
  rownames(dists) <- dists$state
  
  dists <- dists[, dists$state]
}

# read cluster data
if (T) {
  file <- "data/saved_data/29.common_cluster_name_for_merged_2cells_CS_clusters.qs"
  clusters <- qread(file, nthreads=6)
  clusters <- clusters[order(clusters$seurat_clusters), ]
}

# format the data
if (T) {
  dists <- dists[clusters$orig.ident, ]
  dists <- dists[, row.names(dists)]
}

# normalization: max-min-normalization
if (T) {
  values <- c(as.matrix(dists))
  max <- max(values)
  min <- min(values)
  
  dists <- (dists-min)/(max-min)
  dists <- as.matrix(dists)
  
  rm(values, max, min)
}

# color setting
if (T) {
  library(RColorBrewer)
  colors <- brewer.pal(9, "Set1") 
  
  colors <- colors[1:2]
  col_fun <- colorRamp2(c(0, 0.5, 1), c(colors[1], "white", colors[2]))
}

# all state distance
if (T) {
  # heatmap plot
  if (T) {
    p <- Heatmap(dists, name = "Distance", col = col_fun,
                 border_gp = gpar(col = "black"), rect_gp = gpar(col = "black"),
                 show_heatmap_legend = TRUE, heatmap_legend_param = list(border = T, title_position = "lefttop", legend_direction = "horizontal"),
                 column_title = "c%s", row_title = "c%s",
                 row_title_gp = gpar(fontsize = 13.2), column_title_gp = gpar(fontsize = 13.2),
                 row_split = clusters$seurat_clusters, column_split = clusters$seurat_clusters, 
                 row_names_max_width = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), row_names_rot = 0,
                 column_names_max_height = unit(6, "cm"), column_names_gp = gpar(fontsize = 10), column_names_rot = 60,
                 cluster_rows = T, cluster_columns = T, 
                 clustering_distance_columns = "euclidean", clustering_distance_rows = "euclidean", 
                 clustering_method_columns = "ward.D2", clustering_method_rows = "ward.D2")
    draw(p, heatmap_legend_side = "bottom")
  }
  
  # save the results
  if (T) {
    p <- ggplotify::as.ggplot(p)
    
    file <- "results/2.pic/55.cluster_distance_detailed_heatmap.pdf"
    ggsave(filename = file, plot = p, width = 8, height = 6)
  }
  
  # save the dat
  if (T) {
    file <- "data/saved_data/55.cs_distance.qs"
    qsave(dists, file = file, nthreads = 6)
  }
  
  rm(p, file)
}

# all cluster distance
if (T) {
  # calculate the distance for each cluster
  if (T) {
    all_clusters <- sort(unique(clusters$seurat_clusters))
    distance <- lapply(all_clusters, function(c) {
      # c <- 1
      
      dis <- sapply(all_clusters, function(r_c) {
        # r_c <- all_clusters[1]
        r_states <- clusters$orig.ident[clusters$seurat_clusters == r_c]
        c_states <- clusters$orig.ident[clusters$seurat_clusters == c]
        
        mean(dists[c_states, r_states], na.rm = T)
      })
      
      dis <- data.frame(target_c = c, 
                        compared_c = all_clusters, 
                        distance = dis)
      rownames(dis) <- 1:nrow(dis)
      
      return(dis)
    })
    distance <- do.call(rbind, distance)
  }
  
  # long2wide
  if (T) {
    head(distance)
    distance <- dcast(distance, target_c~compared_c, value.var = "distance")
    rownames(distance) <- paste0("C", distance[, 1])
    
    distance <- distance[, -1]
    colnames(distance) <- paste0("C", colnames(distance))
    
    distance <- as.matrix(distance)
  }
  
  # heatmap plot
  if (T) {
    p <- Heatmap(distance, name = "Distance", col = col_fun,
                 border_gp = gpar(col = "black"), rect_gp = gpar(col = "black"),
                 show_heatmap_legend = TRUE, heatmap_legend_param = list(border = T, title_position = "lefttop", legend_direction = "horizontal"),
                 row_title_gp = gpar(fontsize = 13.2), column_title_gp = gpar(fontsize = 13.2),
                 row_names_max_width = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), row_names_rot = 0,
                 column_names_max_height = unit(6, "cm"), column_names_gp = gpar(fontsize = 10), column_names_rot = 60,
                 cluster_rows = T, cluster_columns = T, 
                 clustering_distance_columns = "euclidean", clustering_distance_rows = "euclidean", 
                 clustering_method_columns = "ward.D2", clustering_method_rows = "ward.D2", 
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(round(distance[i, j], 2), x, y, gp = gpar(fontsize = 10))
                 })
    draw(p, heatmap_legend_side = "bottom")
  }
  
  # save the results
  if (T) {
    p <- ggplotify::as.ggplot(p)
    
    file <- "results/2.pic/55.cluster_distance_merged_heatmap.pdf"
    ggsave(filename = file, plot = p, width = 8, height = 6)
  }
  
  # save the dat
  if (T) {
    file <- "data/saved_data/55.cluster_distance_merged_heatmap.qs"
    qsave(distance, file = file, nthreads = 6)
  }
}
