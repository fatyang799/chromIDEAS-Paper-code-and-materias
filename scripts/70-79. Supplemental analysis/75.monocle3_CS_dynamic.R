# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(stringr)
  library(ggplot2)
  library(reshape2)
  library(qs)
  suppressPackageStartupMessages(library(monocle3))
}

# define function
if (T) {
  sankey <- function(x, y, name_x, name_y) {
    # test data
    if (F) {
      x <- sample(1:4, 20, replace = T)
      y <- sample(1:4, 20, replace = T)
      x <- colData(cds)$seurat_clusters
      y <- partitions(cds)
      name_x <- "seurat"
      name_y <- "monocle"
    }
    
    library(networkD3)
    
    # link data prepare
    if (T) {
      df <- as.data.frame.table(table(x, y))
      colnames(df) <- c(name_x, name_y, "value")
      
      df[, name_x] <- paste0(name_x, "_", df[, name_x])
      df[, name_y] <- paste0(name_y, "_", df[, name_y])
    }
    
    # node data prepare
    if (T) {
      df.nodes <- data.frame(
        name=unique(c(df[, name_x], df[, name_y]))
      )
      
      df$IDsource <- match(df[, name_x], df.nodes$name)-1
      df$IDtarget <- match(df[, name_y], df.nodes$name)-1  
    }
    
    # sankey plot
    if (T) {
      Sankey.p <- sankeyNetwork(Links = df, Nodes = df.nodes,
                                Source = "IDsource", Target = "IDtarget", Value = "value", 
                                NodeID = "name", 
                                fontFamily="arial", 
                                nodeWidth=20, nodePadding=10, sinksRight=T, fontSize=20)
      Sankey.p
    }
    
    return(Sankey.p)
  }
}

# load seurat data
if (T) {
  seurat_wnn <- function(input, outfile) {
    # test data
    if (F) {
      input <- "data/saved_data/28.seurat_object.qs"
      outfile <- "data/saved_data/75.seurat_wnn_matrix.qs"
    }
    
    if (file.exists(outfile)) {
      mess <- paste0(outfile, ": existed.\n")
    }
    if (! file.exists(outfile)) {
      suppressPackageStartupMessages(library(Seurat))
      
      dat <- qread(input, nthreads = 6)
      wnn_embedding <- Embeddings(dat, reduction = "wnn.umap")
      
      qsave(wnn_embedding, outfile, nthreads = 6)
      
      mess <- paste0("Seurat object of wnn.umap: ", input, "\nOutput file: ", outfile, "\n******* Successfully Done *******\n")
    }
    
    cat(mess)
  }
  
  seurat_wnn(input = "data/saved_data/28.seurat_merged_object.qs", outfile = "data/saved_data/75.seurat_merged_wnn_matrix.qs")
  seurat_wnn(input = "data/saved_data/24.seurat_thp1_object.qs", outfile = "data/saved_data/75.seurat_thp1_wnn_matrix.qs")
  seurat_wnn(input = "data/saved_data/24.seurat_cd34_object.qs", outfile = "data/saved_data/75.seurat_cd34_wnn_matrix.qs")
}

for (cell in c("cd34", "thp1", "merged")) {
  # cell <- "thp1"
  # cell <- "merged"
  
  # construct monocle3 data
  if (T) {
    monocle3_cds <- function(gene_file, cluster_file, cell) {
      # test data
      if (F) {
        gene_file <- "data/saved_data/27.seurat_gene_CS_percentage_input_mat.mode1.qs"
        cluster_file <- "data/saved_data/29.common_cluster_name_for_merged_2cells_CS_clusters.qs"
      }
      
      # load cs distribution data
      if (T) {
        gene <- qread(gene_file, nthreads = 6)
      }
      
      # get the sorted matrix log transition
      if (T) {
        gene <- gene[, paste0("S", seq(min(as.numeric(gsub("S", "", colnames(gene)))), max(as.numeric(gsub("S", "", colnames(gene))))))]
        gene[1:4, 1:4]
        sapply(gene, function(x) {
          summary(x[x>0])
        })
        gene <- gene*1e4
        sapply(gene, function(x) {
          summary(x[x>0])
        })
        gene <- round(gene)
        gene <- log2(gene+1)
      }
      
      # load cluster info
      if (T) {
        cluster <- qread(cluster_file, nthreads = 6)
        cluster <- cluster[, c("orig.ident", "seurat_clusters", "cell")]
        cluster <- cluster[cluster$cell == cell, ]
        
        cluster$orig.ident <- ifelse(grepl("S", cluster$orig.ident, ignore.case = T), 
                                     cluster$orig.ident, 
                                     paste0("S", cluster$orig.ident))
        
        cluster <- cluster[match(colnames(gene), cluster$orig.ident), ]
        rownames(cluster) <- cluster$orig.ident
      }
      
      # get gene info
      if (T) {
        gene_anno <- rownames(gene)
        gene_anno <- data.frame(str_split(gene_anno, "@", simplify = T))
        if (ncol(gene_anno) == 2) {
          gene_anno <- cbind(cell, gene_anno)
        }
        colnames(gene_anno) <- c("cell", "gene_short_name", "gene_part")
        
        rownames(gene_anno) <- rownames(gene)
      }
      
      # construct the cds object
      if (T) {
        cds <- new_cell_data_set(as.matrix(gene),
                                 cell_metadata = cluster,
                                 gene_metadata = gene_anno)
      }
      
      return(cds)
    }
    
    gene_file <- ifelse(cell == "merged", 
                        "data/saved_data/27.seurat_gene_CS_percentage_input_mat.mode1.qs", 
                        paste0("data/saved_data/23.seurat_gene_CS_percentage_input_mat_", cell, ".qs"))
    cluster_file <- ifelse(cell == "merged", 
                           "data/saved_data/29.common_cluster_name_for_merged_2cells_CS_clusters.qs", 
                           "data/saved_data/25.common_cluster_name_for_2_cell_CS_clusters.qs")
    
    cds <- monocle3_cds(gene_file = gene_file, cluster_file = cluster_file, cell = cell)
  }
  
  # pre-process
  if (T) {
    cds <- preprocess_cds(cds, num_dim = 50, norm_method = c("none"))
  }
  
  # load seurat umap
  if (T) {
    wnn_embedding <- qread(paste0("data/saved_data/75.seurat_", cell, "_wnn_matrix.qs"), nthreads = 6)
  }
  
  # reduce_dimension
  if (T) {
    cds <- reduce_dimension(cds, max_components=2, reduction_method="UMAP", preprocess_method="PCA")
    
    cds.embed <- reducedDims(cds)$UMAP
    wnn_embedding <- wnn_embedding[rownames(cds.embed),]
    reducedDims(cds)$UMAP <- wnn_embedding
    
    colData(cds)
    plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters", group_cells_by="seurat_clusters", group_label_size=4, cell_size=3)
  }
  
  # clustering
  if (T) {
    cds <- cluster_cells(cds, reduction_method="UMAP", cluster_method="leiden", k=round(ncol(cds)*0.1), resolution=0.1)
    
    plot_cells(cds, reduction_method="UMAP", color_cells_by="partition", group_cells_by="partition", group_label_size=4, cell_size=3)
    plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters", group_cells_by="seurat_clusters", group_label_size=4, cell_size=3)
    
    p <- sankey(colData(cds)$seurat_clusters, partitions(cds), "seurat", "monocle")
    
    htmlwidgets::saveWidget(p, file=paste0("results/2.pic/75.Sankey_between_monocle_seurat_CS_clusters_", cell, ".html"))
    webshot::webshot(paste0("results/2.pic/75.Sankey_between_monocle_seurat_CS_clusters_", cell, ".html"), 
                     paste0("results/2.pic/75.Sankey_between_monocle_seurat_CS_clusters_", cell, ".pdf"))
  }
  
  # replace seurat cluster into clustering results
  if (T) {
    part <- ifelse(colData(cds)$seurat_clusters %in% c(1, 2, 5), "1", "2")
    names(part) <- colnames(cds)
    
    cds@clusters@listData$UMAP$clusters <- part
    cds@clusters@listData$UMAP$partitions <- part
    
    plot_cells(cds, reduction_method="UMAP", color_cells_by="partition", group_cells_by="partition", group_label_size=4, cell_size=3)
  }
  
  # structure of cell fate
  if (T) {
    cds <- learn_graph(cds, use_partition=T, learn_graph_control=list(
      "minimal_branch_len"=round(ncol(cds)*0.1),
      "nn.method"="nn2",
      "nn.k"=round(ncol(cds)*0.1),
      "nn.metric"="euclidean"
    ))
    
    p <- plot_cells(cds,
                    reduction_method="UMAP", 
                    color_cells_by = "seurat_clusters",
                    label_groups_by_cluster=F, label_leaves=T, label_branch_points=T, 
                    trajectory_graph_segment_size=1.5,
                    group_label_size=4, cell_size=3)
    file <- paste0("results/2.pic/75.cs_transition_structure_", cell, ".pdf")
    ggsave(filename = file, plot = p, width = 6, height = 5)
  }
  
  # root selection
  if (T) {
    get_earliest_principal_node <- function(cds, select.classify, my_select){
      cell_ids <- which(colData(cds)[, select.classify] == my_select)
      
      closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
      closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
      
      root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
      
      root_pr_nodes
    }
    
    pseudoLocation <- function(cds, root) {
      # root <- "1"
      cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes=get_earliest_principal_node(cds, select.classify="seurat_clusters", my_select=root))
      pseudoLoc <- pseudotime(cds)
      pseudoLoc <- data.frame(state = names(pseudoLoc), 
                              loc = pseudoLoc)
      return(pseudoLoc)
    }
    
    pseudoLoc_active <- pseudoLocation(cds, "1")
    pseudoLoc_repressive <- pseudoLocation(cds, "3")
    
    write.table(pseudoLoc_active, file = paste0("results/1.tab/75.pseudoLoc_active.", cell, ".txt"), quote = F, sep = "\t", col.names = T, row.names = F)
    write.table(pseudoLoc_repressive, file = paste0("results/1.tab/75.pseudoLoc_repressive.", cell, ".txt"), quote = F, sep = "\t", col.names = T, row.names = F)
  }
}




