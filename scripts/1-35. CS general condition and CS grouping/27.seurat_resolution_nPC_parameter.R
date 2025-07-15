# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  set.seed(799)
  library(reshape2)
  library(stringr)
  library(ggplot2)
  library(qs)
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(circlize))
}

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
body_bin_num <- 10
rna_levels <- 10
length_leveles <- 10
cells <- c(cell1, cell2)
resolutions <- c(seq(0.1, 0.9, 0.1), 
                 seq(1, 5, 0.2))

# define functions
if (T) {
  get_cell_specific_mat <- function(cell, gene_body_ID) {
    # test data
    if (F) {
      cell <- cell1
    }
    
    # CS percentage input matrix
    if (T) {
      gene <- qread(paste0("data/saved_data/21.tx_Body(", body_bin_num, "parts)_percentage_based_on_Body(", 
                           length_leveles, "parts)_RNA(", rna_levels, "parts).bs", bin_size, ".", cell, ".qs"))
    }
    
    # get cell specific genes
    if (T) {
      c <- ifelse(cell == cell1, 1, 
                  ifelse(cell == cell2, 2, NA))
      target_gene_pat <- gene_body_ID$tx_id[gene_body_ID[, paste0("cell", c, "_info")]]
    }
    
    # get target gene
    if (T) {
      rowname <- rownames(gene)
      rowname <- str_split(rowname, "@", simplify = T)[, 1]
      torf <- rowname %in% target_gene_pat
      sum(torf)
      
      gene <- gene[torf, ]
    }
    
    # format the data
    if (T) {
      rownames(gene) <- paste0(cell, "@", rownames(gene))
    }
    
    return(gene)
  }
}

# stratified random sampling gene body state percentage preparation
if (T) {
  for (mode in paste0("mode", 1)) {
    # mode <- "mode1"
    
    file <- paste0("data/saved_data/27.seurat_gene_CS_percentage_input_mat.", mode, ".qs")
    if (! file.exists(file)) {
      # get cell specific gene body metadata
      gene_body_ID <- qread(paste0("data/saved_data/26.geneBody_metadata_tx_level_stratified_random_sampling.", mode, ".qs"), nthreads = 6)
      
      gene1 <- get_cell_specific_mat(cell1, gene_body_ID)
      gene2 <- get_cell_specific_mat(cell2, gene_body_ID)
      
      gene <- rbind(gene1, gene2)
      
      qsave(gene, file, nthreads = 6)
    }
  }
}

# test all PC and resolution parameters
if (T) {
  for (mode in paste0("mode", 1)) {
    # mode <- "mode1"
    
    file <- paste0("data/saved_data/27.seurat_PC_resolution_parameters.", mode, ".qs")
    if (! file.exists(file)) {
      # WNN calculate the cluster data
      if (T) {
        # read the raw data
        if (T) {
          # emission data
          emission_file <- "data/raw_data/1.emission/chromIDEAS.emission.txt"
          emission <- read.table(emission_file, header = T, sep = "\t")
          
          # gene body state percentage
          gene <- qread(paste0("data/saved_data/27.seurat_gene_CS_percentage_input_mat.", mode, ".qs"), nthreads = 6)
        }
        
        total_n_pc <- ncol(gene)
        pc_res <- sapply(seq(10, total_n_pc), function(pc) {
          # pc <- 10
          # Create Seurat Object
          if (T) {
            # gene parts' percentage
            if (T) {
              # get the sorted matrix
              gene <- gene[, paste0("S", seq(min(as.numeric(gsub("S", "", colnames(gene)))), max(as.numeric(gsub("S", "", colnames(gene))))))]
              sapply(gene, function(x) {
                summary(x[x>0])
              })
              gene <- gene*1e4
              sapply(gene, function(x) {
                summary(x[x>0])
              })
              gene <- round(gene)
              gene <- log2(gene+1)
              
              gene <- CreateSeuratObject(counts = gene,
                                         assay = "gene",
                                         names.delim = "+",
                                         min.cells = 0, min.features = 0)
              gene@meta.data$orig.ident <- colnames(gene)
            }
            
            # emission
            if (T) {
              emission[1:4, 1:4]
              rownames(emission) <- emission[, 1]
              emission <- emission[, -c(1:2)]
              emission <- data.frame(t(emission))
              emission <- emission[, colnames(gene)]
              
              sapply(emission, function(x) {
                summary(x[x>0])
              })
              emission <- emission*1e2
              sapply(emission, function(x) {
                summary(x[x>0])
              })
              emission <- round(emission)
              emission <- log2(emission+1)
              
              emission <- CreateSeuratObject(counts = emission,
                                             assay = "emission",
                                             names.delim = "+",
                                             min.cells = 0, min.features = 0)
            }
            
            # merge the data
            if (T) {
              dat <- gene
              dat@assays$emission <- emission@assays$emission
            }
          }
          
          # pre-qc for gene
          if (T) {
            DefaultAssay(dat) <- 'gene'
            
            dat[["gene"]]$data <- dat[["gene"]]$counts
            VariableFeatures(dat) <- rownames(dat[["gene"]])
            dat <- ScaleData(dat, features = rownames(dat))
            dat <- RunPCA(dat, features = VariableFeatures(dat), npcs = ncol(dat), verbose=F, approx=F)
            n_pc_gene <- pc
          }
          
          # pre-qc for emission
          if (T) {
            DefaultAssay(dat) <- 'emission'
            
            dat[["emission"]]$data <- dat[["emission"]]$counts
            VariableFeatures(dat) <- rownames(dat[["emission"]])
            dat <- ScaleData(dat, features = rownames(dat))
            dat <- RunPCA(dat, features = VariableFeatures(dat), npcs = nrow(dat), reduction.name = 'apca', verbose=F, approx=F)
            
            n_pc_emission <- ncol(dat@reductions$apca)
          }
          
          # WNN
          if (T) {
            dat <- FindMultiModalNeighbors(
              dat, 
              reduction.list = list("pca", "apca"), 
              dims.list = list(1:n_pc_gene, 1:n_pc_emission),
              k.nn = round(ncol(dat)*0.1),
              knn.range = round(ncol(dat)*0.5),
            )
          }
          
          # find cluster plot
          if (T) {
            dat <- FindClusters(dat, graph.name = "wsnn", algorithm = 3, resolution = resolutions, verbose = FALSE)
            
            library(clustree)
            clustree(dat@meta.data, prefix = "wsnn_res.")
            ggsave(paste0("results/2.pic/27.PC_resolution/", mode, "/PC", pc, "_cluster.jpeg"), width = 6, height = 10)
          }
          
          # format the number of clusters
          if (T) {
            num <- dat@meta.data
            num <- num[, grep("wsnn_res", colnames(num))]
            colnames(num) <- gsub("wsnn_res.", "res", colnames(num))
            num <- sapply(num, function(x) {
              length(unique(x))
            })
          }
          
          return(num)
        })
        
        pc_res <- data.frame(pc_res)
        colnames(pc_res) <- paste0("PC", seq(10, total_n_pc))
      }
      
      pc_res$res <- paste0("res", resolutions)
      
      qsave(pc_res, file = file, nthreads = 6)
    }
  }
}

# plot the results
if (T) {
  for (mode in paste0("mode", 1)) {
    # mode <- "mode1"
    file <- paste0("data/saved_data/27.seurat_PC_resolution_parameters.", mode, ".qs")
    pc_res <- qread(file, nthreads = 6)
    
    pc_res$res <- gsub("res", "", pc_res$res)
    pc_res$res <- factor(pc_res$res, levels = sort(as.numeric(unique(pc_res$res))))
    
    # format the data
    if (T) {
      ggdat <- melt(pc_res, id.vars = "res")
      ggdat$variable <- factor(ggdat$variable, levels = paste0("PC", sort(as.numeric(gsub("PC", "", levels(ggdat$variable))))))
    }
    
    # get most cluster number
    if (T) {
      stat <- tapply(ggdat$value, ggdat$variable, function(x){
        # x <- ggdat[ggdat$variable == "PC28", 3]
        dat <- as.data.frame.table(table(x))
        colnames(dat) <- c("num_cluster", "num")
        max(as.vector(dat[dat$num == max(dat$num), 1]))
      })
      
      yintercept <- as.numeric(names(sort(table(stat), decreasing = T)))[1]
    }
    
    head(ggdat)
    
    # get ggplot parameters: facet nrow
    if (T) {
      len <- length(unique(ggdat$variable))
      nrow <- ifelse(len %% 6 == 0, 6, 
                     ifelse(len %% 5 == 0, 5, 
                            ifelse(len %% 4 == 0, 4, 
                                   ifelse(len %% 3 == 0, 3, 
                                          ifelse(len %% 2 == 0, 2, 3)))))
    }
    
    # plot
    if (T) {
      p <- ggplot(ggdat) +
        geom_point(aes(x=res, y=value), size=3) +
        geom_line(aes(x=res, y=value, group=variable, color=variable), linewidth=1) +
        scale_y_continuous(breaks = seq(min(ggdat$value), max(ggdat$value), 1)) +
        geom_hline(yintercept = yintercept) +
        xlab("Resolution") +
        ylab("Cluster Number") +
        facet_wrap(~variable, scales = "fixed", nrow = nrow) +
        theme(axis.title = element_text(size = rel(1.2)),
              axis.text.x = element_text(size = rel(1.1)),
              axis.text.y = element_text(size = rel(1.2)),
              legend.position = "none",
              strip.text = element_text(size = rel(1.2)))
    }
    
    # save the figure
    ggsave(filename = paste0("results/2.pic/27.seurat_PC_resolution_parameters.", mode, ".pdf"), 
           plot = p, width = 10, height = 6)
  }
}

# format clustering results when PC=37
if (T) {
  mode <- "mode1"
  resolutions <- c(seq(0.1, 0.9, 0.1), 
                   seq(1, 5, 0.1))
  
  # read the raw data
  if (T) {
    # emission data
    emission_file <- "data/raw_data/1.emission/chromIDEAS.emission.txt"
    emission <- read.table(emission_file, header = T, sep = "\t")
    
    # gene body state percentage
    gene <- qread(paste0("data/saved_data/27.seurat_gene_CS_percentage_input_mat.", mode, ".qs"), nthreads = 6)
  }
  
  # Create Seurat Object
  if (T) {
    # gene parts' percentage
    if (T) {
      # get the sorted matrix
      gene <- gene[, paste0("S", seq(min(as.numeric(gsub("S", "", colnames(gene)))), max(as.numeric(gsub("S", "", colnames(gene))))))]
      sapply(gene, function(x) {
        summary(x[x>0])
      })
      gene <- gene*1e4
      sapply(gene, function(x) {
        summary(x[x>0])
      })
      gene <- round(gene)
      gene <- log2(gene+1)
      
      gene <- CreateSeuratObject(counts = gene,
                                 assay = "gene",
                                 names.delim = "+",
                                 min.cells = 0, min.features = 0)
      gene@meta.data$orig.ident <- colnames(gene)
    }
    
    # emission
    if (T) {
      emission[1:4, 1:4]
      rownames(emission) <- emission[, 1]
      emission <- emission[, -c(1:2)]
      emission <- data.frame(t(emission))
      emission <- emission[, colnames(gene)]
      
      sapply(emission, function(x) {
        summary(x[x>0])
      })
      emission <- emission*1e2
      sapply(emission, function(x) {
        summary(x[x>0])
      })
      emission <- round(emission)
      emission <- log2(emission+1)
      
      emission <- CreateSeuratObject(counts = emission,
                                     assay = "emission",
                                     names.delim = "+",
                                     min.cells = 0, min.features = 0)
    }
    
    # merge the data
    if (T) {
      dat <- gene
      dat@assays$emission <- emission@assays$emission
    }
  }
  
  # pre-qc for gene
  if (T) {
    DefaultAssay(dat) <- 'gene'
    
    dat[["gene"]]$data <- dat[["gene"]]$counts
    VariableFeatures(dat) <- rownames(dat[["gene"]])
    dat <- ScaleData(dat, features = rownames(dat))
    dat <- RunPCA(dat, features = VariableFeatures(dat), npcs = ncol(dat), verbose=F, approx=F)
    n_pc_gene <- ncol(dat@reductions$pca)
  }
  
  # pre-qc for emission
  if (T) {
    DefaultAssay(dat) <- 'emission'
    
    dat[["emission"]]$data <- dat[["emission"]]$counts
    VariableFeatures(dat) <- rownames(dat[["emission"]])
    dat <- ScaleData(dat, features = rownames(dat))
    dat <- RunPCA(dat, features = VariableFeatures(dat), npcs = nrow(dat), reduction.name = 'apca', verbose=F, approx=F)
    
    n_pc_emission <- ncol(dat@reductions$apca)
  }
  
  # WNN
  if (T) {
    dat <- FindMultiModalNeighbors(
      dat, 
      reduction.list = list("pca", "apca"), 
      dims.list = list(1:n_pc_gene, 1:n_pc_emission),
      k.nn = round(ncol(dat)*0.1),
      knn.range = round(ncol(dat)*0.5),
    )
  }
  
  # find cluster
  if (T) {
    dat <- FindClusters(dat, graph.name = "wsnn", algorithm = 3, resolution = resolutions, verbose = FALSE)
    ggdat <- dat@meta.data
    
    ggdat <- ggdat[, c("orig.ident", paste0("wsnn_res.", resolutions))]
  }
  
  file <- paste0("data/saved_data/27.all_PC_all_resolutions_chromIDEAS_clusering.qs")
  qsave(ggdat, file, nthreads = 6)
}

# plot merged emission table and chromHMM emission with identical CS number
if (T) {
  mode <- "mode1"
  resolutions <- c(seq(0.1, 0.9, 0.1), 
                   seq(1, 1.8, 0.2), 
                   seq(2, 5, 1))
  
  # load the clustering data
  if (T) {
    file <- paste0("data/saved_data/27.all_PC_all_resolutions_chromIDEAS_clusering.qs")
    ggdat <- qread(file, nthreads = 6)
    ggdat <- ggdat[, c("orig.ident", paste0("wsnn_res.", resolutions))]
  }
  
  # clustree plot (1)
  if (T) {
    suppressPackageStartupMessages(library(clustree))
    p <- clustree(ggdat, prefix = "wsnn_res.", show_axis = T, edge_width = 1, edge_arrow = F)
    print(p)
    
    rm(p)
  }
  
  # river plot (2)
  if (F) {
    library(networkD3)
    sankey_dat <- ggdat[, -1]
    
    # format the data
    if (T) {
      for (i in 1:ncol(sankey_dat)) {
        # i <- 1
        sankey_dat[, i] <- paste0(gsub("wsnn_res.", "", colnames(sankey_dat)[i]), "C", sankey_dat[, i])
        rm(i)
      }
    }
    
    # links data
    if (T) {
      links <- apply(sankey_dat, 1, paste, collapse="@")
      links <- as.data.frame.table(table(links))
      dat <- data.frame(str_split(links$links, "@", simplify = T))
      dat$value <- links$Freq
      
      links <- lapply(1:(ncol(dat)-2), function(x) {
        aggregate(dat$value, by=list(dat[, x],dat[, (x+1)]),sum)
      })
      links <- do.call(rbind, links)
      colnames(links) <- c('source',"target","value")
    }
    
    # node data
    if (T) {
      nodes <- data.frame(name=unique(c(links$source,links$target)))
    }
    
    # id modification
    if (T) {
      # 0-indexed
      links$IDsource <- match(links$source, nodes$name)-1
      links$IDtarget <- match(links$target, nodes$name)-1
    }
    
    # sankey diagram
    if (T) {
      sankeyNetwork(Links = links, Nodes = nodes,
                    Source = "IDsource", Target = "IDtarget",
                    Value = "value", NodeID = "name", 
                    fontSize = 0, 
                    nodeWidth = 5, 
                    nodePadding = 1, 
                    sinksRight=T)
      
    }
  }
  
  # only save cluster results with highest resolution
  if (T) {
    N_cluster <- sapply(ggdat[, -1], function(x) {
      length(unique(x))
    })
    N_cluster <- data.frame(
      res = names(N_cluster), 
      n_cluster = as.numeric(N_cluster)
    )
    N_cluster <- N_cluster[match(paste0("wsnn_res.", sort(resolutions, decreasing = T)), N_cluster$res),]
    N_cluster <- N_cluster[! duplicated(N_cluster$n_cluster), ]
    ggdat <- ggdat[, c("orig.ident", N_cluster$res)]
  }
  
  # get chromIDEAS emission table
  if (T) {
    emission <- read.table("data/raw_data/1.emission/chromIDEAS.emission.txt", header = T, sep = "\t")
  }
  
  # merge chromIDEAS emission table and clustering results
  if (T) {
    emission <- emission[match(ggdat$orig.ident, emission$State), ]
    
    ggdat <- cbind(ggdat, emission[, -c(1:2)])
    
    rm(emission, N_cluster)
  }
  
  # calculate the correlation across within-group and between-group
  if (T) {
    all_res <- sort(grep("wsnn_res", colnames(ggdat), value = T))
    all_res <- all_res[-1]
    
    intra_cor <- lapply(all_res, function(res) {
      # res <- "wsnn_res.5"
      subdat <- ggdat[, !grepl("wsnn_res|orig.ident", colnames(ggdat))]
      subdat$cluster <- ggdat[, res]
      
      intra_cor <- lapply(sort(unique(subdat$cluster)), function(x) {
        # x <- 0
        mat <- subdat[subdat$cluster == x, 1:8]
        cor_dat <- cor(t(mat), method="spearman")[upper.tri(cor(t(mat)), diag=F)]
        exp <- data.frame(cluster = x, 
                          cor = cor_dat, 
                          res = res, 
                          group = "WG")
        
        set.seed(799)
        ct <- lapply(1:5, function(i_ct) {
          # i_ct <- 1
          n_ct <- sample(1:nrow(subdat), nrow(mat))
          mat <- subdat[n_ct, 1:8]
          cor_dat <- cor(t(mat), method="spearman")[upper.tri(cor(t(mat)), diag=F)]
          ct <- data.frame(cluster = x, 
                            cor = cor_dat, 
                            res = res, 
                            group = "WG_ct")
          
          return(ct)
        })
        ct <- do.call(rbind, ct)
        
        rbind(exp, ct)
      })
      intra_cor <- do.call(rbind, intra_cor)
      
      return(intra_cor)
    })
    intra_cor <- do.call(rbind, intra_cor)
    
    inter_cor <- lapply(all_res, function(res) {
      # res <- "wsnn_res.5"
      subdat <- ggdat[, !grepl("wsnn_res|orig.ident", colnames(ggdat))]
      subdat$cluster <- ggdat[, res]
      
      inter_cor <- lapply(sort(unique(subdat$cluster)), function(x) {
        # x <- 0
        
        mat_intra <- subdat[subdat$cluster == x, 1:8]
        mat_inter <- subdat[subdat$cluster != x, 1:8]
        
        cor_dat <- lapply(data.frame(t(mat_intra)), function(y) {
          # y <- t(mat_intra)[, 1]
          apply(mat_inter, 1, function(inter) {
            cor(inter, y, method="spearman")
          })
        })
        cor_dat <- do.call(c, cor_dat)
        
        exp <- data.frame(cluster = x, 
                          cor = cor_dat, 
                          res = res, 
                          group = "BG")
        
        set.seed(799)
        ct <- lapply(1:10, function(i_ct) {
          # i_ct <- 1
          n_ct <- sample(1:nrow(subdat), nrow(mat_intra))
          mat_intra <- subdat[n_ct, 1:8]
          mat_inter <- subdat[! 1:nrow(subdat) %in% n_ct, 1:8]
          
          cor_dat <- lapply(data.frame(t(mat_intra)), function(y) {
            # y <- t(mat_intra)[, 1]
            apply(mat_inter, 1, function(inter) {
              cor(inter, y, method="spearman")
            })
          })
          cor_dat <- do.call(c, cor_dat)
          
          ct <- data.frame(cluster = x, 
                           cor = cor_dat, 
                           res = res, 
                           group = "BG_ct")
          
          return(ct)
        })
        ct <- do.call(rbind, ct)
        
        rbind(exp, ct)
      })
      inter_cor <- do.call(rbind, inter_cor)
      
      return(inter_cor)
    })
    inter_cor <- do.call(rbind, inter_cor)
    
    cor_dat <- rbind(intra_cor, inter_cor)
    rm(intra_cor, inter_cor)
  }
  
  # format the data
  if (T) {
    head(cor_dat)
    
    unique(cor_dat$group)
    cor_dat$group2 <- ifelse(grepl("WG", cor_dat$group), "WG", "BG")
    cor_dat$group <- factor(cor_dat$group, levels = c("WG", "WG_ct", "BG", "BG_ct"))
    
    cor_dat$res <- factor(cor_dat$res, levels = all_res, labels = gsub("wsnn_", "", all_res))
  }
  
  # plot the correlation figure
  if (T) {
    head(cor_dat)
    
    p <- ggplot(cor_dat) +
      geom_boxplot(aes(x=group, y=cor, fill=group)) +
      facet_grid(.~res) +
      theme(axis.title = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1.2)),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "right",
            strip.text = element_text(size = rel(1.2)))
    print(p)
    
    file <- paste0("results/2.pic/27.cor_withingroup_and_betweengroup_for_all_clustering_with_highest_resolutions_results.", mode, ".pdf")
    ggsave(filename = file, plot = p, width = 15, height = 8)
    
    rm(cor_dat, p)
  }
}

for (res in gsub("wsnn_", "", all_res)) {
  bg <- cor_dat$cor[cor_dat$res == res & cor_dat$group == "BG"]
  bg_ct <- cor_dat$cor[cor_dat$res == res & cor_dat$group == "BG_ct"]
  wg_ct <- cor_dat$cor[cor_dat$res == res & cor_dat$group == "WG_ct"]
  
  if (t.test(bg, bg_ct)$p.value > 0.05 | t.test(bg, wg_ct)$p.value > 0.05) {
    print(res)
  }
}
