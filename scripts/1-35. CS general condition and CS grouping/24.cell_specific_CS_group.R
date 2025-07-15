# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  set.seed(799)
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(stringr))
  library(ggplot2)
  library(qs)
}

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
body_bin_num <- 10
rna_levels <- 10
length_leveles <- 10
emission_file <- "data/raw_data/1.emission/chromIDEAS.emission.txt"
cells <- c(cell1, cell2)

for (cell in cells) {
  # cell <- "thp1"
  # cell <- "cd34"
  
  # read the raw data
  if (T) {
    # emission data
    emission <- read.table(emission_file, header = T, sep = "\t")
    
    # gene body state percentage
    file <- paste0("data/saved_data/23.seurat_gene_CS_percentage_input_mat_", cell, ".qs")
    gene <- qread(file, nthreads = 6)
  }
  
  # Create Seurat Object
  if (T) {
    # gene parts' percentage
    if (T) {
      # get the sorted matrix
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
      
      Assays(dat)
    }
  }
  
  # pre-qc for gene
  if (T) {
    DefaultAssay(dat) <- 'gene'
    
    dat[["gene"]]$data <- dat[["gene"]]$counts
    VariableFeatures(dat) <- rownames(dat[["gene"]])
    dat <- ScaleData(dat, features = rownames(dat))
    dat <- RunPCA(dat, features = VariableFeatures(dat), npcs = ncol(dat), verbose=F, approx=F)
    
    # Scree Plot
    if (T) {
      # get std
      if (T) {
        sd_pca <- Stdev(dat, reduction = "pca")
        
        data <- data.frame(
          PrincipalComponent = 1:length(sd_pca),
          VarianceExplained = sd_pca^2 / sum(sd_pca^2)
        )
        data$CumsumVarianceExplained <- cumsum(data$VarianceExplained)
      }
      
      # get cutoff
      if (T) {
        cutoff <- data$PrincipalComponent[which(data$CumsumVarianceExplained>0.99)[1]]
        y_percentage <- data[data$PrincipalComponent == cutoff, 3]
        y_start <- min(data$CumsumVarianceExplained)
        print(cutoff)
      }
      
      # ggplot
      if (T) {
        ggplot(data) +
          geom_line(aes(x = PrincipalComponent, y = CumsumVarianceExplained)) +
          geom_point(aes(x = PrincipalComponent, y = CumsumVarianceExplained)) +
          # geom_segment(aes(x=cutoff, y=y_start, xend=cutoff, yend=y_percentage)) +
          # geom_label(x=cutoff, y=y_percentage+0.03, label=paste0(round(y_percentage*100, 2), "%")) +
          # scale_x_continuous(breaks = sort(c(cutoff, seq(0, max(data$PrincipalComponent), 10)))) +
          labs(title = paste0("Scree Plot: ", cell), x = "Principal Component", y = "Cumulative Variance Explained") +
          cowplot::theme_cowplot()
      }
    }
    ggsave(paste0("results/2.pic/24.elbowplot_gene_", cell, ".pdf"), width = 7, height = 7)
    
    n_pc_gene <- ncol(dat@reductions$pca)
  }
  
  # pre-qc for emission
  if (T) {
    DefaultAssay(dat) <- 'emission'
    
    dat[["emission"]]$data <- dat[["emission"]]$counts
    VariableFeatures(dat) <- rownames(dat[["emission"]])
    dat <- ScaleData(dat, features = rownames(dat))
    dat <- RunPCA(dat, features = VariableFeatures(dat), npcs = nrow(dat), reduction.name = 'apca', verbose=F, approx=F)
    
    # Scree Plot
    if (T) {
      # get std
      if (T) {
        sd_pca <- Stdev(dat, reduction = "apca")
        
        data <- data.frame(
          PrincipalComponent = 1:length(sd_pca),
          VarianceExplained = sd_pca^2 / sum(sd_pca^2)
        )
        data$CumsumVarianceExplained <- cumsum(data$VarianceExplained)
      }
      
      # ggplot
      if (T) {
        ggplot(data) +
          geom_line(aes(x = PrincipalComponent, y = CumsumVarianceExplained)) +
          geom_point(aes(x = PrincipalComponent, y = CumsumVarianceExplained)) +
          labs(title = paste0("Scree Plot: ", cell), x = "Principal Component", y = "Cumulative Variance Explained") +
          cowplot::theme_cowplot()
      }
    }
    ggsave(paste0("results/2.pic/24.elbowplot_emission_", cell, ".pdf"), width = 7, height = 7)
    
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
    resolutions <- c(seq(0.1, 0.9, 0.1), 
                     seq(1, 1.8, 0.2), 
                     seq(2, 5, 1))
    dat <- FindClusters(dat, graph.name = "wsnn", algorithm = 3, resolution = resolutions, verbose = FALSE)
    
    library(clustree)
    p <- clustree(dat@meta.data, prefix = "wsnn_res.")
    
    cutoff <- ifelse(cell == "thp1", 3, 4)
    
    p + geom_hline(yintercept = cutoff)
    ggsave(paste0("results/2.pic/24.clustree_", cell, ".pdf"), width = 7, height = 10)
    
    resolution <- ifelse(cell == "thp1", 2, 1.8)
    dat <- FindClusters(dat, graph.name = "wsnn", algorithm = 3, resolution = resolution, verbose = FALSE)
  }
  
  # visualization
  if (T) {
    dat <- RunUMAP(dat, nn.name = "weighted.nn", 
                   reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
    
    # custom plot
    if (T) {
      data <- data.frame(dat@reductions$wnn.umap@cell.embeddings)
      data$state <- rownames(data)
      data$cluster <- Idents(dat)
      
      head(data)
      ggplot(data) +
        geom_point(aes(x=wnnUMAP_1, y=wnnUMAP_2, color=cluster), size=10) +
        geom_text(aes(x=wnnUMAP_1, y=wnnUMAP_2, label=state)) +
        theme(axis.title = element_text(size = rel(1.2)),
              axis.text = element_text(size = rel(1.2)),
              legend.text = element_text(size = rel(1.2)),
              legend.title = element_text(size = rel(1.2)),
              strip.text = element_text(size = rel(1.2)))
      ggsave(paste0("results/2.pic/24.states_cluster_", cell, ".pdf"), width = 7, height = 7)
    }
  }
  
  # find different for emission only
  if (T) {
    DefaultAssay(dat) <- 'emission'
    
    # find markers for every cluster compared to all remaining states
    emission.markers <- FindAllMarkers(dat,
                                       test.use = "wilcox", only.pos = F, 
                                       logfc.threshold = 0.1, return.thresh = 1,
                                       slot = "data",
                                       min.pct = 0.1,
                                       min.cells.feature = 1,
                                       min.cells.group = 1)
    write.table(emission.markers, file = paste0("results/1.tab/24.emission_differential_markers_info_", cell, ".csv"), 
                quote = F, sep = ",", col.names = T, row.names = F)
    
    VlnPlot(dat, features = sort(rownames(dat)))
    ggsave(paste0("results/2.pic/24.states_cluster_", cell, "_emission_features_vlnplot.jpeg"), width = 10, height = 10)
    
    FeaturePlot(dat, reduction = 'wnn.umap', features = sort(rownames(dat)))
    ggsave(paste0("results/2.pic/24.states_cluster_", cell, "_emission_features_pointplot.jpeg"), width = 10, height = 10)
    
    DoHeatmap(dat, features = emission.markers$gene)
    ggsave(paste0("results/2.pic/24.states_cluster_", cell, "_emission_features_heatmap.jpeg"), width = 10, height = 10)
  }
  
  # save the data
  if (T) {
    file <- paste0("results/1.tab/24.", cell, "_states_group_info.csv")
    cell_specific_group <- dat@meta.data
    write.table(cell_specific_group, file = file, quote = F, sep = ",", col.names = T, row.names = F)
    
    qsave(dat, file = paste0("data/saved_data/24.seurat_", cell, "_object.qs"), nthreads = 6)
  }
}

# 利用seurat的WNN算法对chromatin states进行聚类
# 同时利用2类数据：
#     - 染色质状态每种marker出现的可能性，数据范围[0-1]
#     - 分层抽样后基因body分成10份后，每份基因区域中不同染色质状态在该区域内的百分比，数据范围[0-1]
# 使用seurat时：
#     1. 无需normalize，因为seurat的normalize是做log转化，为了压缩数据，而自己的2类数据本身已经被压缩了，无需额外处理
#     2. ScaleData：正常进行，本质上是zscale，目的是为了去除量纲的影响。由于2类数据在分布上均存在较大差异，故此步骤正常进行
#     3. VariableFeatures：单细胞RNAseq中，使用2000个高变基因即可进行后续分析，这里已经分层挑选了2000个基因，故不做该处理
#     4. RunPCA：正常进行
#     5. Scree Plot：决定使用多少PC进行后续分析，这步非常重要，选择不同的数值会有不同的结果。
#          使用ElbowPlot函数绘制的图因为是解释方差的绝对值，这在挑选时存在主观性，会导致挑选时很难决断
#          故更改为绘制累积解释方差图进行挑选: 找出PC累计解释方差>99%的第一个PC值
#          由于组蛋白数据较少，故PC选择时使用全部PC进行后续分析
#     6. FindMultiModalNeighbors：由于染色质状态的数目比较少，故在KNN算法中使用round(ncol(dat)*0.1)作为k的值
#     7. FindClusters：分辨率设置上尝试0.1~1, 1.5, 2-5，并对所有分群结果进行可视化，挑出分群结果稳定前提下分辨率最大的值。
#          同时根据官方教程，使用algorithm = 3
# 经过上述步骤，便可以获取chromatin states的聚类结果

