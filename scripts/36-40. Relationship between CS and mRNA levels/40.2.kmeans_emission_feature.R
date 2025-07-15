# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(reshape2)
  suppressPackageStartupMessages(library(stringr))
  library(qs)
  suppressPackageStartupMessages(library(circlize))
  suppressPackageStartupMessages(library(ComplexHeatmap))
}

# copy from script32
cell1 <- "thp1"
cell2 <- "cd34"

# read the chromatin states cluster results
if (T) {
  file <- "results/1.tab/40.emission_kmeans_cluster.csv"
  cluster <- read.table(file, header = T, sep = ",", fill = T, comment.char = "")
  cluster <- cluster[, c("state", "cluster5")]
  colnames(cluster) <- c("orig.ident", "seurat_clusters")
}

# read the emission table
if (T) {
  emission <- read.table("data/raw_data/1.emission/chromIDEAS.emission.txt", 
                         header = T, sep = "\t")
  rownames(emission) <- emission[, 1]
  emission <- emission[, c("ATAC", "H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3", "H3K9me3", "H3K79me2")]
}

# get specific histone signal distribution
if (T) {
  # format the data
  if (T) {
    stat_dat <- emission
    stat_dat <- stat_dat[cluster$orig.ident, ]
    stat_dat$cluster <- cluster$seurat_clusters
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
    
    rm(stat_dat)
  }
  
  file <- "results/1.tab/40.2.emission_kmeans5_cluster_epi_feature.csv"
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
    rownames(ggdat) <- ggdat$mk
    ggdat <- ggdat[, 3:4]
    ggdat <- data.frame(t(ggdat))
    colnames(ggdat) <- gsub("[.]", "@", colnames(ggdat))
    ggdat <- ggdat[rownames(ggdat) == "log2fc", ]
  }
  
  # color setting
  if (T) {
    library(RColorBrewer)
    display.brewer.all()
    colors <- brewer.pal(9, "Set1") 
    scales::show_col(colors, labels=T)
    
    range <- range(as.numeric(matrix(ggdat)), na.rm = T)
    range <- sort(c(range, 0))
    col_fun <- colorRamp2(range, c(colors[2], "white", colors[1]))
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
    
    ggsave(filename = "results/2.pic/40.2.cluster_histone_feature_heatmap.pdf", plot = p, width = 8, height = 4)
  }
  
  rm(p, ggdat, range, col_fun, colors)
}

# replot the histone signal based on the clusters (copy from script29)
if (T) {
  # prepare the data
  if (T) {
    ggdat <- emission
    ggdat <- ggdat[cluster$orig.ident, ]
    ggdat$cluster <- cluster$seurat_clusters
    
    ggdat <- melt(ggdat, id.vars = "cluster", variable.name = "mk")
    ggdat$cluster <- factor(ggdat$cluster, levels = sort(unique(ggdat$cluster)))
    ggdat$mk <- factor(ggdat$mk, levels = c("ATAC", "H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3", "H3K9me3", "H3K79me2"))
  }
  
  # plot
  if (T) {
    head(ggdat)
    
    p <- ggplot(ggdat) +
      geom_violin(aes(x=mk, y=value, group=mk, fill=mk), scale="width") +
      geom_jitter(aes(x=mk, y=value)) +
      ylab("Histone Signal") +
      ggtitle("kmeans5") +
      facet_grid(.~cluster, scales = "free") +
      cowplot::theme_cowplot() +
      theme(axis.title = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1)),
            legend.position = "none",
            strip.background = element_rect(fill = NA), 
            strip.text = element_text(size = rel(1), face="bold"))
    
    ggsave("results/2.pic/40.2.states_cluster_kmeans5_emission_features_vlnplot.pdf", plot = p, width = 10, height = 3)
  }
}
