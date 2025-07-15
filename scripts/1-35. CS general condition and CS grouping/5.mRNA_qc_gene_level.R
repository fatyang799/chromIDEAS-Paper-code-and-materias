# quality control (gene level: count)

# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(reshape2)
  library(stringr)
  library(qs)
}

# setting samples
required_sample <- "cd34_|thp1_"

# load the gene data
if (T) {
  dat <- qread("data/saved_data/4.gene_count_matrix_thp1_cd34.qs")
  geneid <- qread("data/saved_data/4.gene_geneid_thp1_cd34.qs")
}

# select data
if (T) {
  colnames(dat)
  dat <- dat[, grepl(required_sample, colnames(dat))]
  head(dat)
}

# DEG object
if (T) {
  colnames(dat)
  x <- DGEList(counts = dat, 
               group = str_split(colnames(dat), "_", simplify = T)[, 1], 
               genes = geneid)
  rm(dat, geneid)
  
  # 61507     7
  dim(x)
  table(x$samples$group)
  
  # remove all 0 genes
  x <- x[apply(x, 1, sum)>0, ]
  # 30897     7
  dim(x)
  total_n <- nrow(x)
}

# quality control: filter low expression genes
if (T) {
  # define the colour
  if (T) {
    col_multi <- c("#2EC4B6", "#E71D36")
    col <- col_multi[x$samples$group]
  }
  
  filter=5
  
  # raw data expression distribution
  if (T) {
    # cutoff for raw data
    if (T) {
      cpm <- cpm(x)
      M <- median(x$samples$lib.size) * 1e-6
      lcpm.cutoff <- log2(filter/M)
    }
    
    # raw data
    lcpm <- cpm(x, log=TRUE)
    
    # hist plot
    file <- "results/2.pic/5.reads_distribution_before_filter_plot.jpeg"
    jpeg(filename = file, width = 8, height = 6, units = "in", res=300)
    if (T) {
      plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.30), las=2, main="", xlab="")
      title(main="A. Raw data", xlab="Log2(cpm)")
      for (i in 2:ncol(x)){
        den <- density(lcpm[,i])
        lines(den$x, den$y, col=col[i], lwd=2)
      }
      abline(v=lcpm.cutoff, lty=3)
      legend("topright", legend=unique(x$samples$group), text.col=unique(col), bty="n", cex = 1)
    }
    dev.off()
  }
  
  # filter out low expression tx
  if (T) {
    keep.exprs <- filterByExpr(x, group=x$samples$group, min.count = filter)
    sum(keep.exprs)
    
    # remove the low expression genes, so the change of lib size should be small.
    # to make the statistical model more precies, we recommend to set keep.lib.sizes=T to update the size of lib
    x <- x[keep.exprs,, keep.lib.sizes=T]
    dim(x)
    
    percent <- round((total_n-nrow(x))*100/total_n,2)
    print(paste0("filter exclude ",percent,"%"))
  }
  
  # filtered data expression distribution
  if (T) {
    lcpm <- cpm(x, log=TRUE)
    
    file <- "results/2.pic/5.reads_distribution_after_filter_plot.jpeg"
    jpeg(filename = file, width = 8, height = 6, units = "in", res=300)
    if (T) {
      plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.30), las=2, main="", xlab="")
      title(main="B. Filtered data", xlab="Log2(cpm)")
      for (i in 2:ncol(x)){
        den <- density(lcpm[,i])
        lines(den$x, den$y, col=col[i], lwd=2)
      }
      abline(v=lcpm.cutoff, lty=3)
      legend("topright", legend=unique(x$samples$group), text.col=unique(col), bty="n", cex = 1)  
    }
    dev.off()
  }
}

# normalization: TMM
if (T) {
  # raw data expression level for each sample
  if (T) {
    lcpm <- cpm(x, log=TRUE)
    
    file <- paste0("results/2.pic/5.reads_distribution_before_normalization_plot.jpeg")
    jpeg(filename = file, width = 8, height = 6, units = "in", res=300)
    if (T) {
      boxplot(lcpm, las=2, main="",outline = F)
      title(main="A. Unnormalised data",ylab="Log2(cpm+1)")  
    }
    dev.off()
  }
  # normalization
  if (T) {
    x <- calcNormFactors(x, method = "TMM")
  }
  # normalized data expression level for each sample
  if (T) {
    lcpm <- cpm(x, log=TRUE)
    
    file <- paste0("results/2.pic/5.reads_distribution_after_normalization_plot.jpeg")
    jpeg(filename = file, width = 8, height = 6, units = "in", res=300)
    if (T) {
      boxplot(lcpm, las=2, main="",outline = F)
      title(main="B. Normalised data",ylab="Log2(cpm+1)")
    }
    dev.off()
  }
}

# PCA plot
if (T) {
  lcpm <- cpm(x, log=TRUE)
  mds <- plotMDS(lcpm, labels = x$samples$group)
  # ggplot
  if (T) {
    ggdat <- cbind(data.frame(Dim1 = mds$x, Dim2 = mds$y),
                   x$samples, 
                   sample = str_split(colnames(x), "_", simplify = T)[, 2])
    head(ggdat)
    
    p <- ggplot(ggdat) + 
      geom_point(aes(x=Dim1, y=Dim2, color=group), size=4, alpha = 0.65) + 
      geom_text(aes(x=Dim1, y=Dim2, color=group, label=sample)) +
      theme_bw() +
      theme(axis.title = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1.2)),
            legend.text = element_text(size = rel(1.2)),
            legend.title = element_text(size = rel(1.2)),
            strip.text = element_text(size = rel(1.2)))
  }
  
  file <- "results/2.pic/5.MDS_plot.jpeg"
  ggsave(filename = file, plot = p, width = 8, height = 6)
}

# correlation analysis
if (T) {
  # correlation analysis
  lcpm <- cpm(x, log=TRUE)
  cor <- cor(lcpm)
  
  # annotation row
  if (T) {
    annotation_row <- data.frame(Cell = str_split(x$samples$group, "_", simplify = T)[, 1])
    rownames(annotation_row) <- rownames(cor)
  }
  
  # visualization
  file <- "results/2.pic/5.correlation_pheatmap.jpeg"
  pheatmap::pheatmap(cor, scale = "none",
                     breaks = seq(min(cor)-0.1,1,length.out = 100),
                     border_color = "black", colorRampPalette(c("blue", "white", "red"))(100),
                     filename = file, width = 8, height = 6,
                     fontsize_row = 13, fontsize = 14,
                     show_colnames = F, display_numbers = T,fontsize_number = 9,
                     annotation_row = annotation_row)
}

# 移除thp1_rep[123]
# 因为thp1_rep[123]的数据分布boxplot更宽
# 且后续如果分析THP1.tb或者THP1.d，这套个数据与thp1_rep[45]是同一实验室来源
# 故保留thp1_rep[45]

#-------------------------- PASS the quality control --------------------------#


# setting samples
required_sample <- "cd34_|thp1_"
removed_sample <- "thp1_rep[123]"

