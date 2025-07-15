# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(reshape2)
  library(stringr)
  library(qs)
  suppressPackageStartupMessages(library(limma))
  suppressPackageStartupMessages(library(edgeR))
  suppressPackageStartupMessages(library(GenomicFeatures))
}

# setting samples
required_sample <- "cd34_|thp1_"
removed_sample <- "thp1_rep[123]"

# load the gene data
if (T) {
  dat <- qread("data/saved_data/4.gene_count_matrix_thp1_cd34.qs", nthreads = 6)
  geneid <- qread("data/saved_data/4.gene_geneid_thp1_cd34.qs", nthreads = 6)
}

# select samples
if (T) {
  head(dat)
  dat <- dat[, grepl(required_sample, colnames(dat)) & !grepl(removed_sample, colnames(dat))]
  head(dat)
  
  dat <- dat[apply(dat, 1, sum)>0, ]
  geneid <- geneid[match(rownames(dat), geneid$id), ]
}

# DEG object
if (T) {
  colnames(dat)
  x <- DGEList(counts = dat, 
               group = str_split(colnames(dat), "_", simplify = T)[, 1], 
               genes = geneid)
  
  # 29416     4
  dim(x)
  table(x$samples$group)
  total_n <- nrow(x)
  
  rm(dat, geneid, required_sample, removed_sample)
}

# quality control: filter low expression genes
if (T) {
  # define the colour
  if (T) {
    col_multi <- c("#2EC4B6", "#E71D36")
    col <- col_multi[x$samples$group]
  }
  
  filter=1
  
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
    file <- "results/2.pic/59.qc_reads_distribution_before_filter_plot.jpeg"
    jpeg(filename = file, width = 480, height = 480, units = "px")
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
    
    file <- "results/2.pic/59.qc_reads_distribution_after_filter_plot.jpeg"
    jpeg(filename = file, width = 480, height = 480, units = "px")
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
  
  rm(total_n, col_multi, col, filter, cpm, M, lcpm.cutoff, keep.exprs, den, percent, i, lcpm)
}

# normalization: TMM
if (T) {
  # raw data expression level for each sample
  if (T) {
    lcpm <- cpm(x, log=TRUE)
    
    file <- paste0("results/2.pic/59.qc_reads_distribution_before_normalization_plot.jpeg")
    jpeg(filename = file, width = 480, height = 480, units = "px")
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
    
    file <- paste0("results/2.pic/59.reads_distribution_after_normalization_plot.jpeg")
    jpeg(filename = file, width = 480, height = 480, units = "px")
    if (T) {
      boxplot(lcpm, las=2, main="",outline = F)
      title(main="B. Normalised data",ylab="Log2(cpm+1)")
    }
    dev.off()
  }
}

# quality control: PCA and Correlation
if (T) {
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
      
      ggplot(ggdat) + 
        geom_point(aes(x=Dim1, y=Dim2, color=group), size=4, alpha = 0.65) + 
        geom_text(aes(x=Dim1, y=Dim2, color=group, label=sample)) +
        theme_bw() +
        theme(legend.position="right")
      
      file <- paste0("results/2.pic/59.qc_MDS_plot.jpeg")
      ggsave(filename = file, width = 6.0, height = 4.2, units = "in")
    }
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
    file <- "results/2.pic/59.qc_correlation_pheatmap.jpeg"
    pheatmap::pheatmap(cor, scale = "none",
                       breaks = seq(min(cor)-0.1,1,length.out = 100),
                       border_color = "black", colorRampPalette(c("blue", "white", "red"))(100),
                       filename = file, width = 15, height = 10,
                       fontsize_row = 13, fontsize = 14,
                       show_colnames = F, display_numbers = T,fontsize_number = 9,
                       annotation_row = annotation_row)
  }
  
  rm(annotation_row, cor, lcpm, ggdat, mds)
}

# design matrix
if (T) {
  group <- x$samples$group
  
  design <- model.matrix( ~ group)
  
  colnames(design) <- gsub("group", "", colnames(design))
  rownames(design) <- colnames(x)
  head(design)
}

# model the data
if (T) {
  file <- paste0("results/2.pic/59.qc_mean-variance_plot.jpeg")
  jpeg(filename=file, width = 5, height = 5, units = "in", res=300)
  
  par(mfrow=c(1,2))
  v <- voom(x, design, plot=TRUE)
  vfit <- lmFit(v, design)
  efit <- eBayes(vfit)
  plotSA(efit, main="Final model: Mean-variance trend")
  dev.off()
  
  rm(v, vfit, group)
}

# DEG analysis: log2FC>1 & adj.P.Val<0.05 for THP1
if (T) {
  colnames(efit)
  # 经典提取结果的方法
  DEG_limma <- topTable(efit, coef="thp1", n=Inf)
  sum(DEG_limma$adj.P.Val<0.05)
  
  DEG_limma$state <- ifelse(DEG_limma$adj.P.Val<0.05 & abs(DEG_limma$logFC)>=1,
                            ifelse(DEG_limma$logFC>0,"Up","Down"),"NotSig")
  # Down NotSig     Up 
  # 4323   9946   5100
  table(DEG_limma$state)
}

# output the results
if (T) {
  file <- "results/1.tab/59.thp1_vs_cd34_mRNA_deg.csv"
  write.table(DEG_limma, file = file, quote = F, sep = ",", col.names = T, row.names = F)
  
  file <- "data/saved_data/59.norm_expression_mat_thp1_cd34.qs"
  qsave(x, file, nthreads = 6)
}
