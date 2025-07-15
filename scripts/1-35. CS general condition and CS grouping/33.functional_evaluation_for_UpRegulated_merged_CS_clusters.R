# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(reshape2)
  library(stringr)
  library(BiocParallel)
  suppressPackageStartupMessages(library(circlize))
  library(qs)
}

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
body_bin_num <- 10
rna_levels <- 10
length_leveles <- 10
cells <- c(cell1, cell2)

# read the function data
if (T) {
  file <- "data/saved_data/32.filtered_deg_analysis_for_all_cluster_WNNmatrix_genes.qs"
  dat <- qread(file, nthreads = 6)
  
  file <- "data/saved_data/22.geneBody_metadata_tx_level_stratified_random_sampling.qs"
  info <- qread(file, nthreads = 6)
}

# format the function data
if (T) {
  # merge expression level into dat
  if (T) {
    head(dat)
    head(info)
    
    dat$gene_exp <- ifelse(dat$cell == cell1, 
                           info[match(dat$geneID, info$tx_id), cell1], 
                           info[match(dat$geneID, info$tx_id), cell2])
    
    dat$gene_exp_group <- ifelse(dat$cell == cell1, 
                                 as.vector(info[match(dat$geneID, info$tx_id), "rna_type_cell1"]), 
                                 as.vector(info[match(dat$geneID, info$tx_id), "rna_type_cell2"]))
  }
  
  # format the results
  if (T) {
    dat$cluster <- factor(dat$cluster, 
                          levels = paste0("cluster", sort(as.numeric(unique(gsub("cluster", "", dat$cluster))))))
    dat$gene_exp_group <- factor(dat$gene_exp_group, levels = c("subdat0", paste0("Q", 1:rna_levels)))
    dat$geneLoc <- factor(dat$geneLoc, levels = paste0("G", 1:length_leveles))
    
    dat$exp_loc <- paste0(dat$gene_exp_group, "_", dat$geneLoc)
    dat$exp_loc <- factor(dat$exp_loc, levels = paste0(rep(levels(dat$gene_exp_group), each=length(levels(dat$geneLoc))), 
                                                       "_", 
                                                       rep(levels(dat$geneLoc), length(levels(dat$gene_exp_group)))))
  }
  
  # transfer the exp value to log2value
  if (T) {
    non0min <- min(dat$gene_exp[dat$gene_exp>0])
    dat$log2rpkm <- log2(dat$gene_exp+non0min)
    rm(non0min)
  }
}

# 1. Volcano plot
if (T) {
  dat$logp <- -log10(dat$P_value)
  head(dat)
  
  # color setting
  if (T) {
    col_fun <- colorRamp2(c(1, length(levels(dat$gene_exp_group))), c("blue", "red"))
  }
  
  # volcano
  if (T) {
    head(dat)
    
    p <- ggplot(dat) +
      geom_point(aes(x=log2fc_mean, y=logp, color=gene_exp_group), alpha=0.5, size=2) +
      scale_colour_manual(values = col_fun(1:length(levels(dat$gene_exp_group))),
                          breaks = levels(dat$gene_exp_group)) +
      xlab("Log2(FC)") +
      ylab("-Log10(Pvalue)") +
      facet_grid(.~cluster) +
      theme_bw() +
      theme(axis.title = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1.2)),
            legend.text = element_text(size = rel(1.2)),
            legend.title = element_text(size = rel(1.2)),
            strip.text = element_text(size = rel(1.2)), 
            panel.border = element_rect(color="black"), 
            strip.background = element_rect(fill=NA, color=NA), 
            legend.position = "bottom") +
      guides(colour = guide_legend(nrow=1))
    
    file <- "results/2.pic/33.volcano_plot_for_all_gene_parts_differential_chromatin_states.jpeg"
    ggsave(filename = file, plot = p, width = 12, height = 3.5)
  }
  
  rm(p, col_fun)
}

# 2. expression level and gene location plot summary
if (T) {
  # data prepare
  if (T) {
    head(dat)
    
    ggdat <- tapply(dat$exp_loc, dat$Target_Cluster, table)
    ggdat <- lapply(ggdat, function(x) {
      # x <- as.data.frame.table(ggdat[[1]])
      x <- as.data.frame.table(x)
      x$exp_group <- str_split(x$Var1, "_", simplify = T)[, 1]
      x$gene_loc <- str_split(x$Var1, "_", simplify = T)[, 2]
      x <- x[, 2:4]
      colnames(x)[1] <- "value"
      x <- x[, c(2,3,1)]
      return(x)
    })
    ggdat <- data.frame(do.call(rbind, ggdat))
    
    ggdat$cluster <- str_split(rownames(ggdat), "[.]", simplify = T)[, 1]
    rownames(ggdat) <- 1:nrow(ggdat)
    
    head(ggdat)
    
    ggdat$exp_group <- factor(ggdat$exp_group, levels = levels(dat$gene_exp_group))
    ggdat$gene_loc <- factor(ggdat$gene_loc, levels = levels(dat$geneLoc))
    ggdat$cluster <- factor(ggdat$cluster, levels = levels(dat$cluster))
  }
  
  # get the number of diff expression level gene
  if (T) {
    gene_total_num_each_rna <- function(cell, length_leveles) {
      # test data
      if (F) {
        cell <- "cell1"
      }
      
      dat <- info[!is.na(info[, paste0(cell, "_info")]), ]
      dat <- as.data.frame.table(table(dat$rna_type_cell1))
      colnames(dat) <- c("mRNA", "Gene_num")
      
      # multiply the length number to get total gene parts
      dat$Total_num <- dat$Gene_num * length_leveles
      
      return(dat)
    }
    dat1 <- gene_total_num_each_rna("cell1", length_leveles)
    dat2 <- gene_total_num_each_rna("cell2", length_leveles)
    dat3 <- data.frame(
      mRNA = dat1$mRNA, 
      Gene_num = dat1$Gene_num + dat2$Gene_num, 
      Total_num = dat1$Total_num+dat2$Total_num
    )
    
    rm(info, dat1, dat2)
    if (F) {
      gene_num <- data.frame(Gene_type = dat1$Var1, 
                             dat1 = dat1$Freq/sum(dat1$Freq)*100, 
                             dat2 = dat2$Freq/sum(dat2$Freq)*100)
      colnames(gene_num)[2:3] <- c(cell1, cell2)
      head(gene_num)
      gene_num <- melt(gene_num, id.vars = "Gene_type", variable.name = "Cell", value.name = "Percentage")
      gene_num$Gene_type <- factor(gene_num$Gene_type, levels = c("subdat0", paste0("Q", 1:rna_levels)))
      
      ggplot(gene_num) +
        geom_bar(aes(x=Gene_type, y=Percentage, fill=Gene_type), stat="identity") +
        facet_grid(Cell~.) +
        theme_bw() +
        theme(axis.title = element_text(size = rel(1.2)),
              axis.text = element_text(size = rel(1.2)),
              legend.position = "none",
              strip.text = element_text(size = rel(1.2)))
      
      file <- "results/2.pic/47.percentage_of_various_expression_level_genes.png"
      ggsave(file, width = 10, height = 6)
      
      rm(gene_num)
    }
  }
  
  # convert the absolute number into percentage
  if (T) {
    ggdat$percentage <- ggdat$value / dat3$Total_num[match(ggdat$exp_group, dat3$mRNA)] * 100
    
    rm(dat3)
  }
  
  # load thp1 and cd34 data
  if (T) {
    dat2 <- qread(file = "data/saved_data/31.cluster_mRNA_exp_loc_line_plot_mat.qs", nthreads = 6)
    ggdat$cell <- "merged"
    dat2 <- dat2[, colnames(ggdat)]
    
    ggdat <- rbind(ggdat, dat2)
    ggdat$cell <- factor(ggdat$cell, levels = c("cd34", "merged", "thp1"))
    
    rm(dat2)
  }
  
  # ggplot2
  if (T) {
    col_fun <- colorRamp2(c(1, length(levels(ggdat$exp_group))), c("blue", "red"))
    
    head(ggdat)
    
    p <- ggplot(data=ggdat) +
      geom_line(aes(x=gene_loc, y=percentage, group=exp_group, color=exp_group), linewidth=0.9, alpha=0.8) +
      scale_colour_manual(values = col_fun(1:length(levels(ggdat$exp_group))),
                          breaks = levels(ggdat$exp_group)) +
      ylab("Differnetial gene Part Percentage in All Gene Parts") +
      scale_x_discrete(name="Gene Location", breaks=NULL) +
      facet_grid(cell~cluster, scales = "free") +
      cowplot::theme_cowplot() +
      theme(axis.title = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1)),
            legend.position = "bottom",
            panel.border = element_rect(color="black"), 
            strip.background = element_rect(color=NA, fill=NA), 
            strip.text = element_text(size = rel(1), face="bold")) +
      guides(color=guide_legend(title="Expression\nLevel", nrow=2, title.hjust=0.5))
  }
  file <- "results/2.pic/33.expression_level_and_gene_location_for_CS_clusters.pdf"
  ggsave(file, plot = p, width = 10, height = 6)
  
  rm(ggdat, col_fun, p)
}

# 3. expression based on a whole gene unit
if (T) {
  # get all loc mean expression
  if (T) {
    head(dat)
    
    ggdat <- lapply(levels(dat$cluster), function(c) {
      # c <- levels(dat$cluster)[1]
      
      # get cell and cluster specific data
      if (T) {
        subdat <- dat[dat$cluster == c, ]
      }
      
      # get location specific log2rpkm value
      if (T) {
        if (nrow(subdat)>0) {
          loc_exp <- tapply(subdat$log2rpkm, subdat$geneLoc, mean)
          loc_exp <- data.frame(loc = names(loc_exp), 
                                log2rpkm = loc_exp)
          rownames(loc_exp) <- 1:nrow(loc_exp)
        } else {
          loc_exp <- data.frame()
        }
      }
      
      # format the data
      if (T) {
        loc_exp$cluster <- c
        loc_exp$loc <- factor(loc_exp$loc, levels = paste0("G", 1:body_bin_num))
      }
      
      return(loc_exp)
    })
    ggdat <- do.call(rbind, ggdat)
  }
  
  # load cell specific data
  if (T) {
    ggdat2 <- qread("data/saved_data/31.cluster_expresssion_mat.qs", nthreads = 6)
    ggdat$cell <- "merged"
    ggdat2 <- ggdat2[, colnames(ggdat)]
    ggdat <- rbind(ggdat, ggdat2)
    
    rm(ggdat2)
  }
  
  # annotation
  if (T) {
    anno <- ggdat[ggdat$loc == "G10", ]
  }
  
  # expression value plot summary
  if (T) {
    head(ggdat)
    
    p <- ggplot(data=ggdat) +
      geom_line(aes(x=loc, y=log2rpkm, group=cluster, color=cluster), linewidth=1) +
      ggrepel::geom_text_repel(data=anno, aes(x=loc, y=log2rpkm, label=cluster, color=cluster)) +
      scale_x_discrete(name="Transcript Location", breaks=paste0("G", 1:length_leveles), 
                       labels=c("TSS(P1)", paste0("P", 2:(length_leveles-1)), "TES(P10)")) +
      facet_grid(cell~., scales = "free") +
      cowplot::theme_cowplot() +
      theme(panel.border = element_rect(color="black"), 
            strip.background = element_rect(fill=NA, color=NA), 
            legend.position = "none", 
            strip.text = element_text(size = rel(1.2)))
  }
  
  file <- "results/2.pic/33.expression_level_line_plot_for_CS_clusters.pdf"
  ggsave(filename = file, plot = p, width = 10, height = 8)
}


# 对所有gene body长度>=10的基因part进行差异分析后，接着需要根据不同cluster中的差异基因进行评估cluster的功能
# line1-67: 读入数据，并格式化数据
# line68-99: 绘制火山图展示每个cluster的差异基因数目
# line100-193: 统计每个cluster中差异基因的情况，包括基因的每个part数目以及每种表达量基因的数目
#              考虑到基因表达量为subdat0的占比高，以及不同表达量基因数目存在差异，为了消除这部分的影响
#              对于统计的绝对数值，将其除以每个表达量中所有基因的总数，即y轴为百分比
# line194-283: 对于目标cluster的差异基因part而言，如果一个基因>50%(5/10个part)部分都存在于目标cluster的差异基因part中
#              则将该基因完整纳入该目标cluster的差异基因，并以基因为单位绘制表达量的boxplot以统计不同cluster的差异基因表达量

