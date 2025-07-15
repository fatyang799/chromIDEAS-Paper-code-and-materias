# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(UpSetR)
  library(stringr)
  library(reshape2)
  # suppressPackageStartupMessages(library(KEGG.db))
  suppressPackageStartupMessages(library(circlize))
  suppressPackageStartupMessages(library(ComplexHeatmap))
  library(qs)
}

cell1 <- "thp1"
cell2 <- "cd34"

# define function
if (T) {
  # script 58
  if (T) {
    # modify: add a color setting
    display_conversion <- function(tab, prefix="S", log2_conversion=F, color = "#E41A1C", 
                                   show_label="num+percent", n_digit=2, show_diag_lab=F, 
                                   xlab=cell1, ylab=cell2, 
                                   show_leg=F, legend_title="Gene Number") {
      # test data
      if (F) {
        prefix <- "S"
        log2_conversion <- F
        
        # "num+percent" "percent" "num"
        show_label <- "num+percent"
        n_digit <- 2
        
        show_diag_lab <- F
        xlab <- cell2
        ylab <- cell1
        legend_title <- "Gene Number"
        show_leg=F
      }
      
      # format the data
      if (T) {
        mat <- tab
        p_tab <- tab/sum(tab)*100
      }
      
      # format the data
      if (T) {
        ord <- sort(unique(as.numeric(rownames(mat))))
        mat <- mat[, as.character(ord)]
        mat <- mat[as.character(ord), ]
        mat <- as.matrix(mat)
        rownames(mat) <- colnames(mat) <- paste0(prefix, ord)
        
        p_tab <- p_tab[, as.character(ord)]
        p_tab <- p_tab[as.character(ord), ]
        p_tab <- as.matrix(p_tab)
        rownames(p_tab) <- colnames(p_tab) <- paste0(prefix, ord)
      }
      
      # log2 conversion
      if (log2_conversion) {
        mat <- log2(mat+1)
      }
      
      # color setting
      if (T) {
        col_fun <- colorRamp2(c(quantile(mat, 0.1), quantile(mat, 0.9)), 
                              c("white", color))
      }
      
      # complex heatmap
      if (T) {
        p <- Heatmap(mat, name = as.character(legend_title), col = col_fun, 
                     border_gp = gpar(col = "black"), rect_gp = gpar(col = "black"),
                     show_row_names = T, show_column_names = T, 
                     row_title = ylab, column_title = xlab, row_title_side = "left", column_title_side = "top", 
                     show_heatmap_legend = show_leg, heatmap_legend_param = list(border = T, title_position = "lefttop", legend_direction = "horizontal"),
                     row_names_max_width = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), row_names_rot = 0,
                     column_names_max_height = unit(6, "cm"), column_names_gp = gpar(fontsize = 10), column_names_rot = 90,
                     cluster_rows = F, cluster_columns = F, show_row_dend = F, show_column_dend = F, 
                     clustering_distance_columns = "euclidean", clustering_distance_rows = "euclidean", 
                     clustering_method_columns = "ward.D2", clustering_method_rows = "ward.D2", 
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       if (show_diag_lab) {
                         if(i == j) {
                           grid.text(rownames(mat)[i], x = x, y = y)
                         }
                       }
                       if (show_label == "num") {
                         if(i != j) {
                           grid.text(round(mat[i, j], n_digit), x, y, gp = gpar(fontsize = 10))
                         }
                       }
                       if (show_label == "percent") {
                         if(i != j) {
                           grid.text(paste0(round(p_tab[i, j], n_digit), "%"), x, y, gp = gpar(fontsize = 10))
                         }
                       }
                       if (show_label == "num+percent") {
                         if(i != j) {
                           grid.text(paste0(round(mat[i, j], n_digit), "\n", round(p_tab[i, j], n_digit), "%"), x, y, gp = gpar(fontsize = 10))
                         }
                       }
                     })
      }
      
      p <- ggplotify::as.ggplot(p)
      return(p)
    }
    
    plot_dcsg_loc <- function(mat, bin_n, lab, norm_location_number=T, norm_type_total_num=T, facet=T, show_leg=F) {
      # test data
      if (F) {
        bin_n <- 50
        lab <- paste0("B", seq(1, 50, length.out=8))
        norm_location_number=T
        norm_type_total_num=F
        show_leg <- F
      }
      
      # stat the number of each bin
      if (T) {
        head(mat)
        
        mat$bins <- cut(mat$location_p*100, 
                        breaks = seq(0, 100, length.out = bin_n+1), 
                        right = F, include.lowest = T, 
                        labels = paste0("B", 1:bin_n))
        
        ggdat <- mat
      }
      
      # format the data
      if (T) {
        head(ggdat)
        
        ggdat$id <- paste0(ggdat$bins, "@", ggdat$type)
        ggdat <- as.data.frame.table(table(ggdat$id))
        
        ggdat$Bin <- str_split(ggdat$Var1, "@", simplify = T)[, 1]
        ggdat$type <- str_split(ggdat$Var1, "@", simplify = T)[, 2]
        
        ggdat$Bin <- factor(ggdat$Bin, paste0("B", 1:bin_n))
      }
      
      # normalization
      if (T) {
        # location total bin number
        if (T) {
          norm_factors <- as.data.frame.table(table(mat$bins))
          ggdat$norm_factors_bin_total <- norm_factors$Freq[match(ggdat$Bin, norm_factors$Var1)]
        }
        
        # type total bin number
        if (T) {
          norm_factors <- as.data.frame.table(table(mat$type))
          ggdat$norm_factors_type_total <- norm_factors$Freq[match(ggdat$type, norm_factors$Var1)]
        }
        
        scale_factor <- ifelse(norm_location_number+norm_type_total_num == 2, 1e6, 
                               ifelse(norm_location_number+norm_type_total_num == 0, 1/1e5, 
                                      ifelse(norm_location_number, 1, 10)))
        
        ggdat$sig <- ggdat$Freq * scale_factor
        if (norm_location_number) {
          ggdat$sig <- (ggdat$Freq * scale_factor / ggdat$norm_factors_bin_total)
          
          if (norm_type_total_num) {
            ggdat$sig <- (ggdat$sig / ggdat$norm_factors_type_total)
          }
        }
        if (norm_type_total_num) {
          ggdat$sig <- (ggdat$Freq * scale_factor / ggdat$norm_factors_type_total)
          
          if (norm_location_number) {
            ggdat$sig <- (ggdat$sig / ggdat$norm_factors_bin_total)
          }
        }
      }
      
      # ggplot
      if (T) {
        head(ggdat)
        
        p <- ggplot(ggdat, aes(x=Bin, y=sig)) +
          geom_point(size=3) +
          geom_line(aes(group=type, color=type), linewidth=1) +
          scale_y_continuous(name="Occur DCSG Probability") +
          scale_x_discrete(name="Gene Relative Position", breaks=lab) +
          theme_bw() +
          theme(axis.text = element_text(size = rel(1.3)), 
                axis.title = element_text(size = rel(1.3)), 
                strip.text = element_text(size = rel(1.4)), 
                legend.position = ifelse(show_leg, "right", "none"), 
                panel.border = element_rect(color="black"), 
                strip.background = element_rect(fill=NA, color=NA))
        
        if (facet) {
          p <- p + facet_grid(type~., scales = "fixed")
        }
      }
      
      return(p)
    }
    
    mkdir_fun <- function(dir) {
      for (d in dir) {
        if (! dir.exists(d)) {
          dir.create(d, showWarnings = F, recursive = T)
        }
      }
    }
    
    one_step_go <- function(gene, prefix="58.dcsg_genes_go", filetype="pdf", 
                            ont="ALL", Count=10, P_cutoff=1, enrichment=0, showCategory=8, # filter data for go barplot
                            x="GeneRatio", orderBy="GeneRatio", decreasing=T) { # plot go barplot
      # test data
      if (F) {
        gene <- dat$id[1:100]
        prefix <- "58.dcsg_genes_go"
        filetype <- "pdf"
        
        # filter
        ont <- "ALL"
        Count <- 2
        P_cutoff <- 1
        enrichment <- 0
        showCategory <- 20
        
        # sorting and x axis: "GeneRatio", "pvalue", "p.adjust", "qvalue", "Count", 
        x <- "GeneRatio"
        orderBy <- "GeneRatio"
        decreasing <- T
      }
      
      suppressPackageStartupMessages(library(clusterProfiler))
      suppressPackageStartupMessages(library(org.Hs.eg.db))
      
      # GO enrichment analysis
      if (T) {
        file <- paste0("data/saved_data/", prefix, ".qs")
        if (file.exists(file)) {
          go <- qread(file, nthreads = 6)
        }
        if (! file.exists(file)) {
          go <- enrichGO(gene = gene,
                         OrgDb = org.Hs.eg.db,
                         keyType = 'ENSEMBL',
                         ont = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 1,
                         qvalueCutoff = 1)
          qsave(go, file = file, nthreads = 6)
        }
      }
      
      # data export (raw)
      if (! is.null(go)) {
        file <- paste0("results/1.tab/", prefix, ".csv")
        go <- DOSE::setReadable(go, OrgDb=org.Hs.eg.db, keyType='ENSEMBL')
        go_dat <- go@result
        
        # format the data
        go_dat$generatio <- sapply(go_dat$GeneRatio, function(x) {eval(parse(text = x))})
        go_dat$bgratio <- sapply(go_dat$BgRatio, function(x) {eval(parse(text = x))})
        go_dat$enrichment_Level <- log2(go_dat$generatio/go_dat$bgratio)
        
        write.table(go_dat, file = file, quote = F, sep = ",", col.names = T, row.names = F)
      } else {
        go_dat <- data.frame(
          ONTOLOGY = 0,
          ID = 0,
          Description = 0,
          GeneRatio = 0,
          BgRatio = 0,
          pvalue = 0,
          p.adjust = 0,
          qvalue = 0,
          geneID = 0,
          Count = 0,
          generatio = 0,
          bgratio = 0,
          enrichment_Level = 0
        )
      }
      
      # sorting and filtering the data
      if (T) {
        orderBy <- ifelse(orderBy == "GeneRatio", "generatio", orderBy)
        go_dat <- go_dat[order(go_dat[, orderBy], decreasing = decreasing), ]
        
        go_dat_filter <- go_dat[go_dat$Count > Count, ]
        go_dat_filter <- go_dat_filter[go_dat_filter$pvalue<P_cutoff, ]
        go_dat_filter <- go_dat_filter[go_dat_filter$enrichment_Level > enrichment, ]
      }
      
      # get target terms
      if (T) {
        get_dat <- function(go_dat, type, n) {
          # type <- "BP"
          subdat <- go_dat[go_dat$ONTOLOGY == type, ]
          subdat <- head(subdat, n)
          
          return(subdat)
        }
        
        if (ont == "ALL") {
          bp <- get_dat(go_dat_filter, "BP", showCategory)
          cc <- get_dat(go_dat_filter, "CC", showCategory)
          mf <- get_dat(go_dat_filter, "MF", showCategory)
          
          dat1 <- rbind(bp, cc, mf)
        } else {
          dat1 <- get_dat(go_dat_filter, ont, showCategory)
        }
      }
      
      # barplot (filtered)
      if (T) {
        p <- barplot_go(dat1, x="GeneRatio")
        
        file <- paste0("results/2.pic/", prefix, ".", filetype)
        ggsave(file, plot = p, width = 7, height = 10)
      }
      
      res <- list("full_list" = go_dat, "filtered_list" = go_dat_filter)
      
      return(res)
    }
    barplot_go <- function(dat1, x="GeneRatio") {
      if (nrow(dat1)>0) {
        # modify x axis label
        if (T) {
          if (x == "GeneRatio") {
            colnames(dat1)[colnames(dat1) == "generatio"] <- "x"
          } else {
            colnames(dat1)[colnames(dat1) == x] <- "x"
          }
        }
        
        # sorting the terms
        if (T) {
          dat1$Description <- factor(dat1$Description, levels = rev(dat1$Description))
        }
        
        # color setting
        if (T) {
          library(RColorBrewer)
          colors <- brewer.pal(9, "Set1") 
          
          colors <- colors[1:2]
        }
        
        # plot
        if (T) {
          p <- ggplot(dat1) +
            geom_bar(aes(x=x, y=Description, fill=pvalue), stat = "identity") +
            scale_x_continuous(name=x) +
            scale_y_discrete(name=NULL) +
            facet_grid(ONTOLOGY~., scales = "free") +
            scale_fill_gradient(high = colors[2], low = colors[1]) +
            cowplot::theme_cowplot() +
            theme(axis.text = element_text(size = rel(0.7)), 
                  panel.border = element_rect(color="black"), 
                  strip.background = element_rect(fill=NA, color=NA))
        }
      } else {
        p <- plot(x=1, y=1, main="No Data")
      }
      
      return(p)
    }
  }
}

# DEG GO
if (T) {
  # get deg
  if (T) {
    file <- "results/1.tab/59.thp1_vs_cd34_mRNA_deg.csv"
    deg <- read.table(file, header = T, sep = ",", fill = T, comment.char = "")
    
    # 9423
    deg <- deg[deg$state != "NotSig", ]
  }
  
  # format the gene id
  if (T) {
    deg$id <- str_split(deg$id, "[.]", simplify = T)[, 1]
  }
  
  # GO analysis for all significant dcsg: ont="BP" Count=ifelse(length(gene)>1000, 10, 5) showCategory=15 P_cutoff=0.001 enrichment=0
  if (T) {
    for (type in c("Up", "Down", "All")) {
      # type <- "Up"
      
      if (type == "All") {
        gene <- deg$id
      } else {
        gene <- deg$id[deg$state == type]
      }
      
      go <- one_step_go(gene, prefix=paste0("60.deg_", type, "_go"), filetype="pdf", 
                        ont="BP", Count=ifelse(length(gene)>1000, 10, 5), P_cutoff=0.001, enrichment=0, showCategory=15, 
                        x="GeneRatio", orderBy="GeneRatio", decreasing=T)
      
      rm(go, gene, type)
    }
  }
}

# DEG volcano plot
if (T) {
  # get deg
  if (T) {
    file <- "results/1.tab/59.thp1_vs_cd34_mRNA_deg.csv"
    deg <- read.table(file, header = T, sep = ",", fill = T, comment.char = "")
  }
  
  # format the dat
  if (T) {
    deg$log10p <- (-log10(deg$adj.P.Val))
    deg$id <- str_split(deg$id, "[.]", simplify = T)[, 1]
  }
  
  # label
  if (T) {
    # wnt go related genes (Up)
    if (T) {
      go_dat <- qread("data/saved_data/60.deg_Up_go.qs")
      go_dat <- go_dat@result
      wnt_dat <- go_dat[grepl("wnt", go_dat$Description, ignore.case = T), ]
      wnt_dat <- wnt_dat$geneID
      wnt_dat <- sapply(wnt_dat, str_split_1, pattern="/")
      wnt_dat <- unlist(wnt_dat)
      wnt_dat <- unique(wnt_dat)
      
      wnt_dat <- deg[deg$id %in% wnt_dat, ]
      wnt_gene1 <- wnt_dat$id[wnt_dat$logFC>4 & grepl("wnt", wnt_dat$name, ignore.case = T)]
      wnt_gene2 <- wnt_dat$id[wnt_dat$logFC>4 & grepl("tcf", wnt_dat$name, ignore.case = T)]
      wnt_gene3 <- deg$id[match(c("WNT10B", "CTNNB1"), deg$name)]
    }
    
    # immunity related genes (Down)
    if (T) {
      go_dat <- qread("data/saved_data/60.deg_Down_go.qs")
      go_dat <- go_dat@result
      immune_dat <- go_dat[grepl("immun|myeloid|erythrocyte|lymphocyte", go_dat$Description, ignore.case = T), ]
      immune_dat <- immune_dat$geneID
      immune_dat <- sapply(immune_dat, str_split_1, pattern="/")
      immune_dat <- unlist(immune_dat)
      immune_dat <- unique(immune_dat)
      
      immune_dat <- deg[deg$id %in% immune_dat, ]
      immune_dat <- deg[grepl("tal1|gata|mpo|pax5|cd34|cd4|cd8|stat|irf", deg$name, ignore.case = T) & deg$state == "Down", ]
      immune_gene <- c("MPO", "GATA1", "GATA2", "GATA3", "TAL1", "STAT5A", "CD34", "IRF4")
      immune_gene <- deg$id[deg$name %in% immune_gene]
    }
    
    label <- deg[deg$id %in% c(wnt_gene1, wnt_gene2, wnt_gene3, immune_gene), ]
  }
  
  # color setting
  if (T) {
    library(RColorBrewer)
    colors <- brewer.pal(9, "Set1") 
    
    colors <- c(colors[1:2], "grey")
    names(colors) <- c("Up", "Down", "NotSig")
  }
  
  # volcano plot
  if (T) {
    head(deg)
    
    p <- ggplot() +
      geom_point(data = deg, aes(x=logFC, y=log10p, color=state), shape=16, alpha=0.8, size=3) +
      geom_point(data=label, aes(x=logFC, y=log10p), shape=21, fill=NA, color="black", alpha=1, size=3) +
      ggrepel::geom_text_repel(data=label, aes(logFC, y=log10p, label=name), color="black", min.segment.length=0.1, max.overlaps=15) +
      geom_hline(yintercept = (-log10(0.05)), linetype=2) +
      geom_vline(xintercept = c(-1, 1), linetype=2) +
      scale_color_manual(values=colors) +
      scale_x_continuous(name="Log2FC", breaks = sort(c(seq(-15, 15, 5), -1, 1))) +
      ylab("-log10(FDR)") +
      cowplot::theme_cowplot()
  }
  
  file <- "results/2.pic/60.volcano_plot_DEG_cd34_thp1.pdf"
  ggsave(filename = file, plot = p, width = 6, height = 4)
  
  rm(p, deg, file, colors)
}

# DEG vs DCSG overlap
if (T) {
  for (data_type in c("gene")) {
    cat(paste0(data_type, ": \n"))
    
    # data_type <- "gene"
    
    # load dcsg
    if (T) {
      file <- paste0("data/saved_data/58.", data_type, "_level_DCSG_representative_CS_Conversion(chromIDEAS).qs")
      dcsg <- qread(file, nthreads = 6)
    }
    
    # get deg
    if (T) {
      file <- "results/1.tab/59.thp1_vs_cd34_mRNA_deg.csv"
      deg <- read.table(file, header = T, sep = ",", fill = T, comment.char = "")
      
      # 9423
      deg <- deg[deg$state != "NotSig", ]
    }
    
    # compare the overlap between deg and dcsg
    if (T) {
      library(ggVennDiagram)
      
      # all DCSG
      if (T) {
        venn_dat <- list(
          deg = deg$id,
          dcsg = dcsg$id
        )
        
        # color setting
        if (T) {
          library(RColorBrewer)
          colors <- brewer.pal(9, "Set1") 
          
          colors <- colors[1:2]
        }
        
        p <- ggVennDiagram(venn_dat, 
                           category.names = names(venn_dat), 
                           set_size = 8, # 字体大小
                           show_intersect = FALSE,
                           set_color = "black", # 圆圈外缘颜色
                           label = "count",
                           label_alpha = 1,
                           label_size = 6, # 数字字体大小
                           label_geom = "text",
                           label_percent_digit = 2, # 百分比小数点后数字
                           edge_lty = "solid", # 圆圈线类型
                           edge_size = 1, # 线条宽度
                           order.intersect.by = "size",
                           order.set.by = "size") + 
          scale_x_continuous(expand = expansion(mult = .2)) + 
          scale_fill_gradient(high = colors[1], low = colors[2]) +
          theme(legend.position = "none")
        
        file <- paste0("results/2.pic/60.overlap_all_deg_and_all_dcsg_chromIDEAS.pdf")
        ggsave(filename = file, plot = p, width = 8, height = 6)
      }
      
      rm(p, venn_dat, file)
    }
    
    # format the gene id
    if (T) {
      deg$id <- str_split(deg$id, "[.]", simplify = T)[, 1]
      dcsg$id <- str_split(dcsg$id, "[.]", simplify = T)[, 1]
    }
    
    # dcsg specific GO
    if (T) {
      gene <- setdiff(dcsg$id, deg$id)
      dcsg_specific_go <- one_step_go(gene, prefix="60.dcsg_specific_go", filetype="pdf", 
                                      ont="BP", Count=ifelse(length(gene)>1000, 10, 5), P_cutoff=0.001, enrichment=0, showCategory=15, 
                                      x="GeneRatio", orderBy="GeneRatio", decreasing=T)
      dcsg_specific_go_filter <- dcsg_specific_go$filtered_list
      dcsg_specific_go_full <- dcsg_specific_go$full_list
    }
    
    # deg specific GO
    if (T) {
      gene <- setdiff(deg$id, dcsg$id)
      deg_specific_go <- one_step_go(gene, prefix="60.deg_specific_go", filetype="pdf", 
                                      ont="BP", Count=ifelse(length(gene)>1000, 10, 5), P_cutoff=0.001, enrichment=0, showCategory=15, 
                                      x="GeneRatio", orderBy="GeneRatio", decreasing=T)
      deg_specific_go_filter <- deg_specific_go$filtered_list
      deg_specific_go_full <- deg_specific_go$full_list
    }
  }
}

# explain DEG using DCSG (gene number)
if (T) {
  for (data_type in c("gene")) {
    cat(paste0(data_type, ": \n"))
    
    # data_type <- "gene"
    
    # get deg and dcsg info
    if (T) {
      file <- "data/saved_data/60.dcsg_deg_merged_info.qs"
      if (file.exists(file)) {
        deg <- qread(file, nthreads = 6)
      }
      if (! file.exists(file)) {
        # load dcsg
        if (T) {
          dcsg <- qread(paste0("data/saved_data/58.", data_type, "_level_DCSG_representative_CS_Conversion(chromIDEAS).qs"), nthreads = 6)
        }
        
        # get deg
        if (T) {
          deg <- read.table("results/1.tab/59.thp1_vs_cd34_mRNA_deg.csv", header = T, sep = ",", fill = T, comment.char = "")
          
          # 9423
          deg <- deg[deg$state != "NotSig", ]
        }
        
        # merge info 
        if (T) {
          deg$CS_conversion <- dcsg$label[match(deg$id, dcsg$id)]
          deg$CS_conversion[is.na(deg$CS_conversion)] <- "nonDCSG"
        }
        
        qsave(deg, file = file, nthreads = 6)
      }
    }
    
    # heatmap
    if (T) {
      # stat function
      stat_num <- function(deg) {
        stat <- as.data.frame.table(table(deg$CS_conversion))
        stat$cell1 <- str_split(stat$Var1, "to", simplify = T)[, 2]
        stat$cell2 <- str_split(stat$Var1, "to", simplify = T)[, 1]
        stat <- dcast(stat, cell1~cell2, value.var = "Freq")
        rownames(stat) <- stat[, 1]
        stat <- stat[, -1]
        stat[is.na(stat)] <- 0
        
        return(stat)
      }
      
      # stat
      if (T) {
        head(deg)
        # 8764
        deg <- deg[deg$CS_conversion != "nonDCSG", ]
        
        up_stat <- stat_num(deg[deg$state == "Up", ])
        down_stat <- stat_num(deg[deg$state == "Down", ])
      }
      
      # color setting
      if (T) {
        library(RColorBrewer)
        colors <- brewer.pal(9, "Set1") 
        
        colors <- colors[1:2]
      }
      
      p <- display_conversion(up_stat, prefix="C", log2_conversion=F, color = colors[1], 
                         show_label="num+percent", n_digit=2, show_diag_lab=T, 
                         xlab=cell2, ylab=cell1, 
                         show_leg=F, legend_title="Diff Bin %")
      ggsave("results/2.pic/60.DEG_explain_by_DCSG_stat_Up_heatmap.pdf", plot = p, width = 6, height = 5)
      
      p <- display_conversion(down_stat, prefix="C", log2_conversion=F, color = colors[2], 
                         show_label="num+percent", n_digit=2, show_diag_lab=T, 
                         xlab=cell2, ylab=cell1, 
                         show_leg=F, legend_title="Diff Bin %")
      ggsave("results/2.pic/60.DEG_explain_by_DCSG_stat_Down_heatmap.pdf", plot = p, width = 6, height = 5)
    }
    
    # barplot
    if (T) {
      # stat
      if (T) {
        ggdat <- as.data.frame.table(table(deg$state, deg$CS_conversion))
        colnames(ggdat) <- c("State", "CS_conversion", "Num")
        
        ggdat$cell1_CSC <- str_split(ggdat$CS_conversion, "to", simplify = T)[, 2]
        ggdat$cell2_CSC <- str_split(ggdat$CS_conversion, "to", simplify = T)[, 1]
        
        ggdat$cell1_CSC <- paste0(cell1, "-", ggdat$cell1_CSC)
        ggdat$cell1_CSC <- factor(ggdat$cell1_CSC, levels = paste0(cell1, "-", 1:5))
        ggdat$cell2_CSC <- factor(ggdat$cell2_CSC, levels = 1:5)
        
        ggdat$State <- factor(ggdat$State, levels = c("Up", "Down"))
      }
      
      head(ggdat)
      
      p <- ggplot(ggdat) +
        geom_bar(aes(x=cell2_CSC, y=Num, fill=cell2_CSC), stat = "identity") +
        facet_grid(State~cell1_CSC, scales = "free", space = "free") +
        scale_y_continuous(name="Number of DEGs") +
        scale_x_discrete(name=paste0(cell2, " state cluster")) +
        cowplot::theme_cowplot() +
        theme(axis.title = element_text(size = rel(1.2)), 
              axis.text = element_text(size = rel(1.1)), 
              strip.text = element_text(size = rel(1.2)), 
              legend.position = "none", 
              panel.border = element_rect(color="black"), 
              strip.background = element_rect(fill=NA, color=NA))
      
      ggsave("results/2.pic/60.DEG_explain_by_DCSG_stat_barplot.pdf", plot = p, width = 12, height = 6)
    }
  }
}

# explain DEG using DCSG (DCSG location)
if (T) {
  for (data_type in c("gene")) {
    cat(paste0(data_type, ": \n"))
    
    # test data
    if (F) {
      data_type <- "gene"
    }
    
    method <- "chromIDEAS"
    
    # get deg and dcsg info
    if (T) {
      file <- "data/saved_data/60.dcsg_deg_location_merged_info.qs"
      if (file.exists(file)) {
        dcsg <- qread(file, nthreads = 6)
      }
      if (! file.exists(file)) {
        # load dcsg
        if (T) {
          dcsg <- qread(paste0("data/saved_data/57.", data_type, "_level_de_novo_DCSG_CS_cluster_CSdist_(", method, ").qs"), nthreads = 6)
          dcsg <- dcsg[dcsg[, cell1] != dcsg[, cell2], ]
        }
        
        # get deg
        if (T) {
          deg <- qread("data/saved_data/60.dcsg_deg_merged_info.qs")
        }
        
        # merge info
        if (T) {
          dcsg$state <- deg$state[match(dcsg$id, deg$id)]
          dcsg$CS_conversion <- deg$CS_conversion[match(dcsg$id, deg$id)]
        }
        
        # filter dat
        if (T) {
          dcsg <- dcsg[! is.na(dcsg$state), ]
          sum(is.na(dcsg$CS_conversion))
          sum(is.na(dcsg$state))
        }
        
        qsave(dcsg, file = file, nthreads = 6)
      }
    }
    
    # calculate location distribution dat
    if (T) {
      mat <- data.frame(
        location_p = dcsg$location_p, 
        type = paste0(dcsg[, cell2], "to", dcsg[, cell1]), 
        state = dcsg$state, 
        CS_conversion = dcsg$CS_conversion
      )
      rm(dcsg)
    }
    
    bin_n <- 50
    lab <- paste0("B", seq(1, 50, length.out=8))
    norm_location_number=T
    norm_type_total_num=F
    
    # stat the number of each bin
    if (T) {
      mat$bins <- cut(mat$location_p*100, 
                      breaks = seq(0, 100, length.out = bin_n+1), 
                      right = F, include.lowest = T, 
                      labels = paste0("B", 1:bin_n))
    }
    
    # format the data
    if (T) {
      head(mat)
      
      mat$id <- paste0(mat$bins, "@", mat$type, "@", mat$state, "@", mat$CS_conversion)
      mat <- as.data.frame.table(table(mat$id))
      
      mat$Bin <- str_split(mat$Var1, "@", simplify = T)[, 1]
      mat$type <- str_split(mat$Var1, "@", simplify = T)[, 2]
      mat$state <- str_split(mat$Var1, "@", simplify = T)[, 3]
      mat$CS_conversion <- str_split(mat$Var1, "@", simplify = T)[, 4]
      
      mat$Bin <- factor(mat$Bin, paste0("B", 1:bin_n))
      
      mat$cell1 <- str_split(mat$type, "to", simplify = T)[, 2]
      mat$cell2 <- str_split(mat$type, "to", simplify = T)[, 1]
      
      mat$cell1 <- factor(mat$cell1, levels = sort(unique(as.numeric(mat$cell1))))
      mat$cell2 <- factor(mat$cell2, levels = sort(unique(as.numeric(mat$cell2))))
      
      mat$group <- paste0(mat$state, "@", mat$type)
    }
    
    # normalization
    if (T) {
      # location total bin number
      if (T) {
        norm_factors <- as.data.frame.table(tapply(mat$Freq, mat$Bin, sum))
        mat$norm_factors_bin_total <- norm_factors$Freq[match(mat$Bin, norm_factors$Var1)]
      }
      
      # location total bin number
      if (T) {
        norm_factors <- as.data.frame.table(tapply(mat$Freq, mat$type, sum))
        mat$norm_factors_type_total <- norm_factors$Freq[match(mat$type, norm_factors$Var1)]
      }
      
      scale_factor <- ifelse(norm_location_number+norm_type_total_num == 2, 1e6, 
                             ifelse(norm_location_number+norm_type_total_num == 0, 1/1e5, 
                                    ifelse(norm_location_number, 1, 10)))
      
      mat$sig <- mat$Freq * scale_factor
      if (norm_location_number) {
        mat$sig <- (mat$Freq * scale_factor / mat$norm_factors_bin_total)
        
        if (norm_type_total_num) {
          mat$sig <- (mat$sig / mat$norm_factors_type_total)
        }
      }
      if (norm_type_total_num) {
        mat$sig <- (mat$Freq * scale_factor / mat$norm_factors_type_total)
        
        if (norm_location_number) {
          mat$sig <- (mat$sig / mat$norm_factors_bin_total)
        }
      }
    }
    
    # color setting
    if (T) {
      library(RColorBrewer)
      colors <- brewer.pal(9, "Set1") 
      
      colors <- colors[1:2]
      names(colors) <- c("Up", "Down")
    }
    
    # ggplot
    if (T) {
      head(mat)
      ggdat <- mat[mat$type == mat$CS_conversion, ]
      
      p <- ggplot(ggdat, aes(x=Bin, y=sig)) +
        # geom_point(size=1) +
        geom_line(aes(group=group, color=state), linewidth=1) +
        scale_y_continuous(name="Occur DCSG Probability") +
        scale_x_discrete(name="Gene Relative Position", breaks=lab) +
        scale_color_manual(values=colors) +
        facet_grid(cell1~cell2, scales = "fixed") +
        theme_bw() +
        theme(axis.text = element_text(size = rel(1.3)), 
              axis.title = element_text(size = rel(1.3)), 
              strip.text = element_text(size = rel(1.4)), 
              legend.position = "none", 
              panel.border = element_rect(color="black"), 
              strip.background = element_rect(fill=NA, color=NA))
    }
    
    file <- paste0("results/2.pic/60.", data_type, "_level_deg_CS_cluster_bin_loc_facet_by_conversion_", method, ".pdf")
    ggsave(filename = file, plot = p, width = 8, height = 6)
  }
}

# GO for explain DEG using DCSG
if (T) {
  for (data_type in c("gene")) {
    cat(paste0(data_type, ": \n"))
    
    # test data
    if (F) {
      data_type <- "gene"
    }
    
    method <- "chromIDEAS"
    
    # get dat
    if (T) {
      file <- "data/saved_data/60.dcsg_deg_merged_info.qs"
      
      dat <- qread(file, nthreads = 6)
    }
    
    # format the data
    if (T) {
      dat$id <- str_split(dat$id, "[.]", simplify = T)[, 1]
    }
    
    # mkdir
    if (T) {
      prefix <- paste0("60.DEG_explain_by_DCSG/", data_type)
      
      dir <- paste0(c(
        "results/1.tab/", 
        "results/2.pic/", 
        "data/saved_data/"), prefix)
      mkdir_fun(dir)
    }
    
    # GO analysis: Count=ifelse(length(gene)>1000, 10, 5), P_cutoff=0.001, enrichment=0
    if (T) {
      file <- paste0("data/saved_data/60.go_", method, ".qs")
      
      if (! file.exists(file)) {
        all_conversions <- sort(unique(paste0(dat$state, "@", dat$CS_conversion)))
        
        go_summary <- lapply(all_conversions, function(type) {
          # type <- all_conversions[1]
          
          # info print
          if (T) {
            n <- which(all_conversions == type)
            total <- length(all_conversions)
            mess <- paste0("\tNow process GO: (", n, "/", total, ") ", type, "\n")
            cat(mess)
          }
          
          # get info
          if (T) {
            state <- str_split(type, "@", simplify = T)[, 1]
            conversion <- str_split(type, "@", simplify = T)[, 2]
          }
          
          gene <- dat$id[dat$state == state & dat$CS_conversion == conversion]
          prefix <- paste0(prefix, "/60.GO.", method, ".", type)
          res <- one_step_go(gene, prefix=prefix, filetype="pdf", 
                             ont="ALL", Count=ifelse(length(gene)>1000, 10, 5), P_cutoff=0.001, enrichment=0, showCategory=8,
                             x="GeneRatio", orderBy="GeneRatio", decreasing=T)
          
          res <- res$filtered_list
          
          if (nrow(res)>0) {
            res$type <- type
          } else {
            res <- NULL
          }
          
          return(res)
        })
        go_summary <- do.call(rbind, go_summary)
        
        qsave(go_summary, file = file, nthreads = 6)
      }
    }
  }
}

# GO analysis summary
if (T) {
  for (data_type in c("gene")) {
    cat(paste0(data_type, ": \n"))
    
    # test data
    if (F) {
      data_type <- "gene"
    }
    
    method <- "chromIDEAS"
    
    # get simplify GO results
    if (T) {
      file <- paste0("data/saved_data/60.go_", method, ".qs")
      go_filter <- qread(file, nthreads = 6)
    }
    
    # format the data
    if (T) {
      head(go_filter)
      go_filter$log10p <- (-log10(go_filter$qvalue))
      ggdat <- dcast(go_filter, ID~type, value.var = "generatio")
      rownames(ggdat) <- ggdat[, 1]
      ggdat <- ggdat[, -1]
      
      for (i in ncol(ggdat):1) {
        ggdat <- ggdat[order(ggdat[, i], na.last = T, decreasing = T), ]
        rm(i)
      }
      ggdat <- as.matrix(ggdat)
    }
    
    # color setting
    if (T) {
      library(RColorBrewer)
      colors <- brewer.pal(9, "Set1") 
      
      colors <- colors[1:2]
      col_fun <- colorRamp2(c(quantile(ggdat[! is.na(ggdat)], 0.05), quantile(ggdat[! is.na(ggdat)], 0.95)), 
                            c("white", colors[1]))
    }
    
    # pheatmap
    if (T) {
      p <- Heatmap(ggdat, name = "GeneRatio", col = col_fun, na_col = "white", 
                   
                   border_gp = gpar(col = "black"), rect_gp = gpar(col = NA),
                   
                   show_heatmap_legend = TRUE, heatmap_legend_param = list(border = T),
                   
                   show_row_names = T, 
                   
                   column_title = NULL, row_title = NULL,
                   
                   row_names_max_width = unit(6, "cm"), row_names_gp = gpar(fontsize = 3), row_names_rot = 0,
                   column_names_max_height = unit(6, "cm"), column_names_gp = gpar(fontsize = 10), column_names_rot = 90,
                   
                   cluster_rows = F, cluster_columns = F)
      p <- ggplotify::as.ggplot(p)
      
      file <- paste0("results/2.pic/60.DEG_explain_by_DCSG/", data_type, "/60.GO_cluster_GeneRatio_heatmap.", method, ".pdf")
      ggsave(filename = file, plot = p, width = 5, height = 8)
    }
    
    # group and word cloud setting
    if (T) {
      ht_opt$message = FALSE
      
      all_types <- sort(unique(go_filter$type))
      
      for (type in all_types) {
        # type <- all_types[4]
        
        cat(paste0("Now process (", which(type == all_types), "/", length(all_types), ") ", type, "\n"))
        
        file <- paste0("results/2.pic/60.DEG_explain_by_DCSG/", data_type, "/60.GO_cluster_", method, ".", type, ".pdf")
        if (! file.exists(file)) {
          set.seed(799)
          
          go_similarity <- simplifyEnrichment::GO_similarity(go_filter$ID[go_filter$type == type], ont="BP", db="org.Hs.eg.db", measure="Rel")
          
          if (is.matrix(go_similarity)) {
            if (nrow(go_similarity) < 3) {
              cluster_res <- 1:nrow(go_similarity)
              names(cluster_res) <- rownames(go_similarity)
            } else {
              cluster_res <- simplifyEnrichment::cluster_by_pam(go_similarity, max_k = max(1, min(round(nrow(go_similarity)/20), 3)))
            }
            
            qsave(cluster_res, file = paste0("data/saved_data/60.DEG_explain_by_DCSG/", data_type, "/60.GO_cluster_", method, ".", type, ".qs"), nthreads = 6)
            
            p <- simplifyEnrichment::ht_clusters(go_similarity, cluster_res)
            
            pdf(file = file, width = 8, height = 6)
            print(p)
            dev.off()
            
            rm(type, go_similarity, cluster_res, p)
          }
        }
      }
    }
    
    # lineplot to show the top20 GO family
    if (T) {
      # get BP data only
      if (T) {
        ggdat <- go_filter[go_filter$ONTOLOGY == "BP", ]
      }
      
      # get ancestors
      if (T) {
        library(GO.db)
        library(ontologyIndex)
        go_obo <- get_ontology("go-basic.obo", propagate_relationships = c("is_a", "part_of"))
        ancestors <- lapply(ggdat$ID, get_ancestors, ontology=go_obo)
        
        for (i in 1:max(lengths(ancestors))) {
          ggdat[, paste0("parent", i)] <- sapply(ancestors, function(x) {
            ifelse(is.na(x[i]), "root", x[i])
          })
        }
      }
      
      # get top20 most GO ancestors
      if (T) {
        top20 <- names(head(sort(table(ggdat$parent5), decreasing = T), 21))
        top20 <- top20[top20 != "root"]
        
        top20 <- ggdat[ggdat$parent5 %in% top20, ]
        top20$parent5_des <- Term(top20$parent5)
        top20 <- top20[, c(3, 4, 6:14, ncol(top20))]
      }
      
      # statistics of TOP20 GO ancestors
      if (T) {
        top20 <- table(top20$parent5_des, top20$type)
        ggdat <- as.data.frame.table(top20)
        colnames(ggdat) <- c("GO", "Conversion", "Number")
      }
      
      # format the data
      if (T) {
        head(ggdat)
        ggdat$deg <- str_split(ggdat$Conversion, "@", simplify = T)[, 1]
        ggdat$Conversion <- str_split(ggdat$Conversion, "@", simplify = T)[, 2]
      }
      
      # sorting
      if (T) {
        ord <- hclust(dist(top20), method = "complete")$order
        ggdat$GO <- factor(ggdat$GO, levels = rownames(top20)[ord])
      }
      
      # ggplot
      if (T) {
        if (method == "chromIDEAS") {
          breaks <- seq(0, 30, 10)
        }
        if (method == "kmeans") {
          breaks <- seq(0, 60, 20)
        }
        
        head(ggdat)
        
        p <- ggplot(ggdat) +
          geom_bar(aes(x=Number, y=GO, fill=GO), stat="identity", color="black") +
          scale_x_continuous(name="Number of GO Terms", breaks=breaks) +
          facet_grid(deg~Conversion) +
          cowplot::theme_cowplot() +
          theme(panel.border = element_rect(color="black"), 
                legend.position = "none",
                strip.background = element_rect(fill=NA, color=NA))
      }
      
      file <- paste0("results/2.pic/60.DEG_explain_by_DCSG/", data_type, "/60.top20_GO_ancestors_stat_barplot.", method, ".pdf")
      ggsave(filename = file, plot = p, width = 15, height = 6)
    }
  }
}

# 3to4 and 3to1 specific GO
if (T) {
  for (data_type in c("gene")) {
    cat(paste0(data_type, ": \n"))
    
    for (method in c("chromIDEAS")) {
      cat(paste0("\t", method, ": \n"))
      
      # test data
      if (F) {
        data_type <- "gene"
        method <- "chromIDEAS"
      }
      
      # get simplify GO results
      if (T) {
        file <- paste0("data/saved_data/60.go_", method, ".qs")
        go_filter <- qread(file, nthreads = 6)
      }
      
      # filter the data
      if (T) {
        ggdat <- go_filter[go_filter$ONTOLOGY == "BP" & go_filter$type %in% c("Up@3to1", "Up@3to4"), ]
        ggdat <- lapply(c("Up@3to1", "Up@3to4"), function(subtype) {
          # subtype <- "Up@3to4"
          subdat <- ggdat[ggdat$type == subtype, ]
          subdat <- subdat[order(subdat$generatio, decreasing = T), ]
          if (subtype == "Up@3to1") {
            subdat <- head(subdat, 10)
          }
          if (subtype == "Up@3to4") {
            subdat1 <- head(subdat, 8)
            subdat2 <- subdat[grepl("wnt", subdat$Description, ignore.case = T) & subdat$generatio > 0.05, ]
            subdat2 <- head(subdat2, 2)
            subdat <- rbind(subdat1, subdat2)
          }
          
          return(subdat)
        })
        names(ggdat) <- c("Up@3to1", "Up@3to4")
      }
      
      # ggplot2
      if (T) {
        for (i in c("Up@3to1", "Up@3to4")) {
          # i <- "Up@3to4"
          p <- barplot_go(ggdat[[i]], x="GeneRatio")
          ggsave(filename = paste0("results/2.pic/60.DEG_explain_by_DCSG/", data_type, "/60.conversion_", i, ".", method, ".pdf"), 
                 plot = p, width = 7, height = 10)
        }
        
        rm(i, p)
      }
    }
  }
}
