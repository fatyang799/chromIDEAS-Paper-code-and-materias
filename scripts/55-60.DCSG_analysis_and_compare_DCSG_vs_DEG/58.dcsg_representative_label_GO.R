# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(stringr)
  library(reshape2)
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(circlize))
  suppressPackageStartupMessages(library(clusterProfiler))
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  library(qs)
}

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
up_bin_num <- 3
down_bin_num <- 3

# define function
if (T) {
  gene_cluster_conversion <- function(dat, tss_id, up_bin_num=3, down_bin_num=3) {
    # prepare multicore environment
    if (T) {
      library(parallel)
      cl <- makeCluster(6, outfile=paste0("log58.", round(rnorm(1), 2), ".txt"))
      clusterExport(cl, ls(), envir = environment())
    }
    
    # prepare hello info
    if (T) {
      start_mess <- paste0("|", paste(rep("-", 100), collapse = ""), "|\n")
      cat(start_mess)
      
      IDs <- unique(dat$id)
      breaks <- round(seq(1, length(IDs), length.out=100))
      breaks <- IDs[breaks]
    }
    
    # format the data
    if (T) {
      dat$cs_trans <- paste0(dat[, cell2], "to", dat[, cell1])
    }
    
    # all types
    if (T) {
      all_conversion <- sort(unique(dat$cs_trans))
    }
    
    # summary the cs conversion situation
    if (T) {
      cs_trans <- parLapply(cl, IDs, function(gene_id) {
        # gene_id <- IDs[1]
        
        # get subdat
        if (T) {
          subdat <- dat[dat$id == gene_id, ]
        }
        
        # print process
        if (gene_id %in% breaks) {
          if (which(gene_id == breaks) == 1) {
            cat("|*")
          }
          if (which(gene_id == breaks) == 100) {
            cat("*|\n")
          }
          if (! which(gene_id == breaks) %in% c(1, 100)) {
            cat("*")
          }
        }
        
        # stat distance
        if (T) {
          convers_dist <- tapply(subdat$distance, subdat$cs_trans, sum)
          convers_dist <- ifelse(all_conversion %in% names(convers_dist), convers_dist[all_conversion], 0)
          names(convers_dist) <- all_conversion
        }
        
        # stat loction
        if (T) {
          tss <- unique(subdat$start)
          tes <- unique(subdat$end)
          strand <- ifelse(tss<tes, "+", "-")
          tss_up <- ifelse(strand == "+", tss - up_bin_num,
                           ifelse(strand == "-", tss + up_bin_num, NA))
          tes_down <- ifelse(strand == "+", tes + down_bin_num,
                             ifelse(strand == "-", tes - down_bin_num, NA))
          
          gene_bin_id <- (tss_up:tes_down)
          dcsg_bin_id <- gene_bin_id[subdat$location_id]
          
          loc_type <- ifelse(dcsg_bin_id %in% tss_id[[gene_id]], "TSS", "Body")
          convers_tss <- tapply(loc_type, subdat$cs_trans, function(x) {
            sum(x == "TSS")
          })
          convers_tss <- ifelse(all_conversion %in% names(convers_tss), convers_tss[all_conversion], 0)
          names(convers_tss) <- all_conversion
        }
        
        # summary the res
        if (T) {
          max_dist <- max(convers_dist)
          max_dist_type <- names(convers_dist)[convers_dist == max_dist]
          max_dist_N <- length(max_dist_type)
          max_dist_conversion = paste(max_dist_type, collapse = "@")
          
          max_tss <- max(convers_tss)
          max_tss_type <- names(convers_tss)[convers_tss == max_tss]
          max_tss_N <- length(max_tss_type)
          max_tss_conversion = paste(max_tss_type, collapse = "@")
          
          stat <- data.frame(
            id = gene_id, 
            max_dist = max_dist, 
            max_dist_N = max_dist_N, 
            max_dist_conversion = max_dist_conversion, 
            max_tss = max_tss, 
            max_tss_N = max_tss_N, 
            max_tss_conversion = max_tss_conversion)
        }
        
        return(stat)
      })
    }
    
    # end multicore environment
    if (T) {
      stopCluster(cl)
    }
    
    # merge and format the data
    if (T) {
      cs_trans <- do.call(rbind, cs_trans)
    }
    
    return(cs_trans)
  }
  
  display_conversion <- function(tab, prefix="S", log2_conversion=F, 
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
      library(RColorBrewer)
      colors <- brewer.pal(9, "Set1") 
      
      colors <- colors[1]
      col_fun <- colorRamp2(c(quantile(mat, 0.1), quantile(mat, 0.9)), 
                            c("white", colors[1]))
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
        p <- p + facet_grid(type~., scales = "free")
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
  one_step_kegg <- function(gene, prefix="58.dcsg_genes_kegg", filetype="pdf", 
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
    suppressPackageStartupMessages(library(KEGG.db))
    
    # id conversion
    if (T) {
      id <- bitr(geneID=gene, fromType = "ENSEMBL", toType = "ENTREZID", org.Hs.eg.db)
    }
    
    # kegg enrichment analysis
    if (T) {
      file <- paste0("data/saved_data/", prefix, ".qs")
      if (file.exists(file)) {
        kegg <- qread(file, nthreads = 6)
      }
      if (! file.exists(file)) {
        kegg <- enrichKEGG(gene = id$ENTREZID,
                           organism = "hsa",
                           keyType = 'kegg',
                           pAdjustMethod = "BH",
                           pvalueCutoff = 1,
                           qvalueCutoff = 1, 
                           minGSSize = 10,
                           maxGSSize = 5000, 
                           use_internal_data = T)
        qsave(kegg, file = file, nthreads = 6)
      }
    }
    
    # data export (raw)
    if (! is.null(kegg)) {
      file <- paste0("results/1.tab/", prefix, ".csv")
      kegg <- DOSE::setReadable(kegg, OrgDb=org.Hs.eg.db, keyType='ENTREZID')
      kegg_dat <- kegg@result
      
      # format the data
      kegg_dat$generatio <- sapply(kegg_dat$GeneRatio, function(x) {eval(parse(text = x))})
      kegg_dat$bgratio <- sapply(kegg_dat$BgRatio, function(x) {eval(parse(text = x))})
      kegg_dat$enrichment_Level <- log2(kegg_dat$generatio/kegg_dat$bgratio)
      
      write.table(kegg_dat, file = file, quote = F, sep = ",", col.names = T, row.names = F)
    } else {
      kegg_dat <- data.frame(
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
      kegg_dat <- kegg_dat[order(kegg_dat[, orderBy], decreasing = decreasing), ]
      
      kegg_dat_filter <- kegg_dat[kegg_dat$Count > Count, ]
      kegg_dat_filter <- kegg_dat_filter[kegg_dat_filter$pvalue<P_cutoff, ]
      kegg_dat_filter <- kegg_dat_filter[kegg_dat_filter$enrichment_Level > enrichment, ]
    }
    
    # get target terms
    if (T) {
      dat1 <- head(kegg_dat_filter, showCategory)
    }
    
    # barplot (filtered)
    if (T) {
      p <- barplot_kegg(dat1, x="GeneRatio")
      
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
  barplot_kegg <- function(dat1, x="GeneRatio") {
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

# representative state conversion
if (T) {
  for (data_type in c("gene")) {
    cat(paste0(data_type, ": \n"))
    
    for (method in c("chromIDEAS", "kmeans")) {
      cat(paste0("\t", method, ": \n"))
      
      # test data
      if (F) {
        data_type <- "gene"
        method <- "chromIDEAS"
      }
      
      # get DCSG representative state conversion
      if (T) {
        file <- paste0("data/saved_data/58.", data_type, "_level_DCSG_representative_CS_Conversion(", method, ").qs")
        
        if (! file.exists(file)) {
          
          # get TSS bin IDs
          if (T) {
            tss_id <- qread("data/saved_data/11.tx_level_TSS_region_bins.qs", nthreads = 6)
            gtf <- qread("data/saved_data/4.tx_geneid_thp1_cd34.qs", nthreads = 6)
            
            tss_id$gene_id <- gtf$gene_id[match(tss_id$tx_id, gtf$tx_id)]
            
            tss_id <- split(tss_id$Bin_ID, tss_id$gene_id)
            rm(gtf)
          }
          
          # get Gene bin IDs
          if (T) {
            genebody_stat <- qread(paste0("data/saved_data/11.", data_type, "_level_Body_region_length_ge3_bins.qs"), nthreads = 6)
          }
          
          # get distance based on cs cluster
          if (T) {
            dat <- qread(paste0("data/saved_data/57.", data_type, "_level_de_novo_DCSG_CS_cluster_CSdist_(", method, ").qs"), nthreads = 6)
            dat <- dat[dat[, cell1] != dat[, cell2], ]
          }
          
          # merge the data
          if (T) {
            dat$start <- genebody_stat$tss_BinID[match(dat$id, genebody_stat$tx_id)]
            dat$end <- genebody_stat$tes_BinID[match(dat$id, genebody_stat$tx_id)]
            
            rm(genebody_stat)
          }
          
          # calculate the representative conversion
          if (T) {
            cs_conversion <- gene_cluster_conversion(dat, tss_id)
          }
          
          # assign the CS cluster conversion label
          if (T) {
            cs_conversion$label <- ifelse(cs_conversion$max_dist_N == 1, cs_conversion$max_dist_conversion, 
                                          ifelse(cs_conversion$max_tss_N == 1, cs_conversion$max_tss_conversion, "confused"))
            cs_conversion <- cs_conversion[cs_conversion$label != "confused", ]
          }
          
          qsave(cs_conversion, file = file, nthreads = 6)
        }
      }
    }
  }
}

# CS conversion visualization
if (T) {
  for (data_type in c("gene")) {
    cat(paste0(data_type, ": \n"))
    
    for (method in c("chromIDEAS", "kmeans")) {
      cat(paste0("\t", method, ": \n"))
      
      # test data
      if (F) {
        data_type <- "gene"
        method <- "chromIDEAS"
      }
      
      # mkdir
      dir <- paste0("results/2.pic/58.DCSG_conversion/", data_type)
      mkdir_fun(dir)
      
      # ------------- bin based ------------- #
      
      # load differential CS state genes
      if (T) {
        file <- paste0("data/saved_data/57.", data_type, "_level_de_novo_DCSG_CS_(", method, "_dist).qs")
        dat <- qread(file, nthreads = 6)
      }
      
      # CS dcsg location distribution
      if (T) {
        mat <- data.frame(
          location_p = dat$location_p, 
          type = ifelse(dat[, cell1] == dat[, cell2], "same", "dcsg")
        )
        
        p <- plot_dcsg_loc(mat, bin_n=50, 
                           lab=paste0("B", seq(1, 50, length.out=8)), 
                           norm_location_number=T, norm_type_total_num=F, 
                           facet=T, show_leg=T)
        
        file <- paste0(dir, "/58.", data_type, "_level_dcsg_CS_bin_loc_", method, ".pdf")
        ggsave(filename = file, plot = p, width = 10, height = 8)
      }
      
      # CS conversion number heatmap
      if (T) {
        dat <- dat[dat[, cell1] != dat[, cell2], ]
        
        # heatmap
        if (T) {
          tab <- as.data.frame.matrix(table(dat$thp1, dat$cd34))
          
          p <- display_conversion(tab, prefix="S", log2_conversion=F, 
                                  show_label="none", n_digit=2, show_diag_lab=F, 
                                  xlab=cell2, ylab=cell1, 
                                  show_leg=T, legend_title="Diff Bin %")
        }
        
        file <- paste0(dir, "/58.", data_type, "_level_dcsg_CS_bin_number_conversion_", method, ".pdf")
        ggsave(filename = file, plot = p, width = 10, height = 8)
      }
      
      
      
      # load differential CS state genes
      if (T) {
        file <- paste0("data/saved_data/57.", data_type, "_level_de_novo_DCSG_CS_cluster_CSdist_(", method, ").qs")
        dat <- qread(file, nthreads = 6)
      }
      
      # CS cluster dcsg location distribution
      if (T) {
        mat <- data.frame(
          location_p = dat$location_p, 
          type = ifelse(dat[, cell1] == dat[, cell2], "same", "dcsg")
        )
        
        p <- plot_dcsg_loc(mat, bin_n=50, 
                           lab=paste0("B", seq(1, 50, length.out=8)), 
                           norm_location_number=T, norm_type_total_num=F, 
                           facet=T, show_leg=T)
        
        file <- paste0(dir, "/58.", data_type, "_level_dcsg_CS_cluster_bin_loc_", method, ".pdf")
        ggsave(filename = file, plot = p, width = 10, height = 8)
      }
      
      # CS cluster conversion number heatmap
      if (T) {
        dat <- dat[dat[, cell1] != dat[, cell2], ]
        
        # heatmap
        if (T) {
          tab <- table(dat$thp1, dat$cd34)
          
          p <- display_conversion(tab, prefix="C", log2_conversion=F, 
                                  show_label="num+percent", n_digit=2, show_diag_lab=T, 
                                  xlab=cell2, ylab=cell1, 
                                  show_leg=F, legend_title="Diff Bin %")
        }
        
        file <- paste0(dir, "/58.", data_type, "_level_dcsg_CS_cluster_bin_number_conversion_", method, ".pdf")
        ggsave(filename = file, plot = p, width = 10, height = 8)
      }
      
      # ------------- gene based ------------- #
      
      # get CS conversion
      if (T) {
        file <- paste0("data/saved_data/58.", data_type, "_level_DCSG_representative_CS_Conversion(", method, ").qs")
        
        dat <- qread(file, nthreads = 6)
      }
      
      # stat CS conversion
      if (T) {
        tab <- table(dat$label)
      }
      
      # format the data
      if (T) {
        tab <- as.data.frame.table(tab)
        tab[, cell1] <- str_split(tab$Var1, "to", simplify = T)[, 2]
        tab[, cell2] <- str_split(tab$Var1, "to", simplify = T)[, 1] 
        tab <- dcast(tab, thp1~cd34, value.var = "Freq")
        tab <- tab[, -1]
        tab <- as.matrix(tab, rownames.force=T)
        
        tab[is.na(tab)] <- 0
      }
      
      # CS conversion number heatmap
      if (T) {
        # heatmap
        if (T) {
          p <- display_conversion(tab, prefix="C", log2_conversion=F, 
                                  show_label="num+percent", n_digit=2, show_diag_lab=T, 
                                  xlab=cell2, ylab=cell1, 
                                  show_leg=F, legend_title="Diff Bin %")
        }
        
        file <- paste0(dir, "/58.", data_type, "_level_dcsg_gene_number_conversion_", method, ".pdf")
        ggsave(filename = file, plot = p, width = 10, height = 8)
      }
    }
  }
}

# cluster conversion location distribution (based on "CS conversion visualization")
if (T) {
  for (data_type in c("gene")) {
    cat(paste0(data_type, ": \n"))
    
    for (method in c("chromIDEAS", "kmeans")) {
      cat(paste0("\t", method, ": \n"))
      
      # test data
      if (F) {
        data_type <- "gene"
        method <- "chromIDEAS"
      }
      
      # mkdir
      dir <- paste0("results/2.pic/58.DCSG_conversion/", data_type)
      mkdir_fun(dir)
      
      # load differential CS state genes
      if (T) {
        file <- paste0("data/saved_data/57.", data_type, "_level_de_novo_DCSG_CS_cluster_CSdist_(", method, ").qs")
        dat <- qread(file, nthreads = 6)
        dat <- dat[dat[, cell1] != dat[, cell2], ]
      }
      
      # CS cluster dcsg location distribution
      if (T) {
        mat <- data.frame(
          location_p = dat$location_p, 
          type = paste0(dat[, cell2], "to", dat[, cell1])
        )
        rm(dat)
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
        
        mat$id <- paste0(mat$bins, "@", mat$type)
        mat <- as.data.frame.table(table(mat$id))
        
        mat$Bin <- str_split(mat$Var1, "@", simplify = T)[, 1]
        mat$type <- str_split(mat$Var1, "@", simplify = T)[, 2]
        
        mat$Bin <- factor(mat$Bin, paste0("B", 1:bin_n))
        
        mat$cell1 <- str_split(mat$type, "to", simplify = T)[, 2]
        mat$cell2 <- str_split(mat$type, "to", simplify = T)[, 1]
        
        mat$cell1 <- factor(mat$cell1, levels = sort(unique(as.numeric(mat$cell1))))
        mat$cell2 <- factor(mat$cell2, levels = sort(unique(as.numeric(mat$cell2))))
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
      
      # ggplot
      if (T) {
        head(mat)
        
        p <- ggplot(mat, aes(x=Bin, y=sig)) +
          geom_point(size=1) +
          geom_line(aes(group=type, color=type), linewidth=1) +
          scale_y_continuous(name="Occur DCSG Probability") +
          scale_x_discrete(name="Gene Relative Position", breaks=lab) +
          facet_grid(cell1~cell2, scales = "fixed") +
          theme_bw() +
          theme(axis.text = element_text(size = rel(1.3)), 
                axis.title = element_text(size = rel(1.3)), 
                strip.text = element_text(size = rel(1.4)), 
                legend.position = "none", 
                panel.border = element_rect(color="black"), 
                strip.background = element_rect(fill=NA, color=NA))
      }
      
      file <- paste0(dir, "/58.", data_type, "_level_dcsg_CS_cluster_bin_loc_facet_by_conversion_", method, ".pdf")
      ggsave(filename = file, plot = p, width = 8, height = 6)
    }
  }
}

# number statistics
if (T) {
  for (data_type in c("gene")) {
    cat(paste0(data_type, ": \n"))
    
    for (method in c("chromIDEAS", "kmeans")) {
      cat(paste0("\t", method, ": \n"))
      
      # test data
      if (F) {
        data_type <- "gene"
        method <- "kmeans"
      }
      
      # mkdir
      dir <- paste0("results/2.pic/58.DCSG_conversion/", data_type)
      mkdir_fun(dir)
      
      # load differential CS state genes
      if (T) {
        file <- paste0("data/saved_data/57.", data_type, "_level_de_novo_DCSG_CS_(", method, "_dist).qs")
        dat_cs <- qread(file, nthreads = 6)
        dat_cs$type <- ifelse(dat_cs[, cell1] == dat_cs[, cell2], 0, 1)
      }
      
      # load differential CS cluster genes
      if (T) {
        file <- paste0("data/saved_data/57.", data_type, "_level_de_novo_DCSG_CS_cluster_CSdist_(", method, ").qs")
        dat_cluster <- qread(file, nthreads = 6)
        dat_cluster$type <- ifelse(dat_cluster[, cell1] == dat_cluster[, cell2], 0, 1)
      }
      
      # DCSG and non_DCSG number
      if (T) {
        # statistics
        if (T) {
          stat_dcsg <- function(dat) {
            ggdat <- tapply(dat$type, dat$id, sum)
            ggdat <- data.frame(id = names(ggdat), 
                                type = ifelse(ggdat>0, "DCSG", "nonDCSG"))
            
            return(ggdat)
          }
          
          ggdat_cs <- stat_dcsg(dat_cs)
          ggdat_cluster <- stat_dcsg(dat_cluster)
          
          # make the order identical
          ggdat_cs <- ggdat_cs[match(ggdat_cluster$id, ggdat_cs$id), ]
          
          ggdat <- as.data.frame.table(table(paste0(ggdat_cluster$type, "_Cluster"), paste0(ggdat_cs$type, "_CS")))
        }
        
        # output gene ID
        if (T) {
          # 7693 / 12156
          dcsg_nondcscg <- ggdat_cs$id[ggdat_cs$type == "DCSG" & ggdat_cluster$type == "nonDCSG"]
          write.table(dcsg_nondcscg, file = paste0("results/1.tab/58.DCSG_conversion/", data_type, "/58.DCSG_nonDCSCG_geneID.", method, ".txt"), 
                      quote = F, sep = "\t", col.names = F, row.names = F)
          
          # 37151 / 32688
          dcscg <- ggdat_cs$id[ggdat_cluster$type == "DCSG"]
          write.table(dcscg, file = paste0("results/1.tab/58.DCSG_conversion/", data_type, "/58.DCSCG_geneID.", method, ".txt"), 
                      quote = F, sep = "\t", col.names = F, row.names = F)
        }
        
        head(ggdat)
        
        p <- ggplot(ggdat) +
          geom_bar(aes(x=Var1, y=Freq, fill=Var2), stat = "identity") +
          labs(fill = "Type") +
          xlab(NULL) +
          ylab(NULL) +
          cowplot::theme_cowplot()
        
        file <- paste0(dir, "/58.", data_type, "_level_number_of_DCSG_VS_nonDCSG_", method, ".pdf")
        ggsave(filename = file, plot = p, width = 8, height = 6)
      }
      
      # DCSG percentage per gene
      if (T) {
        # statistics
        if (T) {
          stat_dcsg <- function(dat) {
            gene_dcsg_num <- tapply(dat$type, dat$id, sum)
            gene_dcsg_num <- data.frame(id = names(gene_dcsg_num), 
                                        dcsg_length = gene_dcsg_num)
            
            gene_length_num <- dat[dat$location_id == 1, ]
            gene_length_num$gene_length <- gene_length_num$location_id / gene_length_num$location_p
            
            # make order identical
            gene_dcsg_num <- gene_dcsg_num[gene_length_num$id, ]
            
            if (identical(gene_dcsg_num$id, gene_length_num$id)) {
              ggdat <- data.frame(id = gene_dcsg_num$id, 
                                  dcsg_length = gene_dcsg_num$dcsg_length, 
                                  gene_length = gene_length_num$gene_length)
              ggdat$dcsg_length_p <- ggdat$dcsg_length / ggdat$gene_length * 100
              ggdat$type <- ifelse(ggdat$dcsg_length > 0, "DCSG", "nonDCSG")
            }
            
            return(ggdat)
          }
          
          ggdat_cs <- stat_dcsg(dat_cs)
          ggdat_cluster <- stat_dcsg(dat_cluster)
          
          # make the order identical
          ggdat_cs <- ggdat_cs[match(ggdat_cluster$id, ggdat_cs$id), ]
          
          ggdat <- data.frame(
            id = ggdat_cs$id, 
            type_Cluster = ggdat_cluster$type, 
            type_CS = ggdat_cs$type, 
            dcsg_cs_p = ggdat_cs$dcsg_length_p, 
            dcsg_cluster_p = ggdat_cluster$dcsg_length_p
          )
          ggdat$group <- paste0(ggdat$type_Cluster, "@", ggdat$type_CS)
          ggdat <- ggdat[ggdat$dcsg_cs_p>0, ]
        }
        
        head(ggdat)
        
        p <- ggplot(ggdat, aes(x=group, y=dcsg_cs_p, group=group, fill=group)) +
          geom_violin(scale = "width") +
          geom_boxplot(width=0.1, fill=NA) +
          ggpubr::stat_compare_means(method="wilcox.test", paired = F, label="p.signif", comparisons=list(unique(ggdat$group))) +
          labs(fill = "type_CS") +
          xlab(NULL) +
          ylab("DCSG Percentage Per Gene") +
          cowplot::theme_cowplot() +
          theme(axis.title = element_text(size = rel(1.2)),
                axis.text = element_text(size = rel(1.2)),
                legend.position = "none")
        
        file <- paste0(dir, "/58.", data_type, "_level_DCS_percentage_per_Gene_", method, ".pdf")
        ggsave(filename = file, plot = p, width = 6, height = 4)
      }
    }
  }
}

# GO analysis
if (T) {
  for (data_type in c("gene")) {
    cat(paste0(data_type, ": \n"))
    
    for (method in c("chromIDEAS", "kmeans")) {
      cat(paste0("\t", method, ": \n"))
      
      # test data
      if (F) {
        data_type <- "gene"
        method <- "chromIDEAS"
      }
      
      # get CS conversion
      if (T) {
        file <- paste0("data/saved_data/58.", data_type, "_level_DCSG_representative_CS_Conversion(", method, ").qs")
        
        dat <- qread(file, nthreads = 6)
      }
      
      # format the data
      if (T) {
        dat$id <- str_split(dat$id, "[.]", simplify = T)[, 1]
      }
      
      # mkdir
      if (T) {
        prefix <- paste0("58.DCSG_conversion/", data_type)
        
        dir <- paste0(c(
          "results/1.tab/", 
          "results/2.pic/", 
          "data/saved_data/"), prefix)
        mkdir_fun(dir)
      }
      
      # GO analysis: Count=ifelse(length(gene)>1000, 10, 5), P_cutoff=0.001, enrichment=0
      if (T) {
        file <- paste0("data/saved_data/58.go_", method, ".qs")
        
        if (! file.exists(file)) {
          all_conversions <- sort(unique(dat$label))
          go_summary <- lapply(all_conversions, function(type) {
            # type <- all_conversions[14]
            
            # info print
            if (T) {
              n <- which(all_conversions == type)
              total <- length(all_conversions)
              mess <- paste0("\tNow process GO: (", n, "/", total, ") ", type, "\n")
              cat(mess)
            }
            
            gene <- dat$id[dat$label == type]
            prefix <- paste0(prefix, "/58.GO.", method, ".", type)
            res <- one_step_go(gene, prefix=prefix, filetype="pdf", 
                               ont="ALL", Count=ifelse(length(gene)>1000, 10, 5), P_cutoff=0.001, enrichment=0, showCategory=8, # filter data for go barplot
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
      
      # simplify GO results
      if (F) {
        file <- paste0("data/saved_data/58.go_simplify_", method, ".qs")
        
        if (! file.exists(file)) {
          all_conversions <- sort(unique(dat$label))
          go_simplify <- lapply(all_conversions, function(type) {
            # type <- all_conversions[3]
            
            # default setting
            if (T) {
              filetype="pdf"
              ont="ALL"
              Count=10
              P_cutoff=0.001
              enrichment=0
              showCategory=8
              x="GeneRatio"
              orderBy="generatio"
              decreasing=T
            }
            
            # info print
            if (T) {
              n <- which(all_conversions == type)
              total <- length(all_conversions)
              mess <- paste0("\tNow process GO: (", n, "/", total, ") ", type, "\n")
              cat(mess)
            }
            
            # simplify GO enrichment analysis
            if (T) {
              go_simplify_file <- paste0("data/saved_data/", prefix, "/58.GO.simplify.", method, ".", type, ".qs")
              if (file.exists(go_simplify_file)) {
                go_dat_simplify <- qread(go_simplify_file, nthreads = 6)
              }
              if (! file.exists(go_simplify_file)) {
                # read raw go data
                if (T) {
                  go_file <- paste0("data/saved_data/", prefix, "/58.GO.", method, ".", type, ".qs")
                  go_dat <- qread(go_file, nthreads = 6)
                }
                
                # filter the go terms: pvalue<0.001 & Count>10
                if (T) {
                  go_dat@result <- go_dat@result[go_dat@result$pvalue<P_cutoff & go_dat@result$Count>Count, ]
                }
                
                # simplify the go terms
                if (nrow(go_dat@result)>0) {
                  go_dat_simplify <- clusterProfiler::simplify(go_dat, cutoff = 0.7, by = "p.adjust", select_fun = min)
                } else {
                  go_dat_simplify <- go_dat
                }
                
                qsave(go_dat_simplify, file = go_simplify_file, nthreads = 6)
              }
            }
            
            # data export (raw)
            if (nrow(go_dat_simplify@result)>0) {
              go_dat_simplify <- DOSE::setReadable(go_dat_simplify, OrgDb=org.Hs.eg.db, keyType='ENSEMBL')
              go_dat <- go_dat_simplify@result
              
              # format the data
              go_dat$generatio <- sapply(go_dat$GeneRatio, function(x) {eval(parse(text = x))})
              go_dat$bgratio <- sapply(go_dat$BgRatio, function(x) {eval(parse(text = x))})
              go_dat$enrichment_Level <- log2(go_dat$generatio/go_dat$bgratio)
              
              write.table(go_dat, file = paste0("results/1.tab/", prefix, "/58.GO.simplify.", method, ".", type, ".csv"), 
                          quote = F, sep = ",", col.names = T, row.names = F)
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
              
              file <- paste0("results/2.pic/", prefix, "/58.GO.simplify.", method, ".", type, ".pdf")
              ggsave(file, plot = p, width = 7, height = 10)
            }
            
            # add type info
            if (T) {
              go_dat$type <- type
            }
            
            return(go_dat)
          })
          go_simplify <- do.call(rbind, go_simplify)
          
          qsave(go_simplify, file = file, nthreads = 6)
        }
      }
    }
  }
}

# GO analysis summary
if (T) {
  for (data_type in c("gene")) {
    cat(paste0(data_type, ": \n"))
    
    for (method in c("chromIDEAS", "kmeans")) {
      cat(paste0("\t", method, ": \n"))
      
      # test data
      if (F) {
        data_type <- "gene"
        method <- "chromIDEAS"
      }
      
      # get simplify GO results
      if (T) {
        file <- paste0("data/saved_data/58.go_", method, ".qs")
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
        
        file <- paste0("results/2.pic/58.DCSG_conversion/", data_type, "/58.GO_cluster_GeneRatio_heatmap.", method, ".pdf")
        ggsave(filename = file, plot = p, width = 5, height = 8)
      }
      
      # group and word cloud setting
      if (T) {
        ht_opt$message = FALSE
        
        all_types <- sort(unique(go_filter$type))
        
        for (type in all_types) {
          # type <- all_types[7]
          
          cat(paste0("Now process (", which(type == all_types), "/", length(all_types), ") ", type, "\n"))
          
          file <- paste0("results/2.pic/58.DCSG_conversion/", data_type, "/58.GO_cluster_", method, ".", type, ".pdf")
          if (! file.exists(file)) {
            set.seed(799)
            
            go_similarity <- simplifyEnrichment::GO_similarity(go_filter$ID[go_filter$type == type], ont="BP", db="org.Hs.eg.db", measure="Rel")
            if (nrow(go_similarity) < 3) {
              cluster_res <- 1:nrow(go_similarity)
              names(cluster_res) <- rownames(go_similarity)
            } else {
              cluster_res <- simplifyEnrichment::cluster_by_pam(go_similarity, max_k = max(1, min(round(nrow(go_similarity)/20), 3)))
            }
            
            qsave(cluster_res, file = paste0("data/saved_data/58.DCSG_conversion/", data_type, "/58.GO_cluster_", method, ".", type, ".qs"), nthreads = 6)
            
            p <- simplifyEnrichment::ht_clusters(go_similarity, cluster_res)
            
            pdf(file = file, width = 8, height = 6)
            print(p)
            dev.off()
            
            rm(type, go_similarity, cluster_res, p)
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
        
        # sorting
        if (T) {
          ord <- hclust(dist(top20), method = "complete")$order
          ggdat$GO <- factor(ggdat$GO, levels = rownames(top20)[rev(ord)])
        }
        
        # ggplot
        if (T) {
          if (method == "chromIDEAS") {
            breaks <- seq(0, 45, 15)
          }
          if (method == "kmeans") {
            breaks <- seq(0, 60, 20)
          }
          
          p <- ggplot(ggdat) +
            geom_bar(aes(x=Number, y=GO, fill=GO), stat="identity", color="black") +
            scale_x_continuous(name="Number of GO Terms", breaks=breaks) +
            facet_grid(.~Conversion) +
            cowplot::theme_cowplot() +
            theme(panel.border = element_rect(color="black"), 
                  legend.position = "none",
                  strip.background = element_rect(fill=NA, color=NA))
        }
        
        file <- paste0("results/2.pic/58.DCSG_conversion/", data_type, "/58.top20_GO_ancestors_stat_barplot.", method, ".pdf")
        ggsave(filename = file, plot = p, width = 15, height = 6)
      }
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
        file <- paste0("data/saved_data/58.go_", method, ".qs")
        go_filter <- qread(file, nthreads = 6)
      }
      
      # filter the data
      if (T) {
        ggdat <- go_filter[go_filter$ONTOLOGY == "BP" & go_filter$type %in% c("3to1", "3to4"), ]
        ggdat <- lapply(c("3to1", "3to4"), function(subtype) {
          # subtype <- "3to4"
          subdat <- ggdat[ggdat$type == subtype, ]
          subdat <- subdat[order(subdat$generatio, decreasing = T), ]
          if (subtype == "3to1") {
            subdat <- head(subdat, 10)
          }
          if (subtype == "3to4") {
            subdat1 <- head(subdat, 8)
            subdat2 <- subdat[grepl("wnt", subdat$Description, ignore.case = T) & subdat$generatio > 0.05, ]
            subdat2 <- head(subdat2, 2)
            subdat <- rbind(subdat1, subdat2)
          }
          
          return(subdat)
        })
        names(ggdat) <- c("3to1", "3to4")
      }
      
      # ggplot2
      if (T) {
        for (i in c("3to1", "3to4")) {
          # i <- "3to1"
          p <- barplot_go(ggdat[[i]], x="GeneRatio")
          ggsave(filename = paste0("results/2.pic/58.DCSG_conversion/", data_type, "/58.conversion_", i, ".", method, ".pdf"), 
                plot = p, width = 7, height = 10)
        }
        
        rm(i, p)
      }
    }
  }
}

# KEGG analysis***************************
if (F) {
  for (data_type in c("gene")) {
    cat(paste0(data_type, ": \n"))
    
    for (method in c("chromIDEAS", "kmeans")) {
      cat(paste0("\t", method, ": \n"))
      
      # test data
      if (F) {
        data_type <- "gene"
        method <- "chromIDEAS"
      }
      
      # get CS conversion
      if (T) {
        file <- paste0("data/saved_data/58.", data_type, "_level_DCSG_representative_CS_Conversion(", method, ").qs")
        
        dat <- qread(file, nthreads = 6)
      }
      
      # format the data
      if (T) {
        dat$id <- str_split(dat$id, "[.]", simplify = T)[, 1]
      }
      
      # mkdir
      if (T) {
        prefix <- paste0("58.DCSG_conversion/", data_type)
        
        dir <- paste0(c(
          "results/1.tab/", 
          "results/2.pic/", 
          "data/saved_data/"), prefix)
        mkdir_fun(dir)
      }
      
      # KEGG analysis: Count=ifelse(length(gene)>1000, 10, 5), P_cutoff=0.001, enrichment=0
      if (T) {
        file <- paste0("data/saved_data/58.kegg_", method, ".qs")
        
        if (! file.exists(file)) {
          all_conversions <- sort(unique(dat$label))
          go_summary <- lapply(all_conversions, function(type) {
            # type <- all_conversions[14]
            
            # info print
            if (T) {
              n <- which(all_conversions == type)
              total <- length(all_conversions)
              mess <- paste0("\tNow process GO: (", n, "/", total, ") ", type, "\n")
              cat(mess)
            }
            
            gene <- dat$id[dat$label == type]
            prefix <- paste0(prefix, "/58.GO.", method, ".", type)
            res <- one_step_go(gene, prefix=prefix, filetype="pdf", 
                               ont="ALL", Count=ifelse(length(gene)>1000, 10, 5), P_cutoff=0.001, enrichment=0, showCategory=8, # filter data for go barplot
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
}

# KEGG analysis summary***************************
if (F) {
  for (data_type in c("gene")) {
    cat(paste0(data_type, ": \n"))
    
    for (method in c("chromIDEAS", "kmeans")) {
      cat(paste0("\t", method, ": \n"))
      
      # test data
      if (F) {
        data_type <- "gene"
        method <- "chromIDEAS"
      }
      
      # get simplify GO results
      if (T) {
        file <- paste0("data/saved_data/58.go_", method, ".qs")
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
        
        file <- paste0("results/2.pic/58.DCSG_conversion/", data_type, "/58.GO_cluster_GeneRatio_heatmap.", method, ".pdf")
        ggsave(filename = file, plot = p, width = 5, height = 8)
      }
      
      # group and word cloud setting
      if (T) {
        ht_opt$message = FALSE
        
        all_types <- sort(unique(go_filter$type))
        
        for (type in all_types) {
          # type <- all_types[7]
          
          cat(paste0("Now process (", which(type == all_types), "/", length(all_types), ") ", type, "\n"))
          
          file <- paste0("results/2.pic/58.DCSG_conversion/", data_type, "/58.GO_cluster_", method, ".", type, ".pdf")
          if (! file.exists(file)) {
            set.seed(799)
            
            go_similarity <- simplifyEnrichment::GO_similarity(go_filter$ID[go_filter$type == type], ont="BP", db="org.Hs.eg.db", measure="Rel")
            if (nrow(go_similarity) < 3) {
              cluster_res <- 1:nrow(go_similarity)
              names(cluster_res) <- rownames(go_similarity)
            } else {
              cluster_res <- simplifyEnrichment::cluster_by_pam(go_similarity, max_k = max(1, min(round(nrow(go_similarity)/20), 4)))
            }
            
            qsave(cluster_res, file = paste0("data/saved_data/58.DCSG_conversion/", data_type, "/58.GO_cluster_", method, ".", type, ".qs"), nthreads = 6)
            
            p <- simplifyEnrichment::ht_clusters(go_similarity, cluster_res)
            
            pdf(file = file, width = 8, height = 6)
            print(p)
            dev.off()
            
            rm(type, go_similarity, cluster_res, p)
          }
        }
      }
    }
  }
}




# 具体思路
## 对于2个细胞系中的同一个基因
## 先用CS序列计算出每个位点上的距离值
## 然后用Cluster替换CS，根据Cluster序列找出基因上不同的位点。
## 有不同Cluster位点的基因定义为DCSG
## 
## 
## DCSG的状态转换定义
## 根据CS序列上的距离差异，统计每种状态转换的距离贡献值。贡献度最高的状态转换则被定义为该DCSG的状态转换
## 如果存在多个状态转换同时最高，则优先考虑分布在TSS上的状态转换