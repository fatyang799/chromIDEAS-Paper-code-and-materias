# # load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(UpSetR)
  library(stringr)
  library(reshape2)
  library(qs)
}

# define function
if (T) {
  get_line_dat <- function(files, names) {
    # test 
    if (F) {
      files <- list.files(path = "data/raw_data/8.paper/PMC8175737/cpm", pattern = "tss.gz", full.names = T, recursive = T)
      names <- basename(files)
    }
    
    # read dat
    if (T) {
      dat <- lapply(files, function(x) {
        # x <- files[1]
        dat <- read.table(x, header = F, sep = "\t", fill = T, comment.char = "")
        
        return(dat)
      })
      dat <- do.call(rbind, dat)
    }
    
    # split the data
    if (T) {
      value_dat <- dat[seq(3, nrow(dat), 3), -c(1:2)]
      head_dat <- dat[1:2, -c(1:2)]
      
      rm(dat)
    }
    
    # merge the data
    if (T) {
      colnames(value_dat)[1:ncol(value_dat)] <- paste0("B", as.numeric(head_dat[2, 1:ncol(value_dat)]))
      value_dat$id <- names
      
      head_dat <- data.frame(t(head_dat))
      rownames(head_dat) <- paste0("B", as.numeric(head_dat$X2))
      
      head_dat <- apply(head_dat, 1, function(x) {
        ifelse(is.na(x[1]), NA, x)
      })
      
      head_dat <- head_dat[!is.na(head_dat)]
    }
    
    res <- list(value_dat=value_dat, head_dat=head_dat)
    
    return(res)
  }
  get_diff_dat <- function(dat, location, mk) {
    subdat_wt <- dat[dat$cell == "thp1" & dat$location == location & dat$mk == mk, ]
    subdat_ko <- dat[dat$cell == "thp1kdkdm4a" & dat$location == location & dat$mk == mk, ]
    
    subdat_ko <- subdat_ko[match(subdat_wt$Loc, subdat_ko$Loc), ]
    diff <- abs(subdat_wt$Signal - subdat_ko$Signal)
    
    return(diff)
  }
  
  get_gz_dat <- function(files, names) {
    # test 
    if (F) {
      files <- list.files(path = "data/raw_data/8.paper/PMC8175737/gzfiles/s3v2norm", pattern = "body.gz", full.names = T, recursive = T)
      names <- basename(files)
    }
    
    # define function to get auc value
    if (T) {
      auc_fun <- function(dat) {
        apply(dat, 1, function(row) {
          # row <- unlist(dat[1, ])
          value <- as.numeric(row)
          
          a <- value[-1]
          b <- value[-length(value)]
          
          sum((a+b)*1/2)
        })
      }
    }
    
    # read dat
    if (T) {
      dat <- lapply(1:length(files), function(x) {
        # x <- 17
        
        file <- files[x]
        name <- names[x]
        
        dat <- read.table(file, header = F, sep = "\t", fill = T, comment.char = "", skip = 1)
        dat <- dat[, -c(1:3, 5:6)]
        auc <- auc_fun(dat[, -1])
        
        dat <- data.frame(value = auc, 
                          tx = dat$V4, 
                          id = name)
        
        return(dat)
      })
      dat <- data.frame(do.call(rbind, dat))
    }
    
    # format the data
    if (T) {
      dat$tx <- str_split(dat$tx, "[.]", simplify = T)[, 1]
    }
    
    return(dat)
  }
  get_gz_diff_dat <- function(dat, ids, mk, method="-") {
    subdat_wt <- dat[dat$cell == "thp1" & dat$mk == mk, ]
    subdat_ko <- dat[dat$cell == "thp1kdkdm4a" & dat$mk == mk, ]
    
    subdat_wt <- subdat_wt[match(ids, subdat_wt$tx), ]
    subdat_ko <- subdat_ko[match(ids, subdat_ko$tx), ]
    
    if (method == "-") {
      diff <- subdat_ko$value - subdat_wt$value
    }
    if (method == "/") {
      diff <- subdat_ko$value / subdat_wt$value
    }
    if (method == "log2") {
      diff <- log2(subdat_ko$value/subdat_wt$value)
    }
    
    return(diff)
  }
  
  # define function
  if (T) {
    # script 58
    if (T) {
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
}

# plot signal distribution
if (T) {
  for (type in c("body", "tss")) {
    for (norm in c("s3v2norm", "cpm")) {
      # test data
      if (F) {
        type <- "tss"
        norm <- "cpm"
      }
      
      # get line data
      if (T) {
        files <- list.files(path = paste0("data/raw_data/8.paper/PMC8175737/", norm), pattern = paste0(type, ".mat"), full.names = T, recursive = T)
        value_dat <- get_line_dat(files, basename(files))
        
        labels <- value_dat$head_dat
        value_dat <- value_dat$value_dat
        
        rm(files)
      }
      
      # format the data
      if (T) {
        value_dat <- melt(value_dat, id.vars = "id", variable.name = "Loc", value.name = "Signal")
        value_dat$type <- "Line"
        
        value_dat$Signal <- as.numeric(value_dat$Signal)
        value_dat$Loc <- factor(value_dat$Loc, levels = paste0("B", sort(unique(as.numeric(gsub("B", "", value_dat$Loc))))))
        
        value_dat$cell <- str_split(value_dat$id, "_", simplify = T)[, 1]
        value_dat$mk <- str_split(value_dat$id, "_", simplify = T)[, 2]
        value_dat$location <- ifelse(grepl("_direct_up[.]", value_dat$id), "direct_up", 
                                     ifelse(grepl("_direct_down[.]", value_dat$id), "direct_down", 
                                            ifelse(grepl("_indirect_up[.]", value_dat$id), "indirect_up", 
                                                   ifelse(grepl("_indirect_down[.]", value_dat$id), "indirect_down", 
                                                          ifelse(grepl("non-regulated[.]", value_dat$id), "non-regulated", 
                                                                 ifelse(grepl("_bg[.]", value_dat$id), "bg", NA))))))
        
        value_dat$location <- factor(value_dat$location, levels = c("indirect_down", "direct_down", "non-regulated", "direct_up", "indirect_up", "bg"))
      }
      
      # filter data
      if (T) {
        dat <- value_dat
        dat <- value_dat[value_dat$mk != "kdm4a", ]
        
        rm(value_dat)
      }
      
      # color setting
      if (T) {
        library(RColorBrewer)
        colors <- brewer.pal(9, "Set1") 
        scales::show_col(colors, labels=T)
        colors <- c(colors[c(2,4)], colors[c(2,4)])
        names(colors) <- c("thp1kdkdm4a", "thp1", "H3K27me3", "H3K9me3")
      }
      
      # compare the signal in WT and kdm4a ko THP1 to confirm: ko kdm4a, the ability of kdm4a erase is loss
      if (T) {
        head(dat)
        
        # diff region dat
        if (T) {
          wt <- dat[dat$cell == "thp1", ]
          ko <- dat[dat$cell == "thp1kdkdm4a", ]
          
          if (identical(paste0(wt$Loc, wt$mk, wt$location), paste0(ko$Loc, ko$mk, ko$location))) {
            diff_dat <- wt[, c("Loc", "mk", "location")]
            diff_dat$ymin <- ifelse(wt$Signal > ko$Signal, ko$Signal, wt$Signal)
            diff_dat$ymax <- ifelse(wt$Signal > ko$Signal, wt$Signal, ko$Signal)
            diff_dat$group <- paste0(diff_dat$location, diff_dat$mk)
          }
          
          rm(wt, ko)
        }
        
        p <- ggplot() +
          geom_line(data = dat, aes(x=Loc, y=Signal, group=id, color=cell), linewidth=1) +
          geom_ribbon(data = diff_dat, aes(x=Loc, ymin=ymin, ymax=ymax, group=group, fill=mk), alpha = 0.3, color=NA) +
          facet_grid(mk~location, scales = "free") +
          scale_color_manual(values=colors) +
          scale_fill_manual(values=colors) +
          scale_x_discrete(name=NULL, breaks=names(labels), labels=labels) +
          facet_grid(mk~location, scales = "free") +
          cowplot::theme_cowplot() +
          theme(axis.title = element_text(size = rel(1.2)),
                axis.text = element_text(size = rel(1.2)),
                legend.text = element_text(size = rel(1.2)),
                legend.title = element_text(size = rel(1.2)),
                strip.text = element_text(size = rel(1.2)), 
                panel.border = element_rect(color="black"), 
                strip.background = element_rect(fill=NA, color=NA))
        
        # plotly::ggplotly(p)
        ggsave(filename = paste0("results/2.pic/64.diff_regions_signal_distribution_in_thp1_wt_vs_kdkdm4a.", type, ".", norm, ".pdf"), plot=p, width = 8, height=6)
      }
      
      rm(p, dat, diff_dat, labels)
    }
  }
}

# plot diff signal boxplot
if (T) {
  for (type in c("body", "tss")) {
    for (norm in c("s3v2norm", "cpm")) {
      # test data
      if (F) {
        type <- "tss"
        norm <- "s3v2norm"
      }
      
      # get line data
      if (T) {
        file <- paste0("data/saved_data/64.auc_dat_for_profile.", norm, ".", type, ".qs")
        if (file.exists(file)) {
          dat <- qread(file, nthreads = 6)
        }
        if (! file.exists(file)) {
          files <- list.files(path = paste0("data/raw_data/8.paper/PMC8175737/gzfiles/", norm), pattern = paste0(type, ".gz"), full.names = T, recursive = T)
          dat <- get_gz_dat(files, basename(files))
          qsave(dat, file, nthreads = 6)
          
          rm(files)
        }
      }
      
      # format the data
      if (T) {
        # ignore na value
        dat <- dat[!is.na(dat$value), ]
        
        dat$value <- as.numeric(dat$value)
        dat$cell <- str_split(dat$id, "_", simplify = T)[, 1]
        dat$mk <- str_split(dat$id, "_", simplify = T)[, 2]
        
        dat$location <- ifelse(grepl("_direct_up[.]", dat$id), "direct_up", 
                                     ifelse(grepl("_direct_down[.]", dat$id), "direct_down", 
                                            ifelse(grepl("_indirect_up[.]", dat$id), "indirect_up", 
                                                   ifelse(grepl("_indirect_down[.]", dat$id), "indirect_down", 
                                                          ifelse(grepl("non-regulated[.]", dat$id), "non-regulated", 
                                                                 ifelse(grepl("_bg[.]", dat$id), "bg", NA))))))
      }
      
      # filter data
      if (T) {
        dat <- dat[dat$mk != "kdm4a", ]
      }
      
      # color setting
      if (T) {
        library(RColorBrewer)
        colors <- brewer.pal(9, "Set1") 
        scales::show_col(colors, labels=T)
        colors <- c(colors[c(2,4)], colors[c(2,4)])
        names(colors) <- c("thp1kdkdm4a", "thp1", "H3K27me3", "H3K9me3")
      }
      
      # get diff dat
      if (T) {
        mks <- unique(dat$mk)
        locs <- unique(dat$location)
        ids <- unique(dat$tx)
        
        ggdat <- data.frame(id = rep(ids, each=length(mks)), 
                            mk = rep(mks, length(ids)))
        for (mk in mks) {
          # test data
          if (F) {
            mk <- mks[1]
          }
          
          ggdat[ggdat$mk == mk, "diff_substract"] <- get_gz_diff_dat(dat, ids=ggdat[ggdat$mk == mk, "id"], mk, method="-")
        }
        
        rm(ids)
      }
      
      # add region info
      if (T) {
        anno <- dat[, c("tx", "location")]
        anno <- anno[!duplicated(anno), ]
        
        ggdat$region <- anno$location[match(ggdat$id, anno$tx)]
        
        rm(anno)
      }
      
      # format the data
      if (T) {
        ggdat$region <- factor(ggdat$region, levels = c("indirect_down", "direct_down", "non-regulated", "direct_up", "indirect_up", "bg"))
        dat$location <- factor(dat$location, levels = c("indirect_down", "direct_down", "non-regulated", "direct_up", "indirect_up", "bg"))
        
        ggdat$log2diff <- log2(abs(ggdat$diff_substract)+1)
        ggdat$log2diff <- ifelse(ggdat$diff_substract>0, ggdat$log2diff, -ggdat$log2diff)
        
        dat$log2value <- log2(dat$value+1)
      }
      
      # stat
      if (F) {
        head(ggdat)
        
        stat <- data.frame(do.call(rbind, tapply(ggdat$diff_substract, paste0(ggdat$mk, "@", ggdat$region), summary)))
        stat$mk <- str_split(rownames(stat), "@", simplify = T)[, 1]
        stat$loc <- str_split(rownames(stat), "@", simplify = T)[, 2]
      }
      
      # boxplot
      if (T) {
        # compare
        if (T) {
          ref <- "direct_down"
          remain <- setdiff(locs, ref)
          my_comparisons <- lapply(remain, function(x) {
            c(ref, x)
          })
          
          rm(ref, remain, locs, mks)
        }
        
        # ggplot
        if (T) {
          # log2diff value distribution
          if (T) {
            head(ggdat)
            
            p <- ggplot(ggdat, aes(x=region, y=diff_substract, group=region, fill=region)) +
              geom_violin(scale="width") +
              geom_boxplot(fill=NA, width=0.2) +
              ggpubr::stat_compare_means(method="wilcox.test", paired = F, label="p.signif", comparisons = my_comparisons) +
              scale_y_continuous(name="diffSigal KO-WT", breaks = seq(0, 60, 5), limits = c(0, 60)) +
              facet_grid(mk~., scales = "free") +
              cowplot::theme_cowplot() +
              theme(axis.text = element_text(size = rel(0.7)), 
                    panel.border = element_rect(color="black"), 
                    legend.position = "none", 
                    strip.background = element_rect(fill=NA, color=NA))
            
            ggsave(filename = paste0("results/2.pic/64.diff_regions_diff_signal_after_kdkdm4a.", type, ".", norm, ".pdf"),
                   plot = p, width = 8, height = 6)
          }
          
          # raw data distribution
          if (T) {
            head(dat)
            
            p <- ggplot(dat, aes(x=location, y=log2value, fill=cell, group=interaction(location, cell))) +
              geom_violin(scale="width", position=position_dodge(width=0.8)) +
              geom_boxplot(fill=NA, width=0.2, position=position_dodge(width=0.8)) +
              ggpubr::stat_compare_means(aes(group=interaction(location, cell)), method = "wilcox.test", paired = T, label="p.signif") +
              facet_grid(mk~., scales = "free") +
              cowplot::theme_cowplot() +
              theme(axis.text = element_text(size = rel(0.7)), 
                    panel.border = element_rect(color="black"), 
                    legend.position = "none", 
                    strip.background = element_rect(fill=NA, color=NA))
            
            ggsave(filename = paste0("results/2.pic/64.diff_regions_log2signal_after_kdkdm4a.", type, ".", norm, ".pdf"),
                   plot = p, width = 8, height = 6)
          }
        }
      }
    }
  }
}

# rescue genes
if (T) {
  # get cd34 vs thp1 deg
  if (T) {
    file <- "results/1.tab/59.thp1_vs_cd34_mRNA_deg.csv"
    deg <- read.table(file, header = T, sep = ",", fill = T, comment.char = "")
    deg$id <- str_split(deg$id, "[.]", simplify = T)[, 1]
    
    deg_up <- deg$id[deg$state == "Up"]
  }
  
  # kdm4a direct down genes
  if (T) {
    kd <- read.table("data/raw_data/8.paper/PMC8175737/ko_deg.csv", header = T, sep = ",", fill = T, comment.char = "")
    kdm4a_erase_k9me3 <- read.table("results/1.tab/63.kdm4a_direct_down_genes.txt", header = F, sep = "\t", fill = T, comment.char = "")[, 1]
    kdm4a_erase_k9me3 <- kd[kd$Row.names %in% kdm4a_erase_k9me3, ]
  }
  
  # wnt related genes
  if (T) {
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    suppressPackageStartupMessages(library(AnnotationDbi))
    suppressPackageStartupMessages(library(GO.db))
    
    # get go id and info
    if (T) {
      keytypes(GO.db)
      
      go_ids <- keys(GO.db, keytype = "GOID")
      # DEFINITION 是详细介绍GO term的结果
      go_info <- select(GO.db, keys = go_ids, keytype = "GOID", columns = c("DEFINITION", "GOID", "ONTOLOGY", "TERM"))
      
      rm(go_ids)
    }
    
    # get wnt related terms and genes
    if (T) {
      wnt_go <- go_info[grepl("wnt", go_info$TERM, ignore.case = T), ]
      
      keytypes(org.Hs.eg.db)
      
      wnt_go_genes <- select(org.Hs.eg.db, keys = wnt_go$GOID, keytype = "GO", columns = c("ENSEMBL", "ALIAS"))
      
      wnt_go_genes <- wnt_go_genes[! is.na(wnt_go_genes$ENSEMBL), ]
      wnt_go_genes <- wnt_go_genes[!duplicated(wnt_go_genes$ENSEMBL), ]
      wnt_go_genes$name1 <- deg$name[match(wnt_go_genes$ENSEMBL, deg$id)]
      wnt_go_genes$name2 <- kd$external_gene_name[match(wnt_go_genes$ENSEMBL, kd$Row.names)]
      wnt_go_genes <- wnt_go_genes[grepl("wnt", wnt_go_genes$ALIAS, ignore.case = T)|
                                     grepl("wnt", wnt_go_genes$name1, ignore.case = T)|
                                     grepl("wnt", wnt_go_genes$name2, ignore.case = T), ]
    }
  }
  
  # transition genes 4to3/1
  if (T) {
    file <- paste0("data/saved_data/58.gene_level_DCSG_representative_CS_Conversion(chromIDEAS).qs")
    chromideas <- qread(file, nthreads = 6)
    chromideas <- chromideas$id[chromideas$label %in% c("3to1", "3to4")]
    chromideas <- str_split(chromideas, "[.]", simplify = T)[, 1]
  }
  
  # overlap genes
  if (T) {
    # 817
    length(intersect(deg_up, chromideas))
    length(kdm4a_erase_k9me3$Row.names)
    
    length(intersect(intersect(deg_up, chromideas), kdm4a_erase_k9me3$Row.names))
    
    if (F) {
      ids <- intersect(intersect(deg_up, chromideas), kdm4a_erase_k9me3)
      go <- clusterProfiler::enrichGO(gene = ids,
                                      OrgDb = org.Hs.eg.db,
                                      keyType = 'ENSEMBL',
                                      ont = "ALL",
                                      pAdjustMethod = "BH",
                                      pvalueCutoff = 1,
                                      qvalueCutoff = 1)
      go <- DOSE::setReadable(go, OrgDb=org.Hs.eg.db, keyType='ENSEMBL')
      go_dat <- go@result
      
      # format the data
      go_dat$generatio <- sapply(go_dat$GeneRatio, function(x) {eval(parse(text = x))})
      go_dat$bgratio <- sapply(go_dat$BgRatio, function(x) {eval(parse(text = x))})
      go_dat$enrichment_Level <- log2(go_dat$generatio/go_dat$bgratio)
      go_dat <- go_dat[order(go_dat$generatio, decreasing = T), ]
    }
    
    ids <- intersect(intersect(intersect(deg_up, kdm4a_erase_k9me3$Row.names), wnt_go_genes$ENSEMBL), chromideas)
    
    ids <- data.frame(id = ids)
    ids$kd <- kdm4a_erase_k9me3$log2FoldChange[match(ids$id, kdm4a_erase_k9me3$Row.names)]
    ids$deg_cd34 <- deg$logFC[match(ids$id, deg$id)]
    ids$name <- deg$name[match(ids$id, deg$id)]
    
    ggplot(ids) +
      geom_point(aes(x=kd, y=deg_cd34)) +
      geom_text(aes(x=kd, y=deg_cd34, label=name))
  }
}

# venn plot
if (T) {
  library(ggVennDiagram)
  
  # all DCSG
  if (T) {
    venn_dat <- list(
      deg_up = deg_up,
      kdm4a_erase_k9me3 = kdm4a_erase_k9me3$Row.names, 
      wnt_go_genes = wnt_go_genes$ENSEMBL
    )
    
    # color setting
    if (T) {
      library(RColorBrewer)
      colors <- brewer.pal(9, "Set1") 
      
      colors <- colors[1:3]
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
  }
  
  file <- paste0("results/2.pic/64.overlap_kdm4a_eraser_genes_and_deg_Up.pdf")
  ggsave(filename = file, plot = p, width = 8, height = 6)
}

