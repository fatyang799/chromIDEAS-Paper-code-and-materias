# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(UpSetR)
  library(stringr)
  library(reshape2)
  library(qs)
  suppressMessages(library(GO.db))
  suppressMessages(library(org.Hs.eg.db))
  suppressMessages(library(AnnotationDbi))
  suppressPackageStartupMessages(library(clusterProfiler))
}

cell1 <- "thp1"
cell2 <- "cd34"
data_type <- "gene"

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

# get KD KDM4A rna-seq data
if (T) {
  # read raw data
  if (T) {
    file <- "data/raw_data/8.paper/PMC8175737/ko_deg.csv"
    # file <- "results/1.tab/71.thp1_kdm4a_mRNA_deg.csv"
    kd <- read.table(file, header = T, sep = ",", fill = T, comment.char = "")
    # kd$id <- str_split(kd$id, "[.]", simplify = T)[, 1]
  }
  
  # format the data: padj<0.05 & log2FC>0.5 (define as the paper)
  if (T) {
    colnames(kd)[1] <- "id"
    kd <- kd[!is.na(kd$log2FoldChange), ]
    kd <- kd[!is.na(kd$padj), ]
    kd$state <- ifelse(kd$padj<0.05 & abs(kd$log2FoldChange)>0.5,
                       ifelse(kd$log2FoldChange>0, "Up", "Down"), "NotSig")

    kd$logPvalue <- (-log10(kd$padj))
    # kd$logPvalue <- (-log10(kd$adj.P.Val))
  }
}

# get deg with chromIDEAS
if (T) {
  # 3to4_3to1_up_deg overlap KD_KDM4A_down_genes
  if (T) {
    dat <- qread("data/saved_data/60.dcsg_deg_merged_info.qs", nthreads = 6)
    dat$id <- str_split(dat$id, "[.]", simplify = T)[, 1]
    
    genes <- intersect(
      dat$id[dat$state == "Up" & (dat$CS_conversion == "3to4" | dat$CS_conversion == "3to1")], 
      kd$id[kd$state == "Down"]
    )
  }
  
  # go for overlap genes
  if (T) {
    file <- "data/saved_data/62.3to4_3to1_up_deg_overlap_KD_KDM4A_down_go.qs"
    if (file.exists(file)) {
      go <- qread(file, nthreads = 6)
    }
    if (! file.exists(file)) {
      go <- enrichGO(gene = genes,
                     OrgDb = org.Hs.eg.db,
                     keyType = 'ENSEMBL',
                     ont = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 1,
                     qvalueCutoff = 1)
      qsave(go, file = file, nthreads = 6)
    }
  }
  
  # format the results
  if (T) {
    gores <- go@result
    
    # format the data
    gores$generatio <- sapply(gores$GeneRatio, function(x) {eval(parse(text = x))})
    gores$bgratio <- sapply(gores$BgRatio, function(x) {eval(parse(text = x))})
    gores$enrichment_Level <- log2(gores$generatio/gores$bgratio)
    gores <- gores[order(gores[, "generatio"], decreasing = T), ]
    gores <- gores[gores$ONTOLOGY == "BP", ]
  }
  
  # barplot show top10 terms
  if (T) {
    gores_filter <- head(gores, 10)
    p <- barplot_go(gores_filter) +
      theme(axis.title = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1.2)),
            legend.text = element_text(size = rel(1.2)),
            legend.title = element_text(size = rel(1.2)),
            strip.text = element_text(size = rel(1.2)))
  }
  
  ggsave(filename = "results/2.pic/62.3to4_3to1_up_deg_overlap_KD_KDM4A_down_genes.pdf", 
         plot = p, width = 10, height = 10)
  
  rm(p, go, gores, gores_filter, dat)
}

# valcano plot
if (T) {
  # label genes
  if (T) {
    # wnt go related genes (Up)
    if (T) {
      wnt_gene1 <- kd$id[kd$state == "Down" & grepl("wnt", kd$external_gene_name, ignore.case = T)]
      wnt_gene2 <- kd$id[kd$state == "Down" & grepl("tcf", kd$external_gene_name, ignore.case = T)]
      wnt_gene3 <- kd$id[match(c("KDM4A", "PAF1"), kd$external_gene_name)]
    }
    
    label <- kd[kd$id %in% c(wnt_gene1, wnt_gene2, wnt_gene3), ]
  }
  
  # color setting
  if (T) {
    library(RColorBrewer)
    colors <- brewer.pal(9, "Set1") 
    
    colors <- c(colors[1:2], "grey")
    names(colors) <- c("Up", "Down", "NotSig")
  }
  
  # valcano plot
  if (T) {
    head(kd)
    
    p <- ggplot() +
      geom_point(data=kd, aes(x=log2FoldChange, y=logPvalue, color=state), shape=16, alpha=0.8, size=3) +
      geom_point(data=label, aes(x=log2FoldChange, y=logPvalue), shape=21, fill=NA, color="black", alpha=1, size=3) +
      ggrepel::geom_text_repel(data=label, aes(log2FoldChange, y=logPvalue, label=external_gene_name), color="black", min.segment.length=0.1, max.overlaps=15) +
      geom_hline(yintercept = (-log10(0.05)), linetype=2) +
      geom_vline(xintercept = c(-0.5, 0.5), linetype=2) +
      scale_color_manual(values=colors) +
      scale_x_continuous(name="Log2FC", breaks = sort(c(seq(-15, 15, 5), -1, 1))) +
      ylab("-log10(FDR)") +
      cowplot::theme_cowplot()
  }
  
  file <- "results/2.pic/62.kd_kdm4a_deg_volcano.pdf"
  ggsave(p, file=file, width = 8, height = 5)
  
  rm(p, label)
}

# KO KDM4A deg gene GO
if (T) {
  file <- "data/saved_data/62.ko_kmd4a_gene_go.qs"
  
  if (! file.exists(file)) {
    deg_go <- lapply(c("Up", "Down"), function(type) {
      # type <- "Up"
      
      gene <- kd$id[kd$state == type]
      go <- one_step_go(gene, prefix=paste0("62.ko_kmd4a_", type, "_gene_go"), filetype="pdf", 
                        ont="BP", Count=ifelse(length(gene)>1000, 10, 5), P_cutoff=0.001, enrichment=0, showCategory=15, 
                        x="GeneRatio", orderBy="GeneRatio", decreasing=T)
      go <- go$full_list
      
      go$type <- type
      
      return(go)
    })
    deg_go <- do.call(rbind, deg_go)
    
    qsave(deg_go, file, nthreads = 6)
  }
  if (file.exists(file)) {
    deg_go <- qread(file, nthreads = 6)
  }
}

# GSEA
if (T) {
  file <- "data/saved_data/62.ko_kmd4a_gene_go_gsea.qs"
  if (! file.exists(file)) {
    genelist <- kd$log2FoldChange
    names(genelist) <- kd$id
    genelist <- sort(genelist, decreasing = T)
    
    gogsea <- gseGO(genelist, ont="BP", OrgDb = org.Hs.eg.db,
                    keyType = 'ENSEMBL', minGSSize = 5, maxGSSize = 2000,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 1)
    
    qsave(gogsea, file = file, nthreads = 6)
  }
  if (file.exists(file)) {
    gogsea <- qread(file, nthreads = 6)
  }
  
  # ggplot
  if (T) {
    gsea <- gogsea@result
    
    p <- enrichplot::gseaplot2(gogsea,
                               geneSetID = gsea$ID[gsea$pvalue<0.05 & grepl("wnt", gsea$Description, ignore.case = T)],
                               pvalue_table = T, rel_heights = c(0.8,0.2,0.5), color = c("#E64B35FF"))
  }
  
  file <- "results/2.pic/62.ko_kmd4a_gene_go_gsea.pdf"
  ggsave(p, file=file, width = 8, height = 5)
}

# KEGG
if (F) {
  # install.packages("KEGG.db_1.0.tar.gz", repos = NULL)
  
  keggfun <- function(genes) {
    # test data
    if (F) {
      genes <- kd$id[kd$state == "Down"]
    }
    
    suppressMessages(library(KEGG.db))
    
    # id conversion
    if (T) {
      id <- bitr(geneID=genes, fromType = "ENSEMBL", toType = "ENTREZID", org.Hs.eg.db)
    }
    
    # kegg enrichment analysis
    if (T) {
      kegg <- enrichKEGG(gene = id$ENTREZID,
                         organism = "hsa",
                         keyType = 'kegg',
                         pAdjustMethod = "BH",
                         pvalueCutoff = 1,
                         qvalueCutoff = 1, 
                         minGSSize = 10,
                         maxGSSize = 5000, 
                         use_internal_data = T)
    }
    
    return(kegg)
  }
  downkess <- keggfun(kd$id[kd$state == "Down"])
  downkess <- downkess@result
  
  upkess <- keggfun(kd$id[kd$state == "Up"])
  upkess <- upkess@result
}

# get wnt related GO from database
if (F) {
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  suppressPackageStartupMessages(library(AnnotationDbi))
  suppressPackageStartupMessages(library(GO.db))
  
  # get go id and info
  if (T) {
    keytypes(GO.db)
    
    go_ids <- keys(GO.db, keytype = "GOID")
    # DEFINITION 是详细介绍GO term的结果
    go_info <- select(GO.db, keys = go_ids, keytype = "GOID", columns = c("DEFINITION", "GOID", "ONTOLOGY", "TERM"))
  }
  
  # get wnt related terms and genes
  if (T) {
    wnt_go <- go_info[grepl("wnt", go_info$TERM, ignore.case = T) & (! grepl("negative", go_info$TERM, ignore.case = T)), ]
    
    keytypes(org.Hs.eg.db)
    
    wnt_go_genes <- select(org.Hs.eg.db, keys = wnt_go$GOID, keytype = "GO", columns = c("ENSEMBL", "ALIAS"))
    
    wnt_go_genes <- wnt_go_genes[! is.na(wnt_go_genes$ENSEMBL), ]
    wnt_go_genes <- wnt_go_genes[!duplicated(wnt_go_genes$ENSEMBL), ]
  }
  
  # output the wnt related gene id
  if (T) {
    file <- "results/1.tab/62.wnt_related_GO_genes.txt"
    write(wnt_go_genes$ENSEMBL, file = file, sep = "\t")
  }
  
  # 提取wnt通路相关基因坐标
  if (F) {
    # get gtf all gene
    if (T) {
      txdb <- loadDb("data/raw_data/gencode.v40.annotation.sqlite")
      # only focus on the chr1-22.XY
      seqlevels(txdb) <- c(paste0("chr", 1:22), "chrX", "chrY")
      
      tx <- genes(txdb)
    }
    
    # get all gene 
    
    
  }
}

# get wnt related genes from GO database
if (F) {
  go_terms <- as.data.frame(GO.db::GOTERM)
  wnt_go <- go_terms[grepl("wnt", go_terms$Term, ignore.case = T) & (! grepl("negative", go_terms$Term, ignore.case = T)), ]
  
  wnt_related_genes <- select(org.Hs.eg.db, 
                              keys = wnt_go$go_id, 
                              columns = c("ENSEMBL", "GENETYPE", "GENENAME", "ALIAS"), 
                              keytype = "GO")
  wnt_related_genes <- wnt_related_genes[! is.na(wnt_related_genes$ENSEMBL), ]
  
  # 405
  wnt_related_genes <- wnt_related_genes[! duplicated(wnt_related_genes$ENSEMBL), ]
  
  rm(go_terms, wnt_go)
}