# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(stringr)
  library(ggplot2)
  library(reshape2)
  suppressMessages(library(GenomicFeatures))
  suppressMessages(library(qs))
}

# load the txdb file
if (T) {
  file <- "data/raw_data/gencode.v40.annotation.sqlite"
  txdb <- loadDb(file)
  txdb
}

# load the peaks regions
if (T) {
  file <- "data/raw_data/8.paper/kdm4a_thp1_peaks.narrowPeak"
  peak <- read.table(file,header = F,sep = "\t")
  colnames(peak) <- c("chr","start","end","name","score","strand","signalValue","pValue","qValue","summit")
  peak <- peak[peak$chr %in% paste0("chr",c(1:22,"X","Y")),]
  peak <- GRanges(seqnames = peak$chr,
                  ranges = IRanges(start = peak$start+1, end = peak$end),
                  # strand = peak$strand,
                  name = peak$name,
                  score = peak$score,
                  signalValue = peak$signalValue,
                  pValue = peak$pValue,
                  qValue = peak$qValue,
                  summit = peak$summit)
  peak
}

# load the blacklist
if (T) {
  file <- "data/raw_data/2.states/hg38-blacklist.v2.bed"
  
  bl <- read.table(file, header = F, sep = "\t", fill = T, comment.char = "")
  bl <- GRanges(seqnames = bl$V1,
                ranges = IRanges(start = bl$V2+1, end = bl$V3))
  bl
}

# target region setting [TSS,TES]
if (T) {
  # 61544
  gene <- genes(txdb)
  
  # findOverlaps with bl
  if (T) {
    hit <- findOverlaps(query=bl, subject=gene, type="any")
    
    overlap_gene <- subjectHits(hit)
    bl_genes <- gene$gene_id[unique(overlap_gene)]
  }
  
  # remove gene overlap with bl (59413)
  gene <- gene[! gene$gene_id %in% bl_genes]
  
  rm(bl, overlap_gene, bl_genes, hit)
}

# findOverlaps
if (T) {
  hit <- findOverlaps(query=peak, subject=gene,
                      minoverlap=10, type="any")
  
  peak_info <- queryHits(hit)
  peak_info <- data.frame(peak[peak_info])
  
  overlap_gene <- subjectHits(hit)
  peak_info$overlap_gene_ID <- gene$gene_id[overlap_gene]
  
  rm(hit, overlap_gene, peak)
}

# stat
if (T) {
  stat <- as.data.frame.table(table(table(peak_info$overlap_gene_ID)))
  colnames(stat) <- c("N_pk_per_Gene", "N_Gene")
  
  p <- ggplot(stat) +
    geom_bar(aes(x=N_pk_per_Gene, y=N_Gene), stat = "identity") +
    theme_bw()
  print(p)
  
  rm(stat)
}

# format the data
if (T) {
  # 22060
  target_genes <- unique(peak_info$overlap_gene_ID)
  target_genes <- str_split(target_genes, "[.]", simplify = T)[, 1]
}

# overlap between KO KDM4A DEG and KDM4A ChIP binding gene
if (T) {
  # get kd deg
  if (T) {
    file <- "data/raw_data/8.paper/PMC8175737/ko_deg.csv"
    kd <- read.table(file, header = T, sep = ",", fill = T, comment.char = "")
    colnames(kd)[1] <- "id"
    kd <- kd[!is.na(kd$log2FoldChange), ]
    kd <- kd[!is.na(kd$padj), ]
    kd$state <- ifelse(kd$padj<0.05 & abs(kd$log2FoldChange)>0.5, 
                       ifelse(kd$log2FoldChange>0, "Up", "Down"), "NotSig")
  }
  
  library(ggVennDiagram)
  
  # all DCSG
  if (T) {
    venn_dat <- list(
      KO_KDM4A_Up = kd$id[kd$state == "Up"],
      KO_KDM4A_Down = kd$id[kd$state == "Down"],
      KDM4A_Binding = target_genes, 
      All_genes = str_split(gene$gene_id, "[.]", simplify = T)[, 1]
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
  }
  
  file <- paste0("results/2.pic/63.overlap_ko_kdm4a_deg_and_kdm4a_binging_pk.pdf")
  ggsave(filename = file, plot = p, width = 8, height = 6)
  
  rm(p, venn_dat, file)
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

# region type output
if (T) {
  # get all genes
  if (T) {
    all_genes <- gene$gene_id
    all_genes <- str_split(all_genes, "[.]", simplify = T)[, 1]
  }
  
  # go for all type genes
  if (T) {
    file <- "data/saved_data/63.go_for_direct_and_indirected_genes_related_KDM4A.qs"
    
    if (! file.exists(file)) {
      gores <- lapply(c("direct_up", "direct_down", "indirect_up", "indirect_down", "non-regulated", "bg"), function(type) {
        # type <- "direct_up"
        
        # get up and down regulated genes
        if (T) {
          up_deg <- kd$id[kd$state == "Up"]
          down_deg <- kd$id[kd$state == "Down"]
          all_deg <- kd$id[kd$state != "NotSig"]
        }
        
        # get target genes
        if (T) {
          if (type == "direct_up") {
            dat <- intersect(target_genes, up_deg)
          }
          if (type == "direct_down") {
            dat <- intersect(target_genes, down_deg)
          }
          if (type == "indirect_up") {
            dat <- setdiff(up_deg, target_genes)
          }
          if (type == "indirect_down") {
            dat <- setdiff(down_deg, target_genes)
          }
          if (type == "non-regulated") {
            dat <- setdiff(target_genes, all_deg)
          }
          if (type == "bg") {
            dat <- setdiff(all_genes, unique(c(all_deg, target_genes)))
          }
        }
        
        # output gene list
        if (T) {
          write.table(dat, file = paste0("results/1.tab/63.kdm4a_", type, "_genes.txt"), 
                      quote = F, sep = "\t", col.names = F, row.names = F)
        }
        
        # go enrichment
        if (T) {
          go <- one_step_go(dat, prefix=paste0("63.kdm4a_", type, "_go"), filetype="pdf", 
                            ont="BP", Count=ifelse(length(gene)>1000, 10, 5), P_cutoff=0.001, enrichment=0, showCategory=10, 
                            x="GeneRatio", orderBy="GeneRatio", decreasing=T)
          gores <- go$full_list
          
          gores$type <- type
        }
        
        return(gores)
      })
      gores <- data.frame(do.call(rbind, gores))
      
      qsave(gores, file, nthreads = 6)
    }
    if (file.exists(file)) {
      gores <- qread(file, nthreads = 6)
    }
  }
  
  # go for directed down
  if (T) {
    go <- gores[gores$type == "direct_down" & gores$pvalue<0.05 & gores$ONTOLOGY == "BP", ]
    go <- head(go, 10)
    
    p <- barplot_go(go, x="GeneRatio")
    
    file <- "results/2.pic/63.kdm4a_direct_down_go_top10_P0.05.pdf"
    ggsave(file, plot = p, width = 7, height = 10)
  }
}

