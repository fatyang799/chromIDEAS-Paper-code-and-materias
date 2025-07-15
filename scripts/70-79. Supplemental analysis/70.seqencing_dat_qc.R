# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(stringr)
  library(ggplot2)
  library(reshape2)
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(circlize))
}

# function
if (T) {
  prepare_dat <- function(files, group) {
    # read data
    if (T) {
      dat <- lapply(files, function(file) {
        read.table(file, header = F, sep = "\t", fill = T, comment.char = "")
      })
      dat <- do.call(rbind, dat)
    }
    
    # split the data
    if (T) {
      value_dat <- dat[seq(3, nrow(dat), 3), -c(1:2)]
      head_dat <- dat[1:2, -c(1:2)]
      
      rm(dat)
      
      # value_dat
      if (T) {
        colnames(value_dat) <- paste0("B", as.numeric(head_dat[2, ]))
        rownames(value_dat) <- 1:nrow(value_dat)
      }
      
      # head_dat
      if (T) {
        head_dat <- data.frame(t(head_dat))
        rownames(head_dat) <- paste0("B", as.numeric(head_dat$X2))
        
        head_dat <- apply(head_dat, 1, function(x) {
          ifelse(is.na(x[1]), NA, x)
        })
        
        head_dat <- head_dat[!is.na(head_dat)]
      }
    }
    
    # format the data
    if (T) {
      value_dat <- data.frame(
        sapply(value_dat, as.numeric)
      )
    }
    
    # merge the data
    if (T) {
      profile <- cbind(value_dat, group)
    }
    
    # wide2long
    if (T) {
      ggdat <- melt(profile, id.vars = colnames(group), measure.vars = colnames(value_dat), variable.name = "Loc", value.name = "Signal")
      ggdat$Loc <- factor(ggdat$Loc, levels = paste0("B", sort(as.numeric(gsub("B", "", unique(ggdat$Loc))))))
    }
    
    res <- list(ggdat=ggdat, head_dat=head_dat)
    
    return(res)
  }
  
  plot_profile <- function(dat, color_n="rep", group_n="id", facet1_n="cell", facet2_n="mk", color=NA, show_legend=F, assist=F) {
    # test data
    if (F) {
      color_n <- "rep"
      group_n <- "id"
      facet1_n <- "cell"
      facet2_n <- "mk"
      show_legend <- F
      assist <- F
    }
    
    # split data
    if (T) {
      ggdat <- dat$ggdat
      head_dat <- dat$head_dat
    }
    
    # format the data
    if (T) {
      colnames(ggdat)[colnames(ggdat) == color_n] <- "color"
      colnames(ggdat)[colnames(ggdat) == group_n] <- "group"
      colnames(ggdat)[colnames(ggdat) == facet1_n] <- "facet1"
      colnames(ggdat)[colnames(ggdat) == facet2_n] <- "facet2"
    }
    
    # color setting
    if (sum(is.na(color)) == 0) {
      if (length(color) != length(unique(ggdat$color))) {
        stop(paste0("The number of colors you provide is not equal to the number of type for ", color_n, ", err1\n"))
      }
    }
    
    # profile
    if (T) {
      head(ggdat)
      
      p <- ggplot(ggdat) +
        geom_line(aes(x=Loc, y=Signal, group=group, color=color), linewidth=1, alpha=0.6) +
        facet_grid(facet1~facet2, scales = "free") +
        scale_x_discrete(name=NULL, breaks=names(head_dat), labels=head_dat) +
        cowplot::theme_cowplot() +
        theme(panel.border = element_rect(color="black"), 
              strip.background = element_rect(fill=NA, color=NA))
      
      # color settingg
      if (sum(is.na(color)) == 0) {
        p <- p + scale_color_manual(values=colors)
      }
      
      # show assistant line
      if (assist) {
        p <- p + geom_vline(xintercept = c(names(head_dat)), linetype=2, alpha=0.6, color="grey")
      }
      
      # show legend
      if (! show_legend) {
        p <- p + theme(legend.position = "none")
      }
    }
    
    return(p)
  }
  plot_average_barplot <- function(ggdat, value="total_passed", average="mean", ylab="value") {
    # test data
    if (F) {
      value <- "total_passed"
      average <- "mean"
    }
    
    # calculate the average value for barplot
    if (T) {
      bar_dat <- tapply(ggdat[, value], paste0(ggdat$cell, "@", ggdat$mk), get(average))
      bar_dat <- data.frame(
        cell = str_split(names(bar_dat), "@", simplify = T)[, 1], 
        mk = str_split(names(bar_dat), "@", simplify = T)[, 2], 
        value = as.numeric(bar_dat)
      )
    }
    
    # format the data
    if (T) {
      colnames(ggdat)[colnames(ggdat) == value] <- "value"
    }
    
    # ggplot
    if (T) {
      p <- ggplot() +
        geom_bar(data=bar_dat, aes(x=mk, y=value, fill=mk), stat = "identity", width=0.8) +
        geom_jitter(data=ggdat, aes(x=mk, y=value), width=0.2, size=3, alpha=0.5) +
        xlab(NULL) +
        ylab(ylab) +
        facet_grid(cell~.) +
        cowplot::theme_cowplot() +
        theme(panel.border = element_rect(color="black"), 
              legend.position = "none", 
              strip.background = element_rect(fill=NA, color=NA))
    }
    
    return(p)
  }
}

# ggplot for profile
if (T) {
  # norm data
  if (T) {
    # get files
    if (T) {
      files <- list.files(path = "data/raw_data/2.states", pattern = "mat", recursive = T, full.names = T)
    }
    
    # prepare data
    if (T) {
      # get files
      if (T) {
        tss_files <- files[grepl("tss", files, ignore.case = T)]
        body_files <- files[grepl("body", files, ignore.case = T)]
      }
      
      # get metadata
      if (T) {
        group_tss <- data.frame(str_split(basename(tss_files), "_", simplify = T)[, 1:2])
        colnames(group_tss) <- c("cell", "id")
        group_tss$mk <- str_split(group_tss$id, "[.]", simplify = T)[, 2]
        group_tss$rep <- str_split(group_tss$id, "[.]", simplify = T)[, 1]
        
        group_body <- data.frame(str_split(basename(body_files), "_", simplify = T)[, 1:2])
        colnames(group_body) <- c("cell", "id")
        group_body$mk <- str_split(group_body$id, "[.]", simplify = T)[, 2]
        group_body$rep <- str_split(group_body$id, "[.]", simplify = T)[, 1]
      }
      
      # prepare data
      if (T) {
        tss_dat <- prepare_dat(tss_files, group_tss)
        body_dat <- prepare_dat(body_files, group_body)
      }
      
      rm(group_tss, group_body)
    }
    
    p <- plot_profile(tss_dat, color_n="rep", group_n="id", facet1_n="cell", facet2_n="mk", color=NA, show_legend=F, assist=F)
    file <- "results/2.pic/70.norm_dat_profile_TSS.pdf"
    ggsave(filename = file, plot = p, width = 8, height = 5)
    
    p <- plot_profile(body_dat, color_n="rep", group_n="id", facet1_n="cell", facet2_n="mk", color=NA, show_legend=F, assist=F)
    file <- "results/2.pic/70.norm_dat_profile_Body.pdf"
    ggsave(filename = file, plot = p, width = 8, height = 5)
  }
  
  # raw data
  if (T) {
    # get files
    if (T) {
      files <- list.files(path = "data/raw_data/3.raw_fq_qc/", pattern = "mat", recursive = T, full.names = T)
    }
    
    # prepare data
    if (T) {
      # get files
      if (T) {
        tss_files <- files[grepl("tss", files, ignore.case = T)]
        body_files <- files[grepl("body", files, ignore.case = T)]
      }
      
      # get metadata
      if (T) {
        group_tss <- data.frame(str_split(basename(tss_files), "_", simplify = T)[, 1:3])
        colnames(group_tss) <- c("cell", "mk", "rep")
        group_tss$id <- paste0(group_tss$rep, ".", group_tss$mk)
        group_tss <- group_tss[, c("cell", "id", "mk", "rep")]
        
        group_body <- data.frame(str_split(basename(body_files), "_", simplify = T)[, 1:3])
        colnames(group_body) <- c("cell", "mk", "rep")
        group_body$id <- paste0(group_body$rep, ".", group_body$mk)
        group_body <- group_tss[, c("cell", "id", "mk", "rep")]
      }
      
      # prepare data
      if (T) {
        tss_dat <- prepare_dat(tss_files, group_tss)
        body_dat <- prepare_dat(body_files, group_body)
      }
      
      rm(group_tss, group_body)
    }
    
    p <- plot_profile(tss_dat, color_n="rep", group_n="id", facet1_n="cell", facet2_n="mk", color=NA, show_legend=F, assist=F)
    file <- "results/2.pic/70.raw_dat_profile_TSS.pdf"
    ggsave(filename = file, plot = p, width = 8, height = 5)
    
    p <- plot_profile(body_dat, color_n="rep", group_n="id", facet1_n="cell", facet2_n="mk", color=NA, show_legend=F, assist=F)
    file <- "results/2.pic/70.raw_dat_profile_Body.pdf"
    ggsave(filename = file, plot = p, width = 8, height = 5)
  }
}

# complexheatmap for correlation
if (T) {
  # read raw data
  if (T) {
    file <- "data/raw_data/3.raw_fq_qc/pearson_correlation_bs200.tab"
    raw_cor <- read.table(file, header = T, sep = "\t", fill = T, comment.char = "", skip = 1, row.names = 1)
  }
  
  # read norm data
  if (T) {
    file <- "data/raw_data/2.states/pearson_correlation_bs200.tab"
    norm_cor <- read.table(file, header = T, sep = "\t", fill = T, comment.char = "", skip = 1, row.names = 1)
  }
  
  # format the data
  if (T) {
    if (! identical(rownames(raw_cor), colnames(raw_cor))) {
      stop(paste0("The order of mat is not identical, please check, err1"))
    }
    
    head(rownames(raw_cor))
    
    if (identical(rownames(norm_cor), colnames(norm_cor))) {
      name <- rownames(norm_cor)
      head(name)
      name <- gsub("[.]", "_", name)
      name <- str_split(name, "_", simplify = T)
      name <- paste(name[, 1], name[, 3], name[, 2], sep = "_")
      
      rownames(norm_cor) <- colnames(norm_cor) <- name
    }
  }
  
  # make the order identical
  if (T) {
    raw_cor <- raw_cor[name, ]
    raw_cor <- raw_cor[, name]
  }
  
  # annotation dat
  if (T) {
    label <- data.frame(str_split(name, "_", simplify = T))
    head(label)
    
    colnames(label) <- c("cell", "mk", "rep")
  }
  
  # color setting
  if (T) {
    library(RColorBrewer)
    display.brewer.all()
    colors <- brewer.pal(9, "Set1") 
    
    col_cell <- colors[1:2]
    names(col_cell) <- c("cd34", "thp1")
    
    col_mk <- colors[1:8]
    sort(unique(label$mk))
    names(col_mk) <- c("H3K4me3", "H3K27me3", "H3K36me3", 
                       "H3K4me1", "ATAC", "H3K79me2", 
                       "H3K9me3", "H3K27ac")
    
    col_fun <- colorRamp2(c(-1, 0, 1), c(col_cell[2], "white", col_cell[1]))
  }
  
  # complexheatmap
  if (T) {
    # norm data
    if (T) {
      # heatmap
      if (T) {
        p <- Heatmap(as.matrix(norm_cor), name = "cor", col=col_fun, 
                     border = "black",
                     cluster_rows = TRUE, cluster_columns = TRUE,
                     right_annotation = rowAnnotation(df = label[, 1:2],
                                                      col = list(cell = col_cell, 
                                                                 mk = col_mk)))
        p <- ggplotify::as.ggplot(p)
      }
      
      file <- "results/2.pic/70.norm_dat_cor.pdf"
      ggsave(filename = file, plot = p, width = 9, height = 8)
    }
    
    # raw data
    if (T) {
      # heatmap
      if (T) {
        p <- Heatmap(as.matrix(raw_cor), name = "cor", col=col_fun, 
                     border = "black",
                     cluster_rows = TRUE, cluster_columns = TRUE,
                     right_annotation = rowAnnotation(df = label[, 1:2],
                                                      col = list(cell = col_cell, 
                                                                 mk = col_mk)))
        p <- ggplotify::as.ggplot(p)
      }
      
      file <- "results/2.pic/70.raw_dat_cor.pdf"
      ggsave(filename = file, plot = p, width = 9, height = 8)
    }
  }
}

# ggplot for depth and enrichment
if (T) {
  # get data 
  if (T) {
    file <- "data/raw_data/3.raw_fq_qc/1.enrichment_depth.txt"
    ggdat <- read.table(file, header = T, sep = "\t", fill = T, comment.char = "")
  }
  
  # format the data
  if (T) {
    ggdat$cell <- str_split(ggdat$Sample, "_", simplify = T)[, 1]
    ggdat$mk <- ifelse(ggdat$Group == "mock", "input", ggdat$Group)
  }
  
  # depth
  if (T) {
    head(ggdat)
    
    p <- plot_average_barplot(ggdat, value="total_passed", average="mean", ylab="Sequence Depth")
    
    file <- "results/2.pic/70.raw_dat_depth.pdf"
    ggsave(filename = file, plot = p, width = 8, height = 6)
  }
  
  # depth
  if (T) {
    head(ggdat)
    
    p <- plot_average_barplot(ggdat, value="Synthetic.JS.Distance", average="mean", ylab="Signal Enrichment")
    
    file <- "results/2.pic/70.raw_dat_signal_enrichment.pdf"
    ggsave(filename = file, plot = p, width = 8, height = 6)
  }
}
