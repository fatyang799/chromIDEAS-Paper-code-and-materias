# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(reshape2)
  library(stringr)
  library(qs)
}

# setting samples
required_sample <- "cd34_|thp1_"
removed_sample <- "thp1_rep[123]"

# get rpkm expression level (tx level)
if (T) {
  # load the tx RPKM data
  if (T) {
    file <- "data/saved_data/4.tx_fpkm_matrix_thp1_cd34.qs"
    rpkm <- qread(file, nthreads = 6)
    file <- "data/saved_data/4.tx_geneid_thp1_cd34.qs"
    gene2id <- qread(file, nthreads = 6)
  }
  
  # get target samples
  if (T) {
    head(rpkm)
    rpkm <- rpkm[, grepl(required_sample, colnames(rpkm)) & !grepl(removed_sample, colnames(rpkm))]
    head(rpkm)
  }
  
  # summary the data
  if (T) {
    sapply(rpkm, summary)
    non0dat <- rpkm[apply(rpkm, 1, sum)>0, ]
  }
  
  # visualization for data
  if (T) {
    file <- "results/2.pic/7.boxplot_expression_level_non0rpkm_tx.jpeg"
    jpeg(file, width = 8, height = 6, units = "in", res = 300)
    boxplot(non0dat, outline=F, las=2, ylab=paste0("non0RPKM value"))
    dev.off()
    
    pheatmap::pheatmap(cor(non0dat),
                       main = "Correlation analysis with non0RPKM value",
                       show_colnames = F, show_rownames = T, display_numbers = T,
                       breaks = seq(min(cor(rpkm))-0.1,1,length.out = 100),
                       border_color = "black", colorRampPalette(c("blue", "white", "red"))(100),
                       filename = "results/2.pic/7.correlation_non0RPKM_tx.jpeg", width = 8, height = 6,
                       fontsize = 13)
  }
  
  # get average tx expression level: rpkm_med
  if (T) {
    colnames(rpkm)
    
    file <- "data/saved_data/7.median_rpkm_value_for_thp1_cd34.qs"
    if (file.exists(file)) {
      rpkm_med <- qread(file, nthreads = 6)
    }
    if (! file.exists(file)) {
      rpkm_med <- lapply(c("thp1", "cd34"), function(cell){
        data <- rpkm[, grepl(cell, colnames(rpkm), ignore.case = T)]
        data <- apply(data, 1, median)
        return(data)
      })
      rpkm_med <- do.call(cbind, rpkm_med)
      colnames(rpkm_med) <- c("thp1", "cd34")
      rpkm_med <- as.data.frame(rpkm_med)
      
      qsave(rpkm_med, file = file, nthreads = 6)
    }
  }
}

# qc for rpkm expression level (tx level)
if (T) {
  # density plot: with 0
  if (T) {
    summary(rpkm_med)
    boxplot(rpkm_med, outline=F, las=2, ylab=paste0("All_RPKM value"))
    
    # get filter data statistics
    if (T) {
      ggdat <- data.frame()
      
      for (i in c(seq(0,0.9,0.1), 1:15)) {
        # i <- 1
        # get filter data
        ggdat[paste0("cutoff_", i), "thp1"] <- sum(rpkm_med[,1]<=i)*100 / nrow(rpkm_med)
        ggdat[paste0("cutoff_", i), "cd34"] <- sum(rpkm_med[,2]<=i)*100 / nrow(rpkm_med)
      }
      
      ggdat$filter <- c(seq(0,0.9,0.1), 1:15)
    }
    
    # ggplot2
    if (T) {
      ggdat <- melt(ggdat, id.vars = "filter")
      colnames(ggdat)
      
      p <- ggplot(data = ggdat)+
        geom_line(aes(x=filter, y=value, group=variable, color=variable), linewidth=1, alpha=0.7)+
        scale_x_continuous(name="rpkm cutoff", breaks = 0:15)+
        scale_y_continuous(name="Cumulative Percentage (%)", breaks = seq(55,100,5)) +
        theme_bw() +
        theme(axis.title = element_text(size = rel(1.2)),
              axis.text = element_text(size = rel(1.2)),
              legend.text = element_text(size = rel(1.2)),
              legend.title = element_text(size = rel(1.2)),
              strip.text = element_text(size = rel(1.2)))
      
      file <- "results/2.pic/7.tx_expression_level_with0_cumulative_curve_plot.jpeg"
      ggsave(file, plot = p, width = 8, height = 6)
    }
  }
  
  # density plot: without 0
  if (T) {
    summary(rpkm_med)
    
    # get filter data statistics
    if (T) {
      ggdat <- data.frame()
      
      for (i in c(seq(0,0.9,0.1), 1:15)) {
        # i <- 1
        # get filter data
        ggdat[paste0("cutoff_", i), "thp1"] <- sum(rpkm_med[,1]<=i & rpkm_med[,1]!=0)*100 / sum(rpkm_med[,1]!=0)
        ggdat[paste0("cutoff_", i), "cd34"] <- sum(rpkm_med[,2]<=i & rpkm_med[,2]!=0)*100 / sum(rpkm_med[,2]!=0)
      }
      
      ggdat$filter <- c(seq(0,0.9,0.1), 1:15)
    }
    
    # ggplot2
    if (T) {
      ggdat <- melt(ggdat, id.vars = "filter")
      colnames(ggdat)
      
      p <- ggplot(data = ggdat)+
        geom_line(aes(x=filter, y=value, group=variable, color=variable), linewidth=1, alpha=0.7)+
        scale_x_continuous(name="rpkm cutoff", breaks = 0:17)+
        scale_y_continuous(name="Cumulative Percentage (%)", breaks = seq(0,100,5)) +
        theme_bw() +
        theme(axis.title = element_text(size = rel(1.2)),
              axis.text = element_text(size = rel(1.2)),
              legend.text = element_text(size = rel(1.2)),
              legend.title = element_text(size = rel(1.2)),
              strip.text = element_text(size = rel(1.2)))
      
      file <- "results/2.pic/7.tx_expression_level_without0_cumulative_curve_plot.jpeg"
      ggsave(file, plot = p, width = 8, height = 6)
    }
  }
  
  rm(ggdat, p, file, i)
}

# classify into 4 class (tx level)
for (cell in colnames(rpkm_med)) {
  # cell <- colnames(rpkm_med)[1]
  cat(paste0("\n", cell, ":\n"))
  
  # get data
  if (T) {
    subdat <- rpkm_med[, cell]
    names(subdat) <- rownames(rpkm_med)
    subdat <- sort(subdat, decreasing = T)
    
    sorted_decreasing_tx <- subdat  
  }
  
  # group data
  if (T) {
    subdat0 <- names(subdat)[subdat == 0]
    
    subdat <- subdat[subdat != 0]
    subdat_low <- names(subdat)[subdat<=quantile(subdat,0.1)]
    subdat_med <- names(subdat)[subdat>quantile(subdat,0.1) & subdat<=quantile(subdat,0.95)]
    subdat_hig <- names(subdat)[subdat>quantile(subdat,0.95)]
  }
  
  # the expression level
  if (T) {
    ggdat <- data.frame(x = 1:length(subdat), 
                        y = log(subdat), 
                        type = ifelse(subdat<=quantile(subdat,0.1), "Low", 
                                      ifelse(subdat>quantile(subdat,0.95), "High", "Med")))
    p <- ggplot(data = ggdat) +
      geom_point(aes(x=x, y=y, color=type)) +
      geom_hline(yintercept=log(quantile(subdat,0.1))) +
      geom_hline(yintercept=log(quantile(subdat,0.95))) +
      ggtitle(cell) +
      scale_x_continuous(name="Genes Expression Level Rank", breaks=NULL) +
      scale_y_continuous(name="log(non0rpkm)", 
                         breaks = c(log(quantile(subdat,0.1)), log(quantile(subdat,0.95))), 
                         labels=c("Q10", "Q95")) +
      cowplot::theme_cowplot() +
      theme(axis.title = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1.2)),
            legend.position = "bottom", 
            legend.text = element_text(size = rel(1.2)),
            legend.title = element_text(size = rel(1.2)),
            strip.text = element_text(size = rel(1.2))) +
      guides(color=guide_legend(nrow=1, title="Expression Level"))
    
    file <- paste0("results/2.pic/7.4_type_tx_expression_level_classification_", cell, ".jpeg")
    ggsave(file, plot = p, width = 8, height = 6)
  }
  
  # save the data
  if (T) {
    dat<- list(sorted_decreasing_tx=sorted_decreasing_tx,
               subdat0=subdat0,
               subdat_low=subdat_low,
               subdat_med=subdat_med,
               subdat_hig=subdat_hig)
    
    # save the classification
    file <- paste0("data/saved_data/7.4_type_tx_expression_level_classification_", cell, ".qs")
    qsave(dat, file = file, nthreads = 6)
  }
  
  # plot the distribution for each group data
  if (T) {
    ggdat <- data.frame(value = rev(sorted_decreasing_tx),
                        group = c(rep("0", length(subdat0)),
                                  rep("Low", length(subdat_low)),
                                  rep("Med", length(subdat_med)),
                                  rep("High", length(subdat_hig))))
    non0min <- min(subdat)
    ggdat$value <- log(ggdat$value+non0min)
    ggdat$group <- factor(ggdat$group, levels = c("0", "Low", "Med", "High"))
    
    # message 
    if (T) {
      mess <- tapply(ggdat$value, ggdat$group, median)
      mess <- paste0(names(mess), " group: ", round(mess, 4))
      mess <- paste(mess, collapse = "\n  ")
      mess <- paste0("The median for each group:\n  ", mess, "\n")
      cat(mess)
    }
    
    # plot 
    if (T) {
      p <- ggplot() + 
        geom_boxplot(data=ggdat, aes(x=group, y=value, fill=group)) +
        geom_text(aes(x=0.5, y=max(ggdat[ggdat$group == "High", 1]), label=mess, hjust="left", vjust="top")) +
        labs(x = "Gruop", y = "log(rpkm)") +
        ggtitle(cell) +
        cowplot::theme_cowplot() +
        theme(axis.title = element_text(size = rel(1.2)),
              axis.text = element_text(size = rel(1.2)),
              legend.position = "none")
    }
    
    # save the figure
    file <- paste0("results/2.pic/7.4_type_tx_expression_level_classification_(", cell, ").jpeg")
    ggsave(file, plot = p, width = 8, height = 6)
  }
}

