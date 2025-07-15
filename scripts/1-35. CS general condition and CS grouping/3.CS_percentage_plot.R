# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(reshape2)
  library(ggplot2)
  library(qs)
}

file <- "data/raw_data/2.states/chromIDEAS.state"

# load the state file
if (T) {
  state <- data.table::fread(file, sep = " ", header = T, data.table = F)
}

# common value
if (T) {
  all_states <- sort(unique(as.numeric(state$cd34)))
}

# format the state data
if (T) {
  head(state)
  state <- state[, c(1, 5, 6)]
  colnames(state)[1] <- "ID"
  head(state)
}

# statistics for rand index
if (T) {
  # all states
  if (T) {
    # 0.720895932311724
    all_ri <- flexclust::randIndex(table(state$cd34, state$thp1),correct = F)
    # 0.427793255610639
    all_ari <- flexclust::randIndex(table(state$cd34, state$thp1),correct = T)
    
    print("Chromatin state similarity of all regions in 2 cells:")
    message <- paste0("RI value: ", all_ri)
    print(message)
    message <- paste0("ARI value: ", all_ari)
    print(message)
  }
  
  # non0state
  if (T) {
    non0state <- state[state$cd34 != 0 & state$thp1 != 0, ]
    
    # 0.931303522587957
    non0_ri <- flexclust::randIndex(table(non0state$cd34, non0state$thp1),correct = F)
    # 0.664525548250809
    non0_ari <- flexclust::randIndex(table(non0state$cd34, non0state$thp1),correct = T)
    
    print("Chromatin state similarity of simultaneous non-state0 regions in 2 cells:")
    message <- paste0("RI value: ", non0_ri)
    print(message)
    message <- paste0("ARI value: ", non0_ari)
    print(message)  
  }
  
  # ggplot
  if (T) {
    ggdat <- data.frame(all_ari=all_ari, 
                        non0_ari=non0_ari)
    ggdat <- melt(ggdat, id.vars = NULL, variable.name = "Type", value.name = "ARI")
    
    p <- ggplot(ggdat, aes(x=Type, y=ARI)) +
      geom_bar(aes(group=Type, fill=Type), stat = "identity", color="black") +
      scale_x_discrete(name="Region", breaks=c("all_ari", "non0_ari"), labels=c("All Genomics", "Non state0 Genomics")) +
      scale_y_continuous(limit=c(0, 1), breaks=seq(0, 1, 0.2)) +
      cowplot::theme_cowplot() +
      theme(legend.position = "None")
    
    file <- "results/2.pic/3.thp1_cd34_CS_similarity_ARI.pdf"
    ggsave(file, plot = p, width = 4, height = 6)
  }
  
  rm(all_ri, all_ari, message, non0state, non0_ri, non0_ari, ggdat, p, file)
}

# get order
if (T) {
  ord <- qread("data/saved_data/2.order_of_emission_table.qs", nthreads = 6)
  ord <- ord[[1]]
}

# statistics for percentage
if (T) {
  percentage <- data.frame(
    sapply(state[, 2:3], function(x) {
      table(x)/nrow(state)*100
    })
  )
  
  percentage$State <- rownames(percentage)
}

# color setting
if (T) {
  display.brewer.all()
  colors <- brewer.pal(9, "Set1") 
  scales::show_col(colors, labels=T)
  
  col_cell <- colors[1:2]
  names(col_cell) <- c("cd34", "thp1")
}

# plot (absolute number)
if (T) {
  # wide2long
  if (T) {
    ggdat <- melt(percentage, id.var="State", value.name = "Percentage", variable.name = "Cell")
    
    # sorting
    head(ggdat)
    ggdat$State <- factor(ggdat$State, levels = all_states)
  }
  
  # ggplot
  if (T) {
    head(ggdat)
    
    p <- ggplot(data=ggdat) +
      geom_bar(aes(x=State, y=Percentage, fill=Cell), position="dodge", stat="identity") +
      scale_fill_manual(values = col_cell) +
      ylab("State_Genomic_Percentage") +
      scale_y_continuous(breaks = c(seq(1,15,2), seq(20,70,5))) +
      ggbreak::scale_y_break(c(18, 50), space=0.1, scales=0.8) + 
      ggbreak::scale_y_break(c(7.5, 15), space=0.1, ticklabels=c(15, 18), scales=0.1) +
      theme_bw() +
      theme(axis.title = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1.2)),
            legend.text = element_text(size = rel(1.2)),
            legend.title = element_text(size = rel(1.2)),
            legend.position = "none",
            strip.text = element_text(size = rel(1.2)))
  }
  
  file <- "results/2.pic/3.thp1_cd34_CS_absolute_percentage.pdf"
  ggsave(file, plot = p, width = 10, height = 8, onefile=F)
  rm(p, ggdat, file)
}

# plot (logFC of percentage)
if (T) {
  # calculating
  if (T) {
    ggdat <- percentage
    ggdat$log2FC <- log2(ggdat$thp1/ggdat$cd34)
  }
  
  # sorting
  if (T) {
    ggdat$State <- factor(ggdat$State, levels = all_states)
    head(ggdat) 
  }
  
  # plot
  if (T) {
    ggdat$main_cell <- ifelse(ggdat$log2FC>0, "thp1", "cd34")
    head(ggdat)
    
    p <- ggplot(data=ggdat) +
      geom_hline(yintercept=c(-1, 1), linetype=2, linewidth=0.8, alpha=0.8, colour=colors[9]) + 
      geom_bar(aes(x=State, y=log2FC, fill=main_cell), stat="identity", color="black") +
      scale_fill_manual(values=col_cell) +
      xlab("State") +
      ylab("log2(THP1/CD34)") +
      theme_bw() +
      theme(axis.title = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1.2)),
            legend.position = "none",
            strip.text = element_text(size = rel(0.8)))
  }
  
  file <- "results/2.pic/3.thp1_cd34_CS_log2fc_percentage.pdf"
  ggsave(file, plot = p, width = 10, height = 3)
}

# plot merged percentage
if (T) {
  # format the data
  if (T) {
    percentage <- c(as.matrix(state[, 2:3]))
    percentage <- table(percentage) / (nrow(state)*2) * 100
    percentage <- as.data.frame.table(percentage)
    colnames(percentage) <- c("State", "Percentage")
    percentage$State <- paste0("S", percentage$State)
    rownames(percentage) <- percentage$State
    
    percentage <- percentage[ord, ]
    percentage$State <- factor(percentage$State, levels = percentage$State)
  }
  
  # ggplot2: barplot
  if (T) {
    head(percentage)
    
    p <- ggplot(data=percentage) +
      geom_bar(aes(x=State, y=Percentage), fill="#3953A4", stat="identity") +
      ylab("State_Genomic_Percentage") +
      scale_y_continuous(breaks = c(seq(1,15,2), seq(20,70,5))) +
      ggbreak::scale_y_break(c(15, 55), scales=0.3) + 
      ggbreak::scale_y_break(c(6, 10), ticklabels=c(10, 15), scales=0.1) +
      theme_bw() +
      theme(axis.title = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1.2)),
            legend.text = element_text(size = rel(1.2)),
            legend.title = element_text(size = rel(1.2)),
            legend.position = "none",
            strip.text = element_text(size = rel(1.2)))
    
    file <- "results/2.pic/3.merged_CS_absolute_percentage_bar.pdf"
    ggsave(file, plot = p, width = 10, height = 8, onefile=F)
    rm(p, file)
  }
  
  # heatmap
  if (T) {
    col_fun = colorRamp2(c(0, ceiling(quantile(c(percentage$Percentage), 0.95))),
                         c("#FFFFFF", "#3953A4"))
    percentage <- as.matrix(percentage[, 2])
    rownames(percentage) <- ord
    colnames(percentage) <- "Percentage"
    
    p <- Heatmap(percentage, name = "Percentage", col = col_fun, 
                 border_gp = gpar(col = "black"), rect_gp = gpar(col = "black"),
                 show_heatmap_legend = TRUE, heatmap_legend_param = list(border = T, title_position = "lefttop", legend_direction = "horizontal"),
                 row_names_max_width = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), row_names_rot = 0,
                 show_column_names = F, 
                 cluster_rows = F, cluster_columns = F, show_row_dend = F, show_column_dend = F)
    
    file <- "results/2.pic/3.merged_CS_absolute_percentage_heatmap.pdf"
    p <- ggplotify::as.ggplot(p)
    ggsave(filename = file, plot = p, width = 2.5, height = 10)
  }
}

