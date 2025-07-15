# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(stringr)
  library(reshape2)
  library(qs)
  suppressPackageStartupMessages(library(ComplexHeatmap))
}

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
up_bin_num <- 5
down_bin_num <- 5
body_bin_num <- 10
rna_levels <- 10

# define functions
if (T) {
  calculate_cell_specific_auc <- function(input, rm_up_down) {
    # test dat
    if (F) {
      input <- "data/saved_data/36.tx_level_profile_Body_5UP_5DW_10Body.bs200.format_genomic_percentage_mRNA_classification_Q10_mat.cd34.qs"
      rm_up_down <- F
    }
    
    # get parameters
    if (T) {
      type <- str_extract(basename(input), "TSS|Body")
      cell <- str_extract(basename(input), paste0(cell1, "|", cell2))
      scale_type <- str_split(str_extract(basename(input), "genomic_percentage|state_percentage"), "_", simplify = T)[, 1]
    }
    
    # get distribution data
    if (T) {
      dat <- qread(input, nthreads = 6)
    }
    
    # format the data
    if (T) {
      colnames(dat) <- c("State", "Loc", "Percentage", "mRNA")
      dat$group <- paste0(dat$mRNA, "@", dat$State)
      mat <- dcast(dat, group~Loc, value.var = "Percentage")
      mat$State <- str_split(mat$group, "@", simplify = T)[, 2]
      mat$mRNA <- str_split(mat$group, "@", simplify = T)[, 1]
      
      if (type == "Body") {
        order <- c(
          paste0("U", up_bin_num:1), 
          "TSS", 
          paste0("B", 1:body_bin_num), 
          "TES", 
          paste0("D", 1:down_bin_num), 
          "State", "mRNA"
        )
      }
      if (type == "TSS") {
        order <- c(
          paste0("U", up_bin_num:1), 
          "TSS", 
          paste0("D", 1:down_bin_num), 
          "State", "mRNA"
        )
      }
      mat <- mat[, order]
    }
    
    # common value
    if (T) {
      state_order <- paste0("S", sort(as.numeric(gsub("S", "", unique(mat$State)))))
      mrna_type <- c("subdat0", paste0("Q", 1:rna_levels))
    }
    
    # remove up and down stream for geneBody region
    pattern <- ifelse(rm_up_down & type == "Body", "TSS|^B|TES", "^U|TSS|^B|TES|^D")
    
    # calculation the AUC value
    if (T) {
      auc <- data.frame(t(sapply(state_order, function(s) {
        # s <- state_order[1]
        auc <- sapply(mrna_type, function(r) {
          # r <- mrna_type[1]
          line_dat <- mat[mat$State == s & mat$mRNA == r, grep(pattern, colnames(mat))]
          line_dat <- as.numeric(unlist(line_dat))
          
          # assume the x gap is 1
          area <- sapply(2:length(line_dat), function(n) {
            # n <- 2
            n2 <- line_dat[n]
            n1 <- line_dat[n-1]
            
            # The formula for calculating the area of a trapezoid
            auc <- (n1+n2) * 1 / 2
            return(auc)
          })
          auc <- sum(area)
          
          return(auc)
        })
        
        names(auc) <- mrna_type
        return(auc)
      })))
      
      auc$State <- state_order
    }
    
    # format the data
    if (F) {
      ggdat <- melt(auc, id.vars = "State", variable.name = "mRNA", value.name = "area")
      ggdat$mRNA <- factor(ggdat$mRNA, levels = mrna_type)
      ggdat$State <- factor(ggdat$State, levels = state_order)
    }
    
    return(auc)
  }
  plot_line_add_merged_line <- function(mat, cluster, state_order, mrna_type) {
    # wide2long
    if (T) {
      mat1 <- mat
      mat1$State <- rownames(mat1)
      mat1 <- melt(mat1, id.vars = "State", variable.name = "mRNA")
      
      mat1$State <- factor(mat1$State, levels = state_order)
      mat1$mRNA <- factor(mat1$mRNA, levels = mrna_type)
    }
    
    # add clustering data
    if (T) {
      mat1$cluster <- cluster$seurat_clusters[match(mat1$State, cluster$orig.ident)]
      mat1$cluster <- factor(mat1$cluster, levels = sort(unique(mat1$cluster)))
    }
    
    # plot the line plot
    if (T) {
      head(mat1)
      
      p <- ggplot(mat1) +
        geom_line(aes(x=mRNA, y=value, group=State, color=State), linewidth=1, alpha=0.3) +
        xlab(NULL) +
        ylab("Relative Preference") +
        facet_grid(.~cluster) +
        cowplot::theme_cowplot() +
        theme(legend.position = "none", 
              panel.border = element_rect(color="black"), 
              strip.background = element_rect(fill=NA, color=NA))
    }
    
    # get merged line
    if (T) {
      merged_line <- tapply(mat1$value, paste0(mat1$cluster, "@", mat1$mRNA), median)
      merged_line <- data.frame(id = names(merged_line), 
                                value = as.numeric(merged_line))
      merged_line$cluster <- str_split(merged_line$id, "@", simplify = T)[, 1]
      merged_line$mRNA <- str_split(merged_line$id, "@", simplify = T)[, 2]
      
      merged_line$cluster <- factor(merged_line$cluster, levels = levels(mat1$cluster))
      merged_line$mRNA <- factor(merged_line$mRNA, levels = levels(mat1$mRNA))
    }
    
    # add merged line into the ggplot
    if (T) {
      head(merged_line)
      merged_p <- p + ggnewscale::new_scale_color() +
        geom_line(data = merged_line, aes(x=mRNA, y=value, group=cluster, color=cluster), linewidth=2)
    }
    
    return(merged_p)
  }
}

# plot the AUC distribution
if (T) {
  data_type <- "tx"
  scale <- "free"
  rms0 <- F
  
  for (type in c("Body")) {
    cat(paste0(type, ":\n"))
    for (scale_type in c("genomic")) {
      cat(paste0("\t", scale_type, ":\n"))
      for (rm_up_down in c(T, F)) {
        name <- ifelse(rm_up_down, "only", "UD")
        cat(paste0("\t\t", name, ":\n"))
        
        # test data
        if (F) {
          type <- "Body"
          scale_type <- "genomic"
          rm_up_down <- F
          name <- ifelse(rm_up_down, "only", "UD")
        }
        
        # mkdir
        if (T) {
          dir <- paste0("results/2.pic/38.stata_distribution_AUC/Body", name)
          if (! dir.exists(dir)) {
            dir.create(dir, showWarnings = F, recursive = T)
          }
        }
        
        # for each single cell: line plot absolute number
        if (T) {
          # calculation
          if (T) {
            ggdat <- lapply(c(cell1, cell2, "merged"), function(cell) {
              # cell <- "merged"
              
              # get input
              if (T) {
                input <- ifelse(type == "TSS", 
                                paste0("data/saved_data/36.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                       "DW.bs", bin_size, ".format_", scale_type, 
                                       "_percentage_mRNA_classification_Q", rna_levels, "_mat.", cell, ".qs"), 
                                paste0("data/saved_data/36.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                       "DW_", body_bin_num, "Body.bs", bin_size, ".format_", scale_type, 
                                       "_percentage_mRNA_classification_Q", rna_levels, "_mat.", cell, ".qs"))
              }
              
              # get auc value data
              if (T) {
                auc <- calculate_cell_specific_auc(input, rm_up_down)
              }
              
              # common value
              if (T) {
                state_order <- paste0("S", sort(as.numeric(gsub("S", "", unique(auc$State)))))
                mrna_type <- c("subdat0", paste0("Q", 1:rna_levels))
              }
              
              # format the data
              if (T) {
                ggdat <- melt(auc, id.vars = "State", variable.name = "mRNA", value.name = "area")
                ggdat$mRNA <- factor(ggdat$mRNA, levels = mrna_type)
                ggdat$State <- factor(ggdat$State, levels = state_order)
                ggdat$Cell <- cell
              }
              
              return(ggdat)
            })
          }
          
          # merge the data
          if (T) {
            ggdat <- do.call(rbind, ggdat)
          }
          
          # get output
          if (T) {
            file <- paste0(dir, "/38.CS_distribution_", type, name, "_merged_", data_type, "_level_", 
                           scale_type, "_AUC_mRNA_Q", rna_levels, ".pdf")
          }
          
          # setting
          if (T) {
            title <- paste0(data_type, "+", scale_type)
            ggdat$group <- paste0(ggdat$Cell, "@", ggdat$State)
            
            # color setting
            if (T) {
              library(RColorBrewer)
              colors <- brewer.pal(9, "Set1") 
              
              col_cell <- c(colors[1:2], "black")
              names(col_cell) <- c("cd34", "thp1", "merged")
            }
          }
          
          # ggplot
          if (T) {
            head(ggdat)
            
            p <- ggplot(ggdat) +
              geom_line(aes(x=mRNA, y=area, group=group, color=Cell), linewidth=1) +
              geom_point(aes(x=mRNA, y=area), size=1, alpha=0.8) +
              ggtitle(title) +
              scale_color_manual(values = col_cell) +
              scale_x_discrete(breaks=NULL) +
              xlab(paste("Expression Level: 0", paste(paste0("Q", 1:10), collapse = "  ->  "), sep = "  ->  ")) +
              ylab("AUC (Area Under Distribution Curve)") +
              ggh4x::facet_grid2(State~Cell, scale = scale, independent = T) +
              cowplot::theme_cowplot() +
              theme(legend.position = c(0.8, 0.1), 
                    panel.border = element_rect(color="black"), 
                    strip.background = element_rect(fill=NA, color=NA))
            
            ggsave(filename = file, plot = p, width = 4, height = 20)
          }
          
          # ggplot for each cell
          if (T) {
            head(ggdat)
            
            plot_dis <- function(cell) {
              p <- ggplot(ggdat[ggdat$Cell == cell, ]) +
                geom_line(aes(x=mRNA, y=area, group=group), linewidth=1) +
                geom_point(aes(x=mRNA, y=area), size=1, alpha=0.8) +
                ggtitle(cell) +
                scale_x_discrete(breaks=NULL) +
                xlab(paste("Expression Level: 0", paste(paste0("Q", 1:10), collapse = "  ->  "), sep = "  ->  ")) +
                ylab("AUC (Area Under Distribution Curve)") +
                facet_wrap(~State, scale = scale, ncol = 5) +
                cowplot::theme_cowplot() +
                theme(legend.position = c(0.8, 0.1), 
                      panel.border = element_rect(color="black"), 
                      strip.background = element_rect(fill=NA, color=NA))
              
              return(p)
            }
            
            p <- plot_dis(cell1)
            ggsave(filename = paste0(dir, "/38.CS_distribution_", type, name, "_merged_", data_type, "_level_", 
                                     scale_type, "_AUC_mRNA_Q", rna_levels, ".", cell1, ".pdf"), 
                   plot = p, width = 10, height = 12)
            
            p <- plot_dis(cell2)
            ggsave(filename = paste0(dir, "/38.CS_distribution_", type, name, "_merged_", data_type, "_level_", 
                                     scale_type, "_AUC_mRNA_Q", rna_levels, ".", cell2, ".pdf"), 
                   plot = p, width = 10, height = 12)
            
            p <- plot_dis("merged")
            ggsave(filename = paste0(dir, "/38.CS_distribution_", type, name, "_merged_", data_type, "_level_", 
                                     scale_type, "_AUC_mRNA_Q", rna_levels, ".", "merged", ".pdf"), 
                   plot = p, width = 10, height = 12)
          }
        }
        
        # for merged situation (paper): line plot absolute number
        if (T) {
          # get input
          if (T) {
            input <- ifelse(type == "TSS", 
                            paste0("data/saved_data/36.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                   "DW.bs", bin_size, ".format_", scale_type, 
                                   "_percentage_mRNA_classification_Q", rna_levels, "_mat.merged.qs"), 
                            paste0("data/saved_data/36.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                   "DW_", body_bin_num, "Body.bs", bin_size, ".format_", scale_type, 
                                   "_percentage_mRNA_classification_Q", rna_levels, "_mat.merged.qs"))
          }
          
          # get auc value data
          if (T) {
            ggdat <- calculate_cell_specific_auc(input, rm_up_down)
          }
          
          # scale the data
          if (T) {
            mat <- ggdat[, colnames(ggdat) != "State"]
            mat <- data.frame(t(scale(t(mat))))
          }
          
          # common value
          if (T) {
            state_order <- paste0("S", sort(as.numeric(gsub("S", "", unique(ggdat$State)))))
            mrna_type <- c("subdat0", paste0("Q", 1:rna_levels))
          }
          
          # load clustering data
          if (T) {
            cluster <- qread("data/saved_data/29.common_cluster_name_for_merged_2cells_CS_clusters.qs", nthreads = 6)
            cluster$seurat_clusters <- paste0("C", cluster$seurat_clusters)
          }
          
          # plot
          if (T) {
            p <- plot_line_add_merged_line(mat, cluster, state_order, mrna_type)
          }
          
          file <- paste0("results/2.pic/38.CS_distribution_", type, name, "_AUC_mRNA_Q", rna_levels, "_", data_type, "_", 
                         scale_type, "_lineplot.pdf")
          ggsave(filename = file, plot = p, width = 10, height = 3)
          
          # kmeans merged figure
          if (T) {
            # load clustering data
            if (T) {
              cluster <- read.table("results/1.tab/40.emission_kmeans_cluster.csv", header = T, sep = ",", fill = T, comment.char = "")
              cluster <- cluster[, c("state", "cluster5")]
              colnames(cluster) <- c("orig.ident", "seurat_clusters")
              cluster$seurat_clusters <- paste0("C", cluster$seurat_clusters)
            }
            
            # plot
            if (T) {
              p <- plot_line_add_merged_line(mat, cluster, state_order, mrna_type)
            }
            
            file <- paste0("results/2.pic/38.CS_distribution_", type, name, "_AUC_mRNA_Q", rna_levels, "_", data_type, "_", 
                           scale_type, "_lineplot.kmeans5.pdf")
            ggsave(filename = file, plot = p, width = 10, height = 3)
          }
        }
        
        # for merged situation (paper): heatmap
        if (T) {
          # get input
          if (T) {
            input <- ifelse(type == "TSS", 
                            paste0("data/saved_data/36.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                   "DW.bs", bin_size, ".format_", scale_type, 
                                   "_percentage_mRNA_classification_Q", rna_levels, "_mat.merged.qs"), 
                            paste0("data/saved_data/36.", data_type, "_level_profile_", type, "_", up_bin_num, "UP_", down_bin_num, 
                                   "DW_", body_bin_num, "Body.bs", bin_size, ".format_", scale_type, 
                                   "_percentage_mRNA_classification_Q", rna_levels, "_mat.merged.qs"))
          }
          
          # get auc value data
          if (T) {
            ggdat <- calculate_cell_specific_auc(input, rm_up_down)
          }
          
          # format the data
          if (T) {
            mat <- as.matrix(ggdat[, colnames(ggdat) != "State"])
            mat <- t(scale(t(mat)))
          }
          
          # color setting
          if (T) {
            library(RColorBrewer)
            colors <- brewer.pal(9, "Set1") 
            
            colors <- colors[1:2]
            range <- sort(c(range(mat), 0))
            col_fun <- colorRamp2(range, c(colors[2], "white", colors[1]))
          }
          
          # load clustering data
          if (T) {
            cluster <- qread("data/saved_data/29.common_cluster_name_for_merged_2cells_CS_clusters.qs", nthreads = 6)
            mat <- mat[cluster$orig.ident, ]
          }
          
          # plot the heatmap
          if (T) {
            p <- Heatmap(mat, name = "Enrichment Level", col = col_fun,
                         
                         border_gp = gpar(col = "black"), rect_gp = gpar(col = "black"),
                         
                         show_heatmap_legend = TRUE, heatmap_legend_param = list(border = T),
                         
                         column_title = "Expression Level", row_title = "C%s",
                         row_title_gp = gpar(fontsize = 13.2), column_title_gp = gpar(fontsize = 13.2),
                         
                         row_split = cluster$seurat_clusters, 
                         
                         row_names_max_width = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), row_names_rot = 0,
                         column_names_max_height = unit(6, "cm"), column_names_gp = gpar(fontsize = 10), column_names_rot = 90,
                         
                         cluster_rows = F, cluster_columns = F)
            p <- ggplotify::as.ggplot(p)
          }
          
          # save the fig
          file <- paste0("results/2.pic/38.CS_distribution_", type, name, "_AUC_mRNA_Q", rna_levels, "_", data_type, "_", 
                         scale_type, "_heatmap.pdf")
          ggsave(filename = file, plot = p, width = 6, height = 6)
          
          # merge line plot
          if (T) {
            # merge data
            if (T) {
              ggdat <- lapply(split(data.frame(mat), cluster$seurat_clusters), function(subdat) {
                # subdat <- mat[ggdat$seurat_clusters == 1, ]
                apply(subdat, 2, mean)
              })
              names(ggdat) <- paste0("C", names(ggdat))
              ggdat <- data.frame(do.call(rbind, ggdat))
            }
            
            # format the data
            if (T) {
              ggdat$cluster <- rownames(ggdat)
              ggdat <- melt(ggdat, id.vars = "cluster", variable.name = "mRNA", value.name = "z_value")
              
              ggdat$mRNA <- factor(ggdat$mRNA, levels = c("subdat0", paste0("Q", 1:rna_levels)))
              ggdat$cluster <- factor(ggdat$cluster, levels = sort(unique(ggdat$cluster)))
            }
            
            # ggplot line
            if (T) {
              head(ggdat)
              
              p <- ggplot(ggdat) +
                geom_line(aes(x=mRNA, y=z_value, group=cluster, color=cluster), linewidth=1) +
                geom_point(aes(x=mRNA, y=z_value), size=2, alpha=0.8) +
                scale_x_discrete(breaks=NULL) +
                xlab(paste("Expression Level: 0", paste(paste0("Q", 1:10), collapse = "  ->  "), sep = "  ->  ")) +
                ylab("AUC (Area Under Distribution Curve)") +
                facet_grid(cluster~., scales = "fixed") +
                cowplot::theme_cowplot() +
                theme(legend.position = "none", 
                      panel.border = element_rect(color="black"), 
                      strip.background = element_rect(fill=NA, color=NA))
            }
            
            # save the fig
            file <- paste0("results/2.pic/38.CS_distribution_", type, name, "_AUC_mRNA_Q", rna_levels, "_", data_type, "_", 
                           scale_type, "_line.pdf")
            ggsave(filename = file, plot = p, width = 3, height = 6)
          }
        }
      }
    }
  }
}

