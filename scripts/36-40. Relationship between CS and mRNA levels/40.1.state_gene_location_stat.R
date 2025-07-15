# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(stringr)
  library(reshape2)
  suppressPackageStartupMessages(library(fmsb))
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(circlize))
  library(qs)
}

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
up_bin_num <- 5
down_bin_num <- 5
body_bin_num <- 10
rna_levels <- 10
n_ct <- 5

# define functions
if (T) {
  # gpat pattern seq stat (modify based on script12)
  gene_body_state_percent_stat <- function(dat) {
    # get total number for each state: nonbody
    if (T) {
      colnames(dat)
      
      nonbody <- dat[, grep("^U|^D|TSS|TES", colnames(dat))]
      
      nonbody <- lapply(all_states, function(s){
        # s <- 0
        nonbody_dat <- sapply(nonbody, function(c_s) {
          sum(s==c_s)
        })
        return(nonbody_dat)
      })
      
      nonbody <- data.frame(do.call(rbind, nonbody))
      nonbody$State <- paste0("S", all_states)
    }
    
    # get total number for each state: body
    if (T) {
      body <- dat[, grep("^B", colnames(dat))]
      
      body <- melt(body, id.vars = NULL, value.name = "value", variable.name = "ID")
      body$Loc <- str_split(body$ID, "_", simplify = T)[, 1]
      body$State <- str_split(body$ID, "_", simplify = T)[, 2]
      body <- dcast(body, State~Loc, value.var = "value", fun.aggregate = sum)
      body <- body[, c(paste0("B", 1:body_bin_num), "State")]
    }
    
    # merge the info 
    if (T) {
      # make sure the order of body and nonbody are identical
      if (T) {
        body <- body[match(paste0("S", all_states), body$State), ]
        nonbody <- nonbody[match(paste0("S", all_states), nonbody$State), ]
      }
      
      # merge the state number data
      if (T) {
        state_profile <- merge(body, nonbody, by="State")
        state_profile <- state_profile[, c(order, "State")]
        state_profile <- state_profile[match(paste0("S", all_states), state_profile$State), ]
        rownames(state_profile) <- 1:nrow(state_profile)
      }
    }
    
    # convert into percentage matrix
    if (T) {
      state_profile[, order] <- sapply(state_profile[, order], function(x) {
        x/sum(x)*100
      })
    }
    
    return(state_profile)
  }
  
  # genomic percentage across gene body for states
  state_percentage_across_body <- function(ggdat, title) {
    # test data
    if (F) {
      title <- "title"
    }
    
    # ggplot2
    if (T) {
      p <- ggplot(ggdat, aes(x=Loc, y=Percentage)) +
        geom_line(aes(group=State), linewidth=1) +
        scale_x_discrete(name = NULL, breaks=labels) +
        scale_y_continuous(name = "Genomic Percentage") +
        ggtitle(title) +
        facet_wrap(~State, ncol = 5, scales = "free") +
        cowplot::theme_cowplot() +
        theme(panel.border = element_rect(color="black"), 
              strip.background = element_rect(color=NA, fill=NA))
    }
    
    return(p)
  }
  
  # script53
  if (T) {
    # correlation (modify based on script39)
    plot_state_auc_cor_heatmap <- function(cor_dat, cluster) {
      # get cell specific cluster
      if (T) {
        group <- cluster[match(rownames(cor_dat), cluster$orig.ident), "seurat_clusters"]
        group <- paste0("C", group)
      }
      
      # color setting
      if (T) {
        library(RColorBrewer)
        colors <- brewer.pal(9, "Set1") 
        
        colors <- colors[1:2]
        col_fun <- circlize::colorRamp2(c(-1, 0, 1), c(colors[2], "white", colors[1]))
      }
      
      # heatmap
      if (T) {
        p <- Heatmap(cor_dat, name = "Pearson", col = col_fun, 
                     border_gp = gpar(col = "black"), rect_gp = gpar(col = "black"),
                     show_heatmap_legend = TRUE, heatmap_legend_param = list(border = T),
                     row_names_max_width = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), row_names_rot = 0,
                     column_names_max_height = unit(6, "cm"), column_names_gp = gpar(fontsize = 10), column_names_rot = 90,
                     cluster_rows = T, cluster_columns = T, show_row_dend = T, show_column_dend = F, 
                     clustering_distance_columns = "euclidean", clustering_distance_rows = "euclidean", 
                     clustering_method_columns = "ward.D2", clustering_method_rows = "ward.D2", 
                     right_annotation = rowAnnotation(Group = group, border = T, 
                                                      col = list(Group = c("C1" = "#e76f51", 
                                                                           "C2" = "#8ecae6", 
                                                                           "C3" = "#f4a261", 
                                                                           "C4" = "#219ebc", 
                                                                           "C5" = "#a7c957")), 
                                                      annotation_legend_param = list(border = "black")))
      }
      
      p <- ggplotify::as.ggplot(p)
      return(p)
    }
    
    # k-means (script39)
    kmeans_determine_k <- function(mat, title) {
      # test dat
      if (F) {
        title <- "test"
      }
      
      # scale data
      if (T) {
        df <- mat
      }
      
      # determine the k value: Within groups sum of squares
      if (T) {
        wssplot <- function(df, nc=15, title, seed=799){
          wss <- c()
          for (i in 1:nc){
            set.seed(seed)
            wss[i] <- sum(kmeans(df, centers=i, nstart=25)$withinss)
          }
          ggdat <- data.frame(k=1:nc, y=wss)
          
          p <- ggplot(ggdat) +
            geom_point(aes(x=k, y=wss), size=3) +
            geom_line(aes(x=k, y=wss, group=1)) +
            scale_x_continuous(breaks = 1:nc) +
            xlab("Number of Clusters") +
            ylab("Within groups sum of squares") +
            ggtitle(title) +
            theme_bw() +
            theme(axis.title = element_text(size = rel(1.2)),
                  axis.text = element_text(size = rel(1.2)),
                  legend.text = element_text(size = rel(1.2)),
                  legend.title = element_text(size = rel(1.2)),
                  strip.text = element_text(size = rel(1.2)))
          return(p)
        }
        p <- wssplot(df, nc=15, title, seed=799)
      }
      
      return(p)
    }
    
    # k-means (script50)
    cluster_kmeans <- function(mat, k, scale) {
      # test dat
      if (F) {
        k <- 4
        scale <- "row"
      }
      
      df <- mat
      
      # scale data
      if (T) {
        if (scale == "row") {
          df <- data.frame(t(scale(t(df))))
          torf <- apply(df, 1, function(x) {
            sum(is.na(x)) == 0
          })
          df <- df[torf, ]
        }
        if (scale == "col") {
          df <- data.frame(scale(df))
          torf <- apply(df, 1, function(x) {
            sum(is.na(x)) == 0
          })
          df <- df[torf, ]
        }
        if (scale == "none") {
          df <- df
        }
      }
      
      # k-means
      if (T) {
        set.seed(799)
        fit.km <- kmeans(df, centers=k, nstart=25, iter.max=50, algorithm="Lloyd")
      }
      
      return(fit.km)
    }
    
    # line plot clustered by chromIDEAS (script38)
    plot_line_add_merged_line <- function(mat, cluster) {
      # wide2long
      if (T) {
        mat1 <- mat_row
        mat1$State <- rownames(mat1)
        mat1 <- melt(mat1, id.vars = "State", variable.name = "Loc", measure.vars = order, value.name = "percentage")
        
        mat1$State <- factor(mat1$State, levels = paste0("S", all_states))
        mat1$Loc <- factor(mat1$Loc, levels = order)
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
          geom_line(aes(x=Loc, y=percentage, group=State, color=State), linewidth=1, alpha=0.3) +
          scale_x_discrete(name = NULL, breaks = labels) +
          ylab("Relative Preference") +
          facet_grid(.~cluster) +
          cowplot::theme_cowplot() +
          theme(legend.position = "none", 
                panel.border = element_rect(color="black"), 
                strip.background = element_rect(fill=NA, color=NA))
      }
      
      # get merged line
      if (T) {
        merged_line <- tapply(mat1$percentage, paste0(mat1$cluster, "@", mat1$Loc), median)
        merged_line <- data.frame(id = names(merged_line), 
                                  value = as.numeric(merged_line))
        merged_line$cluster <- str_split(merged_line$id, "@", simplify = T)[, 1]
        merged_line$Loc <- str_split(merged_line$id, "@", simplify = T)[, 2]
        
        merged_line$cluster <- factor(merged_line$cluster, levels = levels(mat1$cluster))
        merged_line$Loc <- factor(merged_line$Loc, levels = levels(mat1$Loc))
      }
      
      # add merged line into the ggplot
      if (T) {
        head(merged_line)
        merged_p <- p + ggnewscale::new_scale_color() +
          geom_line(data = merged_line, aes(x=Loc, y=value, group=cluster, color=cluster), linewidth=2)
      }
      
      return(merged_p)
    }
    
    # kmeans plot (modify based on script39)
    plot_state_p_across_body_distribution_kmeans <- function(mat, k, title) {
      # test dat
      if (F) {
        title <- "test"
        k <- 3
      }
      
      # scale data
      if (T) {
        df <- mat
      }
      
      # k-means
      if (T) {
        cluster_kmeans <- function(df, k) {
          set.seed(799)
          fit.km <- kmeans(df, centers=k, nstart=25)
          cluster <- fit.km$cluster
          cluster <- data.frame(state = names(cluster), cluster = as.numeric(cluster))
          return(cluster)
        }
        
        clusterK <- cluster_kmeans(df, k)
      }
      
      # add cluster info into ggdat
      if (T) {
        df$State <- rownames(df)
        ggdat <- melt(df, id.vars = "State", value.name = "percentage", variable.name = "Loc")
        ggdat$cluster <- paste0("C", clusterK[match(ggdat$State, clusterK$state), 2])
      }
      
      # format the data
      if (T) {
        ggdat$State <- factor(ggdat$State, levels = paste0("S", all_states))
        ggdat$cluster <- factor(ggdat$cluster, levels = paste0("C", 1:k))
        ggdat$Loc <- factor(ggdat$Loc, levels = order)
      }
      
      # ggplot
      if (T) {
        head(ggdat)
        
        p <- ggplot(ggdat) +
          geom_line(aes(x=Loc, y=percentage, group=State, color=State), linewidth=1, alpha=0.3) +
          ggtitle(title) +
          scale_x_discrete(breaks=NULL) +
          xlab(paste("Expression Level: 0", paste(paste0("Q", 1:10), collapse = "  ->  "), sep = "  ->  ")) +
          ylab("AUC (Area Under Distribution Curve)") +
          facet_grid(.~cluster, scale = "fixed") +
          cowplot::theme_cowplot() +
          theme(legend.position = "none", 
                panel.border = element_rect(color="black"), 
                strip.background = element_rect(fill=NA, color=NA))
      }
      
      # get merged line
      if (T) {
        merged_line <- tapply(ggdat$percentage, paste(ggdat$cluster, ggdat$Loc), median)
        merged_line <- data.frame(id = names(merged_line), 
                                  value = as.numeric(merged_line))
        merged_line$cluster <- str_split(merged_line$id, " ", simplify = T)[, 1]
        merged_line$Loc <- str_split(merged_line$id, " ", simplify = T)[, 2]
        
        merged_line$cluster <- factor(merged_line$cluster, levels = levels(ggdat$cluster))
        merged_line$Loc <- factor(merged_line$Loc, levels = levels(ggdat$Loc))
      }
      
      # add merged line into the ggplot
      if (T) {
        head(merged_line)
        merged_p <- p + ggnewscale::new_scale_color() +
          geom_line(data = merged_line, aes(x=Loc, y=value, group=cluster, color=cluster), linewidth=2)
      }
      
      # get state label
      if (T) {
        state_cluster <- tapply(clusterK$state, paste0("C", clusterK$cluster), function(x) {
          # x <- clusterK$state[clusterK$cluster==1]
          x <- gsub("S", "", x)
          if (length(x) > 20) {
            x1 <- paste0(paste(x[1:10], collapse = ", "), "\n")
            x2 <- paste0(paste(x[11:20], collapse = ", "), "\n")
            x_remin <- paste(x[21:length(x)], collapse = ", ")
            x <- paste0(x1, x2, x_remin)
          } else if (length(x) > 10) {
            x1 <- paste0(paste(x[1:10], collapse = ", "), "\n")
            x_remin <- paste(x[11:length(x)], collapse = ", ")
            x <- paste0(x1, x_remin)
          } else {
            x <- paste(x, collapse = ", ")
          }
          paste0("S", x)
        })
        state_cluster <- data.frame(
          cluster = names(state_cluster), 
          label = state_cluster
        )
      }
      
      # add merged line into the ggplot
      if (T) {
        head(state_cluster)
        label_p <- merged_p + 
          geom_text(data = state_cluster, aes(x=1, y=max(ggdat$percentage), label=label), hjust = 0, vjust = 1)
      }
      
      out <- list(label_p, clusterK)
      return(out)
    }
  }
  
  # script54
  if (T) {
    # calculate the similarity between 2 cluster result (script40)
    similarity_between_clusters <- function(cluster1, cluster2) {
      # calculate the ari and ri
      similarity <- data.frame(t(sapply(cluster2, function(c) {
        ari <- flexclust::randIndex(table(c, cluster1),correct = T)
        ri <- flexclust::randIndex(table(c, cluster1),correct = F)
        res <- c(ari, ri)
        names(res) <- c("ari", "ri")
        
        return(res)
      })))
      
      return(similarity)
    }
    plot_radar_figure_with_similarity <- function(similarity, title, color, order=NA) {
      # test dat
      if (F) {
        similarity <- rnorm(8)
        names(similarity) <- c(paste0("cluster", 4:10), "chromIDEAS")
        title <- "Consistency of clustering: ARI"
        order <- c(paste0("cluster", 7:10), "chromIDEAS", paste0("cluster", 4:6))
      }
      
      # format the data
      if (T) {
        ggdat <- data.frame(rbind(
          matrix(rep(max(similarity)*1.1, length(similarity)), ncol = length(similarity)), 
          matrix(rep(min(similarity)*0.9, length(similarity)), ncol = length(similarity)), 
          similarity
        ))
        rownames(ggdat) <- c("max", "min", "value")
        colnames(ggdat) <- names(similarity)
        
        if (sum(is.na(order)) == 0 & length(order) == length(similarity)) {
          ggdat <- ggdat[, order]
        }
      }
      
      # radarchart
      if (T) {
        op <- par(mar = c(1, 1, 1, 1))
        
        radarchart(ggdat, 
                   axistype = 1, # 0:5分别表示6类坐标轴方案
                   axislabcol = "black", # 坐标轴标签颜色
                   caxislabels = paste0(seq(0, 100, 25), "%"), # axistype=1时，自定义坐标轴标签，数量=seg+1
                   calcex = 1, # axistype=1时，坐标轴标签字体大小
                   paxislabels = NULL, # axistype=2时，自定义坐标轴标签，数量=ncol(df)
                   palcex = 1, # axistype=2时，坐标轴标签字体大小
                   
                   pty = 16, # 数据点形状
                   pcol = color, # 数据点之间连线颜色
                   plty = 1, # 数据点之间连线类型，1实线
                   plwd = 2, # 数据点之间连线宽度，默认为1
                   pdensity = NULL, # 数据点围成图形中，以斜线绘制密度，null表示不绘制斜线（不能用0，否则pfcol参数失效），1表示绘制1条斜线
                   pangle = 45, # 数据点围成图形中，以斜线绘制密度，设置斜线的角度
                   pfcol = scales::alpha(color, 0.5), # 数据点围成图形中，填充颜色
                   
                   seg = 4, # 网络圈数
                   cglty = 1, # 雷达背景线条类型，2表示虚线
                   cglwd = 0.8, # 雷达背景线条宽度
                   cglcol = "grey", # 雷达背景线条颜色
                   
                   title = title, 
                   na.itp = F, # 对于NA是否自动插值，默认F
                   centerzero = , # 是否scale，默认F
                   
                   vlabels = NULL, # 自定义每个特征的名字，默认为colnames(df)
                   vlcex = 1 # 自定义每个特征名字的字体大小
        )
        par(op)
      }
    }
  }
}

# data prepare
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir1 <- paste0("data/saved_data/40.1.state_gene_location/", data_type)
      if (! dir.exists(dir1)) {
        dir.create(dir1, showWarnings = F, recursive = T)
      }
      
      dir2 <- paste0("results/2.pic/40.1.state_gene_location/", data_type)
      if (! dir.exists(dir2)) {
        dir.create(dir2, showWarnings = F, recursive = T)
      }
    }
    
    for (type in c("merged", "single")) {
      cat(paste0("\t", type, ": \n"))
      
      # cell specific data
      if (type == "single") {
        for (cell in c(cell1, cell2)) {
          cat(paste0("\t\t", cell, ": \n"))
          
          # test data
          if (F) {
            data_type <- "tx"
            type <- "single"
            cell <- cell1
          }
          
          # load data
          if (T) {
            cell_dat <- function(cell) {
              input <- paste0("data/saved_data/12.", data_type, "_level_profile_Body_", up_bin_num, "UP_", 
                              down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".raw_mat.", cell, ".qs")
              
              dat <- qread(input, nthreads = 6)
              dat$gene_id <- paste0(dat$gene_id, "@", cell)
              
              return(dat)
            }
            
            dat <- cell_dat(cell)
          }
          
          # common value
          if (T) {
            order <- c(
              paste0("U", up_bin_num:1), 
              "TSS", 
              paste0("B", 1:body_bin_num), 
              "TES", 
              paste0("D", 1:down_bin_num)
            )
            
            all_states <- sort(as.numeric(unique(dat$TSS)))
          }
          
          # get profile
          if (T) {
            file <- paste0(dir1, "/40.1.state_percentage_in_", data_type, "s_", type, "_", cell, ".qs")
            
            if (! file.exists(file)) {
              mat <- gene_body_state_percent_stat(dat)
              
              qsave(mat, file, nthreads = 6)
            }
          }
        }
      }
      
      # merged data
      if (type == "merged") {
        cat(paste0("\t\tmerged: \n"))
        
        # test data
        if (F) {
          data_type <- "tx"
          type <- "merged"
        }
        
        # load data
        if (T) {
          cell_dat <- function(cell) {
            input <- paste0("data/saved_data/12.", data_type, "_level_profile_Body_", up_bin_num, "UP_", 
                            down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".raw_mat.", cell, ".qs")
            
            dat <- qread(input, nthreads = 6)
            dat$gene_id <- paste0(dat$gene_id, "@", cell)
            
            return(dat)
          }
          dat1 <- cell_dat(cell1)
          dat2 <- cell_dat(cell2)
          
          dat <- rbind(dat1, dat2)
          
          rm(dat1, dat2)
        }
        
        # common value
        if (T) {
          order <- c(
            paste0("U", up_bin_num:1), 
            "TSS", 
            paste0("B", 1:body_bin_num), 
            "TES", 
            paste0("D", 1:down_bin_num)
          )
          
          all_states <- sort(as.numeric(unique(dat$TSS)))
        }
        
        # get profile
        if (T) {
          file <- paste0(dir1, "/40.1.state_percentage_in_", data_type, "s_", type, "_merged.qs")
          
          if (! file.exists(file)) {
            mat <- gene_body_state_percent_stat(dat)
            
            qsave(mat, file, nthreads = 6)
          }
        }
      }
    }
  }
}

# plot state percentage across gene body (modify based on script52+53)
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir1 <- paste0("data/saved_data/40.1.state_gene_location/", data_type)
      if (! dir.exists(dir1)) {
        dir.create(dir1, showWarnings = F, recursive = T)
      }
      
      dir2 <- paste0("results/2.pic/40.1.state_gene_location/", data_type)
      if (! dir.exists(dir2)) {
        dir.create(dir2, showWarnings = F, recursive = T)
      }
    }
    
    for (type in c("merged", "single")) {
      cat(paste0("\t", type, ": \n"))
      
      # cell specific data
      if (type == "single") {
        for (cell in c(cell1, cell2)) {
          cat(paste0("\t\t", cell, ": \n"))
          
          # test data
          if (F) {
            data_type <- "tx"
            type <- "single"
            cell <- cell1
          }
          
          # load the data
          if (T) {
            file <- paste0(dir1, "/40.1.state_percentage_in_", data_type, "s_", type, "_", cell, ".qs")
            mat <- qread(file)
          }
          
          # common value
          if (T) {
            order <- c(
              paste0("U", up_bin_num:1), 
              "TSS", 
              paste0("B", 1:body_bin_num), 
              "TES", 
              paste0("D", 1:down_bin_num)
            )
            labels <- c(paste0("U", up_bin_num), "TSS", "TES", paste0("D", down_bin_num))
            
            all_states <- sort(unique(as.numeric(gsub("S", "", mat$State))))
          }
          
          # format the data to plot state percentage across gene body
          if (T) {
            ggdat <- melt(mat, id.vars = c("State"), measure.vars = order, variable.name = "Loc", value.name = "Percentage")
            
            head(ggdat)
            
            ggdat$State <- factor(ggdat$State, levels = paste0("S", all_states))
            ggdat$Loc <- factor(ggdat$Loc, levels = order)
          }
          
          # state percentage change across gene body
          if (T) {
            head(ggdat)
            
            title <- paste0(cell, ": ", type)
            p <- state_percentage_across_body(ggdat, title)
            
            file <- paste0(dir2, "/40.1.CS_percentage_across_Body_in_", data_type, "s.", cell, ".", type, ".pdf")
            ggsave(filename = file, plot = p, width = 10, height = 12)
          }
          
          # read the cluster data
          if (T) {
            file <- "data/saved_data/25.common_cluster_name_for_2_cell_CS_clusters.qs"
            cluster <- qread(file, nthreads = 6)
            cluster$orig.ident <- paste0("S", cluster$orig.ident)
            cluster <- cluster[cluster$cell == cell, ]
          }
          
          # cor heatmap
          if (T) {
            ggdat <- mat
            rownames(ggdat) <- ggdat$State
            ggdat <- ggdat[, order]
            ggdat <- t(ggdat)
            
            cor_dat <- cor(ggdat, method = "pearson")
            
            p <- plot_state_auc_cor_heatmap(cor_dat, cluster)
            
            file <- paste0(dir2, "/40.1.CS_gene_body_percentage_heatmap_", data_type, "s.", cell, ".", type, ".pdf")
            ggsave(filename = file, plot = p, width = 15, height = 10)
          }
          
          # plot merged line clustered by chromIDEAS
          if (T) {
            mat_row <- mat[, order]
            rownames(mat_row) <- mat$State
            mat_row <- data.frame(t(scale(t(mat_row))))
            
            p <- plot_line_add_merged_line(mat_row, cluster)
            
            file <- paste0(dir2, "/40.1.CS_gene_body_percentage_clustered_by_chromIDEAS_", data_type, "s.", cell, ".", type, ".pdf")
            ggsave(filename = file, plot = p, width = 10, height = 3)
          }
          
          # plot merged line clustered by kmeans5
          if (T) {
            # read the cluster data
            if (T) {
              cluster <- read.table("results/1.tab/40.emission_kmeans_cluster.csv", header = T, sep = ",", fill = T, comment.char = "")
              cluster <- cluster[, c("state", "cluster5")]
              colnames(cluster) <- c("orig.ident", "seurat_clusters")
            }
            
            p <- plot_line_add_merged_line(mat_row, cluster)
            
            file <- paste0(dir2, "/40.1.CS_gene_body_percentage_clustered_by_kmeans5cluster_", data_type, "s.", cell, ".", type, ".pdf")
            ggsave(filename = file, plot = p, width = 10, height = 3)
          }
          
          # kmeans
          if (T) {
            # determine the k value
            if (T) {
              p <- kmeans_determine_k(mat_row, title=paste0("State Percentage across Gene Body kmeans K"))
              
              file <- paste0(dir2, "/40.1.CS_gene_body_percentage_determine_K_", data_type, "s.", cell, ".", type, ".pdf")
              ggsave(filename = file, plot = p, width = 8, height = 6)
            }
            
            k_auc <- 3:6
            
            # kmeans distribution
            if (T) {
              kmeans_c <- lapply(k_auc, function(k) {
                # k <- 5
                dat <- plot_state_p_across_body_distribution_kmeans(mat_row, k, title=paste0("kmeans: k=", k))
                
                # get figure
                file <- paste0(dir2, "/40.1.CS_gene_body_percentage_clustered_by_kmeans", k, "_", data_type, "s.", cell, ".", type, ".pdf")
                ggsave(filename = file, plot = dat[[1]], width = 10, height = 3)
                
                # get cluster
                c <- dat[[2]]
                colnames(c)[2] <- paste0("K", k)
                
                return(c)
              })
              kmeans_c <- do.call(cbind, kmeans_c)
              
              apply(kmeans_c[, seq(1, ncol(kmeans_c), 2)], 1, table)
              
              kmeans_c <- kmeans_c[, c(1, seq(2, ncol(kmeans_c), 2))]
              
              write.table(kmeans_c, file = paste0("results/1.tab/40.1.CS_gene_body_percentage_clustered_by_kmeans_", data_type, "s.", cell, ".", type, ".csv"),
                          quote = F, sep = ",", col.names = T, row.names = F)
            }
          }
        }
      }
      
      # merged data
      if (type == "merged") {
        cat(paste0("\t\tmerged: \n"))
        
        # test data
        if (F) {
          data_type <- "tx"
          type <- "merged"
        }
        
        # load the auc value
        if (T) {
          file <- paste0(dir1, "/40.1.state_percentage_in_", data_type, "s_", type, "_merged.qs")
          mat <- qread(file)
        }
        
        # common value
        if (T) {
          order <- c(
            paste0("U", up_bin_num:1), 
            "TSS", 
            paste0("B", 1:body_bin_num), 
            "TES", 
            paste0("D", 1:down_bin_num)
          )
          labels <- c(paste0("U", up_bin_num), "TSS", "TES", paste0("D", down_bin_num))
          
          all_states <- sort(unique(as.numeric(gsub("S", "", mat$State))))
        }
        
        # format the data to plot state percentage across gene body
        if (T) {
          ggdat <- melt(mat, id.vars = c("State"), measure.vars = order, variable.name = "Loc", value.name = "Percentage")
          
          head(ggdat)
          
          ggdat$State <- factor(ggdat$State, levels = paste0("S", all_states))
          ggdat$Loc <- factor(ggdat$Loc, levels = order)
        }
        
        # state percentage change across gene body
        if (T) {
          head(ggdat)
          
          title <- paste0("merged: ", type)
          p <- state_percentage_across_body(ggdat, title)
          
          file <- paste0(dir2, "/40.1.CS_percentage_across_Body_in_", data_type, "s.merged.", type, ".pdf")
          ggsave(filename = file, plot = p, width = 10, height = 12)
        }
        
        # read the cluster data
        if (T) {
          file <- "data/saved_data/29.common_cluster_name_for_merged_2cells_CS_clusters.qs"
          cluster <- qread(file, nthreads = 6)
        }
        
        # cor heatmap
        if (T) {
          ggdat <- mat
          rownames(ggdat) <- ggdat$State
          ggdat <- ggdat[, order]
          ggdat <- t(ggdat)
          
          cor_dat <- cor(ggdat, method = "pearson")
          
          p <- plot_state_auc_cor_heatmap(cor_dat, cluster)
          
          file <- paste0(dir2, "/40.1.CS_gene_body_percentage_heatmap_", data_type, "s.merged.", type, ".pdf")
          ggsave(filename = file, plot = p, width = 15, height = 10)
        }
        
        # plot merged line clustered by chromIDEAS
        if (T) {
          mat_row <- mat[, order]
          rownames(mat_row) <- mat$State
          mat_row <- data.frame(t(scale(t(mat_row))))
          
          p <- plot_line_add_merged_line(mat_row, cluster)
          
          file <- paste0(dir2, "/40.1.CS_gene_body_percentage_clustered_by_chromIDEAS_", data_type, "s.merged.", type, ".pdf")
          ggsave(filename = file, plot = p, width = 10, height = 3)
        }
        
        # plot merged line clustered by kmeans5
        if (T) {
          # read the cluster data
          if (T) {
            cluster <- read.table("results/1.tab/40.emission_kmeans_cluster.csv", header = T, sep = ",", fill = T, comment.char = "")
            cluster <- cluster[, c("state", "cluster5")]
            colnames(cluster) <- c("orig.ident", "seurat_clusters")
          }
          
          p <- plot_line_add_merged_line(mat_row, cluster)
          
          file <- paste0(dir2, "/40.1.CS_gene_body_percentage_clustered_by_kmeans5cluster_", data_type, "s.merged.", type, ".pdf")
          ggsave(filename = file, plot = p, width = 10, height = 3)
        }
        
        # kmeans
        if (T) {
          # determine the k value
          if (T) {
            p <- kmeans_determine_k(mat_row, title=paste0("State Percentage across Gene Body kmeans K"))
            
            file <- paste0(dir2, "/40.1.CS_gene_body_percentage_determine_K_", data_type, "s.merged.", type, ".pdf")
            ggsave(filename = file, plot = p, width = 8, height = 6)
          }
          
          k_auc <- 3:6
          
          # kmeans distribution
          if (T) {
            kmeans_c <- lapply(k_auc, function(k) {
              # k <- 5
              dat <- plot_state_p_across_body_distribution_kmeans(mat_row, k, title=paste0("kmeans: k=", k))
              
              # get figure
              file <- paste0(dir2, "/40.1.CS_gene_body_percentage_clustered_by_kmeans", k, "_", data_type, "s.merged.", type, ".pdf")
              ggsave(filename = file, plot = dat[[1]], width = 10, height = 3)
              
              # get cluster
              c <- dat[[2]]
              colnames(c)[2] <- paste0("K", k)
              
              return(c)
            })
            kmeans_c <- do.call(cbind, kmeans_c)
            
            apply(kmeans_c[, seq(1, ncol(kmeans_c), 2)], 1, table)
            
            kmeans_c <- kmeans_c[, c(1, seq(2, ncol(kmeans_c), 2))]
            
            write.table(kmeans_c, file = paste0("results/1.tab/40.1.CS_gene_body_percentage_clustered_by_kmeans_", data_type, "s.merged.", type, ".csv"),
                        quote = F, sep = ",", col.names = T, row.names = F)
          }
        }
      }
    }
  }
}

# adjust rand index radar to assess the superiority of various cluster results
if (T) {
  k_auc <- 3:6
  
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir1 <- paste0("data/saved_data/40.1.state_gene_location/", data_type)
      if (! dir.exists(dir1)) {
        dir.create(dir1, showWarnings = F, recursive = T)
      }
      
      dir2 <- paste0("results/2.pic/40.1.state_gene_location/", data_type)
      if (! dir.exists(dir2)) {
        dir.create(dir2, showWarnings = F, recursive = T)
      }
    }
    
    for (type in c("merged", "single")) {
      cat(paste0("\t", type, ": \n"))
      
      # cell specific data
      if (type == "single") {
        for (cell in c(cell1, cell2)) {
          cat(paste0("\t\t", cell, ": \n"))
          
          # test data
          if (F) {
            data_type <- "tx"
            type <- "single"
            cell <- cell1
          }
          
          # get state kmeans cluster as the golden standard
          if (T) {
            file <- paste0("results/1.tab/40.1.CS_gene_body_percentage_clustered_by_kmeans_", data_type, "s.", cell, ".", type, ".csv")
            cluster_golden <- read.table(file, header = T, sep = ",", fill = T, comment.char = "")
          }
          
          # using kmeans to cluster emission table
          if (T) {
            cluster_emission_kmeans <- read.table("results/1.tab/40.emission_kmeans_cluster.csv", header = T, sep = ",")
          }
          
          # read the chromIDEAS cluster
          if (T) {
            cluster_chrom <- qread("data/saved_data/25.common_cluster_name_for_2_cell_CS_clusters.qs", nthreads = 6)
            cluster_chrom$orig.ident <- paste0("S", cluster_chrom$orig.ident)
            cluster_chrom <- cluster_chrom[cluster_chrom$cell == cell, ]
          }
          
          # get random group as control
          if (T) {
            file <- paste0("data/saved_data/40.random_cluster_control_n", n_ct, ".", cell, ".qs")
            cluster_random <- qread(file, nthreads = 6)
          }
          
          # make the order of all clusters are identical
          if (T) {
            cluster_emission_kmeans <- cluster_emission_kmeans[match(cluster_golden$state, cluster_emission_kmeans$state), ]
            cluster_chrom <- cluster_chrom[match(cluster_golden$state, cluster_chrom$orig.ident), ]
            cluster_random <- cluster_random[match(cluster_golden$state, cluster_random$state), ]
            
            cluster_compare <- cbind(cluster_emission_kmeans, cluster_chrom, cluster_random)
            cluster_compare <- cluster_compare[, c("cluster5", paste0("ct", 1:n_ct), "seurat_clusters")]
            colnames(cluster_compare) <- c("kmeans5", paste0("CT", 1:n_ct), "chromIDEAS")
            
            rm(cluster_emission_kmeans, cluster_chrom, cluster_random)
          }
          
          # compare the similarity compared with auc distribution line cluster
          if (T) {
            # setting
            if (T) {
              file <- paste0(dir2, "/40.1.CS_gene_body_percentage_ARI_in_", data_type, "s.", cell, ".", type, ".pdf")
              pdf(file = file, width = 8, height = 6)
              colors <- c("#00AFBB", "#E7B800", "#FC4E07", "#a7c957")
              par(mfrow = c(2,2))
            }
            
            # statistics ari
            if (T) {
              ari_stat <- lapply(k_auc, function(kk) {
                similarity <- similarity_between_clusters(cluster_golden[, paste0("K", kk)], cluster_compare)
                ari <- similarity$ari
                names(ari) <- rownames(similarity)
                
                return(ari)
              })
              names(ari_stat) <- paste0("K", k_auc)
            }
            
            # plot
            if (T) {
              for (kk in k_auc) {
                # setting
                title <- paste0("State % in Gene Body cluster", kk)
                
                plot_radar_figure_with_similarity(ari_stat[[paste0("K", kk)]], 
                                                  title, 
                                                  colors[max(k_auc)+1-kk], 
                                                  c(paste0("CT", c(3:5)), "chromIDEAS", "kmeans5", paste0("CT", c(1:2))))
              }
            }
            
            # save the results
            if (T) {
              ari <- data.frame(t(do.call(rbind, ari_stat)))
              ari$method <- rownames(ari)
              rownames(ari) <- 1:nrow(ari)
              
              file <- paste0("results/1.tab/40.1.CS_gene_body_percentage_ARI_in_", data_type, "s.", cell, ".", type, ".csv")
              write.table(ari, file = file, quote = F, sep = ",", col.names = T, row.names = F)
            }
            
            Sys.sleep(0.05)
            dev.off()
          }
        }
      }
      
      # merged data
      if (type == "merged") {
        cat(paste0("\t\tmerged: \n"))
        
        # test data
        if (F) {
          data_type <- "tx"
          type <- "merged"
        }
        
        # get state kmeans cluster as the golden standard
        if (T) {
          file <- paste0("results/1.tab/40.1.CS_gene_body_percentage_clustered_by_kmeans_", data_type, "s.merged.", type, ".csv")
          cluster_golden <- read.table(file, header = T, sep = ",", fill = T, comment.char = "")
        }
        
        # using kmeans to cluster emission table
        if (T) {
          cluster_emission_kmeans <- read.table("results/1.tab/40.emission_kmeans_cluster.csv", header = T, sep = ",")
        }
        
        # read the chromIDEAS cluster
        if (T) {
          cluster_chrom <- qread("data/saved_data/29.common_cluster_name_for_merged_2cells_CS_clusters.qs", nthreads = 6)
        }
        
        # get random group as control
        if (T) {
          file <- paste0("data/saved_data/40.random_cluster_control_n", n_ct, ".merged.qs")
          cluster_random <- qread(file, nthreads = 6)
        }
        
        # make the order of all clusters are identical
        if (T) {
          cluster_emission_kmeans <- cluster_emission_kmeans[match(cluster_golden$state, cluster_emission_kmeans$state), ]
          cluster_chrom <- cluster_chrom[match(cluster_golden$state, cluster_chrom$orig.ident), ]
          cluster_random <- cluster_random[match(cluster_golden$state, cluster_random$state), ]
          
          cluster_compare <- cbind(cluster_emission_kmeans, cluster_chrom, cluster_random)
          cluster_compare <- cluster_compare[, c("cluster5", paste0("ct", 1:n_ct), "seurat_clusters")]
          colnames(cluster_compare) <- c("kmeans5", paste0("CT", 1:n_ct), "chromIDEAS")
          
          rm(cluster_emission_kmeans, cluster_chrom, cluster_random)
        }
        
        # compare the similarity compared with auc distribution line cluster
        if (T) {
          # setting
          if (T) {
            file <- paste0(dir2, "/40.1.CS_gene_body_percentage_ARI_in_", data_type, "s.merged.", type, ".pdf")
            pdf(file = file, width = 8, height = 6)
            colors <- c("#00AFBB", "#E7B800", "#FC4E07", "#a7c957")
            par(mfrow = c(2,2))
          }
          
          # statistics ari
          if (T) {
            ari_stat <- lapply(k_auc, function(kk) {
              similarity <- similarity_between_clusters(cluster_golden[, paste0("K", kk)], cluster_compare)
              ari <- similarity$ari
              names(ari) <- rownames(similarity)
              
              return(ari)
            })
            names(ari_stat) <- paste0("K", k_auc)
          }
          
          # plot
          if (T) {
            for (kk in k_auc) {
              # setting
              title <- paste0("State % in Gene Body cluster", kk)
              
              plot_radar_figure_with_similarity(ari_stat[[paste0("K", kk)]], 
                                                title, 
                                                colors[max(k_auc)+1-kk], 
                                                c(paste0("CT", c(3:5)), "chromIDEAS", "kmeans5", paste0("CT", c(1:2))))
            }
          }
          
          # save the results
          if (T) {
            ari <- data.frame(t(do.call(rbind, ari_stat)))
            ari$method <- rownames(ari)
            rownames(ari) <- 1:nrow(ari)
            
            file <- paste0("results/1.tab/40.1.CS_gene_body_percentage_ARI_in_", data_type, "s.merged.", type, ".csv")
            write.table(ari, file = file, quote = F, sep = ",", col.names = T, row.names = F)
          }
          
          Sys.sleep(0.05)
          dev.off()
        }
      }
    }
  }
}
