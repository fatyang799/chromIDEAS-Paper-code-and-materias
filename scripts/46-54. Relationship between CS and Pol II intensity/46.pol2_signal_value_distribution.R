# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(ggplot2)
  library(stringr)
  library(reshape2)
  library(qs)
}

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
up_bin_num <- 5
down_bin_num <- 5
body_bin_num <- 10
rna_levels <- 10
state_file <- "data/raw_data/2.states/chromIDEAS.state"

# read the states data
if (T) {
  state <- data.table::fread(state_file, sep = " ", header = T, data.table = F)
  state <- state[, c(1, 5, 6)]
  colnames(state)[1] <- "ID"
  head(state)
}

# define the funtion
if (T) {
  # get cell specific state pol2 signal
  cell_specific_state_pol2_signal <- function(cell, state, mk) {
    # test dat
    if (F) {
      cell <- cell1
      mk <- "pol2"
    }
    
    # get cell specific pol2 data
    if (T) {
      path <- ifelse(mk == "pol2", 
                     paste0("data/raw_data/6.pol2/"), 
                     paste0("data/raw_data/7.histone/"))
      file <- paste0(path, "bs", bin_size, "/", cell, ".", mk, ".S3V2.bedgraph.NBP.txt")
      file <- ifelse(mk == "pol2", file, paste0(file, ".gz"))
      
      pol2 <- data.table::fread(file, sep = "\t", header = F, data.table = F)[, 1]
    }
    
    # merge the state and pol2 data
    if (T) {
      ggdat <- data.frame(state[, grep(paste0("ID|", cell), colnames(state))], 
                          pol2 = pol2)
      rm(pol2)
    }
    
    return(ggdat)
  }
  plot_state_pol2_signal <- function(ggdat, title, ylab="signal") {
    # test dat
    if (F) {
      title <- "title"
      ylab <- "signal"
    }
    
    # get boxplot statistics
    if (T) {
      boxstat <- tapply(ggdat$pol2, ggdat$state, summary)
      boxstat <- data.frame(do.call(rbind, boxstat))
      colnames(boxstat) <- c("ymin", "lower", "middle", "mean", "upper", "ymax")
      
      boxstat$state <- rownames(boxstat)
    }
    
    # plot the basic figure
    if (T) {
      head(ggdat)
      
      p <- ggplot() +
        geom_hline(yintercept = 0) +
        geom_violin(data = ggdat, aes(x=state, y=pol2, group=state, fill=state), width=1, scale = "width", na.rm=T) +
        geom_boxplot(data = boxstat, aes(x = state, ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax), width=0.5, fill=NA, color="black", alpha=0.5, stat = "identity", na.rm=T) + 
        ylab(ylab) +
        ggtitle(title) + 
        cowplot::theme_cowplot() +
        theme(legend.position = "none",
              strip.background = element_rect(fill = NA))
    }
    
    return(p)
  }
}

# the mk signal within state location across 2 cell data
if (T) {
  mks <- c("pol2", "ATAC", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K79me2", "H3K9me3")
  
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir <- paste0("results/2.pic/46.mk_signal_distribution/", data_type)
      if (! dir.exists(dir)) {
        dir.create(dir, showWarnings = F, recursive = T)
      }
    }
    
    for (type in c("merged", "single")) {
      cat(paste0("\t", type, ": \n"))
      
      # cell specific data
      if (type == "single") {
        for (cell in c(cell1, cell2)) {
          cat(paste0("\t\t", cell, ": \n"))
          
          for (mk in mks) {
            cat(paste0("\t\t\t", mk, ": \n"))
            
            # test data
            if (F) {
              data_type <- "tx"
              type <- "single"
              cell <- cell1
              mk <- "H3K27ac"
            }
            
            # read cell specific signal
            if (T) {
              ggdat <- cell_specific_state_pol2_signal(cell, state, mk)
            }
            
            # format the data
            if (T) {
              colnames(ggdat)[2] <- "state"
              ggdat$cell <- cell
              ggdat$state <- factor(ggdat$state, levels = sort(as.numeric(unique(ggdat$state))))
            }
            
            # plot pol2 distribution based on chromatin state
            if (T) {
              head(ggdat)
              
              p <- plot_state_pol2_signal(ggdat, 
                                          title=paste0(toupper(mk), " signal at each Chromatin State Location"), 
                                          ylab=paste0(toupper(mk), " Signal"))
              file <- paste0(dir, "/46.", mk, "_signal_within_CS.", type, ".", cell, ".jpeg")
              ggsave(filename = file, plot = p, width = 15, height = 10)
            }
          }
        }
      }
      
      # merged data
      if (type == "merged") {
        cat(paste0("\t\tmerged: \n"))
        
        for (mk in mks) {
          cat(paste0("\t\t\t", mk, ": \n"))
          
          # test data
          if (F) {
            data_type <- "tx"
            type <- "merged"
            mk <- "H3K27ac"
          }
          
          # read and merge signal
          if (T) {
            # read 2 cell pol2 data
            if (T) {
              dat1 <- cell_specific_state_pol2_signal(cell1, state, mk)
              dat2 <- cell_specific_state_pol2_signal(cell2, state, mk)
              
              dat1$cell <- cell1
              dat2$cell <- cell2
            }
            
            # merge the data
            if (T) {
              colnames(dat1)[2] <- colnames(dat2)[2] <- "state"
              
              ggdat <- rbind(dat1, dat2)
              rm(dat1, dat2)
            }
          }
          
          # format the data
          if (T) {
            ggdat$state <- factor(ggdat$state, levels = sort(as.numeric(unique(ggdat$state))))
          }
          
          # plot pol2 distribution based on chromatin state
          if (T) {
            head(ggdat)
            
            p <- plot_state_pol2_signal(ggdat, 
                                        title=paste0(toupper(mk), " signal at each Chromatin State Location"), 
                                        ylab=paste0(toupper(mk), " Signal"))
            file <- paste0(dir, "/46.", mk, "_signal_within_CS.", type, ".merged.jpeg")
            ggsave(filename = file, plot = p, width = 15, height = 10)
          }
        }
      }
    }
  }
}

# 这里single是每个细胞单独信号的分布图
# 这里merged是2个细胞信号合并在一起绘制的一个分布图
# 需要综合2个细胞画出一条曲线