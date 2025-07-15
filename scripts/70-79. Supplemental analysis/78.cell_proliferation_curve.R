# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(stringr)
  library(ggplot2)
  library(reshape2)
}

# KD_KDM4A_THP1
if (T) {
  # read the data
  if (T) {
    file <- "data/raw_data/11.experiments/KD_KDM4A_THP1.csv"
    dat <- read.table(file, header = T, sep = ",", fill = T, comment.char = "")
  }
  
  # format the data
  if (T) {
    head(dat)
    dat <- melt(dat, id.vars = "Cell", variable.name = "Time", value.name = "CellCount")
    
    summary_dat1 <- tapply(dat$CellCount, paste0(dat$Cell, "@", dat$Time), mean)
    summary_dat2 <- tapply(dat$CellCount, paste0(dat$Cell, "@", dat$Time), sd)
    summary_dat3 <- tapply(dat$CellCount, paste0(dat$Cell, "@", dat$Time), length)
    
    if (identical(names(summary_dat1), names(summary_dat2)) & identical(names(summary_dat1), names(summary_dat3))) {
      summary_dat <- data.frame(
        id = names(summary_dat1), 
        mean = summary_dat1, 
        sd = summary_dat2, 
        se = summary_dat2/sqrt(summary_dat3)
      )
    }
    
    summary_dat$Cell <- str_split(summary_dat$id, "@", simplify = T)[, 1]
    summary_dat$Time <- str_split(summary_dat$id, "@", simplify = T)[, 2]
    
    summary_dat$Cell <- factor(summary_dat$Cell, levels = unique(summary_dat$Cell))
    summary_dat$Time <- factor(summary_dat$Time, levels = paste0("Day", 0:7))
  }
  
  # ggplot
  if (T) {
    p <- ggplot(summary_dat, aes(x=Time, y=mean, color=Cell, group=Cell)) +
      geom_line(linewidth=1.2) +
      geom_point(size=2.5) +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, linewidth=1) +
      labs(x="Time (days)", y="Cell Count", color="Cell Type") +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=16), 
            legend.text=element_text(size=12),
            legend.title=element_text(size=14),
            panel.background=element_rect(color="black", fill=NA), 
            panel.grid=element_line(color="#EBEBEB", linewidth=0.6))
  }
  
  file <- "results/2.pic/78.cell_proliferation_curve_kd_kdm4a_thp1.pdf"
  ggsave(file, width = 8, height = 6)
  
  # statistics
  if (T) {
    tapply(dat$CellCount, paste0(dat$Cell, "@", dat$Time), mean)
    p <- lapply(split(dat, dat$Time), function(subdat) {
      # subdat <- dat[dat$Time == "Day7", ]
      
      if (all(subdat$Time != "Day0")) {
        cells <- c("WT", "KDM4A_1", "KDM4A_2", "KDM4A_3")
        sapply(cells, function(cell) {
          t.test(subdat$CellCount[subdat$Cell == cell], subdat$CellCount[subdat$Cell == "EV"])$p.value
        })
      }
    })
    
    do.call(cbind, p)
  }
}

# KD_KDM4A_THP1 OE WNT10B
if (T) {
  # read the data
  if (T) {
    file <- "data/raw_data/11.experiments/KD_KDM4A_OE_Wnt10b.csv"
    dat <- read.table(file, header = T, sep = ",", fill = T, comment.char = "")
  }
  
  # filter the data
  if (T) {
    dat <- dat[!grepl("WNT11|KD_merge", dat$Cell), ]
  }
  
  # read kd data
  if (T) {
    dat1 <- read.table("data/raw_data/11.experiments/KD_KDM4A_THP1.csv", header = T, sep = ",", fill = T, comment.char = "")
    dat1 <- dat1[grepl("KDM4A_", dat1$Cell), ]
    
    dat <- rbind(dat, dat1)
  }
  
  # format the data
  if (T) {
    head(dat)
    dat <- melt(dat, id.vars = "Cell", variable.name = "Time", value.name = "CellCount")
    
    summary_dat1 <- tapply(dat$CellCount, paste0(dat$Cell, "@", dat$Time), mean)
    summary_dat2 <- tapply(dat$CellCount, paste0(dat$Cell, "@", dat$Time), sd)
    summary_dat3 <- tapply(dat$CellCount, paste0(dat$Cell, "@", dat$Time), length)
    
    if (identical(names(summary_dat1), names(summary_dat2)) & identical(names(summary_dat1), names(summary_dat3))) {
      summary_dat <- data.frame(
        id = names(summary_dat1), 
        mean = summary_dat1, 
        sd = summary_dat2, 
        se = summary_dat2/sqrt(summary_dat3)
      )
    }
    
    summary_dat$Cell <- str_split(summary_dat$id, "@", simplify = T)[, 1]
    summary_dat$Time <- str_split(summary_dat$id, "@", simplify = T)[, 2]
    
    summary_dat$Cell <- factor(summary_dat$Cell, levels = unique(summary_dat$Cell))
    summary_dat$Time <- factor(summary_dat$Time, levels = paste0("Day", 0:7))
  }
  
  # ggplot
  if (T) {
    p <- ggplot(summary_dat, aes(x=Time, y=mean, color=Cell, group=Cell)) +
      geom_line(linewidth=1.2) +
      geom_point(size=2.5) +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, linewidth=1) +
      labs(x="Time (days)", y="Cell Count", color="Cell Type") +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=16), 
            legend.text=element_text(size=12),
            legend.title=element_text(size=14),
            panel.background=element_rect(color="black", fill=NA), 
            panel.grid=element_line(color="#EBEBEB", linewidth=0.6))
  }
  
  file <- "results/2.pic/78.cell_proliferation_curve_kd_kdm4a_oe_wnt10b_thp1.pdf"
  ggsave(file, width = 8, height = 6)
  
  # statistics
  if (T) {
    p <- lapply(split(dat, dat$Time), function(subdat) {
      # subdat <- dat[dat$Time == "Day7", ]
      
      if (all(subdat$Time != "Day0")) {
        cells <- c("WNT10B", "KDM4A KD + WNT10B", "KDM4A_1", "KDM4A_2", "KDM4A_3")
        sapply(cells, function(cell) {
          t.test(subdat$CellCount[subdat$Cell == cell], subdat$CellCount[subdat$Cell == "EV"])$p.value
        })
      }
    })
    
    do.call(cbind, p)
  }
}

