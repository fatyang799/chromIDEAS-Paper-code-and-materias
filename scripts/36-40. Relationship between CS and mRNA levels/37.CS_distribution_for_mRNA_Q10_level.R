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

# define functions
if (T) {
  plot_cell_specifc_CS_distribution_for_mRNA_Q10 <- function(input, scale, rms0, title) {
    # test dat
    if (F) {
      input <- "data/saved_data/36.tx_level_profile_Body_5UP_5DW_10Body.bs200.format_genomic_percentage_mRNA_classification_Q10_mat.cd34.qs"
      input <- "data/saved_data/36.tx_level_profile_TSS_5UP_5DW.bs200.format_genomic_percentage_mRNA_classification_Q10_mat.cd34.qs"
      scale <- "free"
      rms0 <- F
      title <- "test"
    }
    
    # get parameters
    if (T) {
      # type <- str_extract(basename(input), "TSS|Body")
      # cell <- str_extract(basename(input), paste0(cell1, "|", cell2))
      scale_type <- str_split(str_extract(basename(input), "genomic_percentage|state_percentage"), "_", simplify = T)[, 1]
    }
    
    # get distribution data
    if (T) {
      dat <- qread(input, nthreads = 6)
    }
    
    if (rms0) {
      dat <- dat[dat$State != "S0", ]
    }
    
    # format the data
    if (T) {
      colnames(dat) <- c("State", "Loc", "Percentage", "mRNA")
      ylabel <- ifelse(scale_type == "genomic", "Genomic_Percentage", 
                       ifelse(scale_type == "state", "State_Percentage", NA))
      
      dat$group <- paste0(dat$State, "@", dat$mRNA)
    }
    
    # color setting
    if (T) {
      col_fun <- colorRamp2(c(1, length(levels(dat$mRNA))), c("blue", "red"))
      cols <- col_fun(1:length(levels(dat$mRNA)))
      names(cols) <- levels(dat$mRNA)
    }
    
    # ggplot
    if (T) {
      head(dat)
      p <- ggplot(dat) +
        geom_line(aes(x=Loc, y=Percentage, group=group, color=mRNA), alpha=0.6, linewidth=1) +
        scale_x_discrete(name = NULL, breaks=c(paste0("U", up_bin_num), "TSS", "TES", paste0("D", down_bin_num))) +
        scale_color_manual(values = cols) +
        ylab(ylabel) +
        ggtitle(title) +
        facet_wrap(.~State, scales = scale, ncol = 5) +
        cowplot::theme_cowplot() +
        theme(legend.position = c(0.90, 0.05), 
              panel.border = element_rect(color="black"), 
              strip.background = element_rect(fill=NA, color=NA))
    }
    
    return(p)
  }
}

# plot the CS distribution
if (T) {
  data_type <- "tx"
  scale <- "free"
  rms0 <- F
  
  for (type in c("Body")) {
    cat(paste0(type, ":\n"))
    for (scale_type in c("genomic")) {
      cat(paste0("\t", scale_type, ":\n"))
      for (cell in c(cell1, cell2, "merged")) {
        cat(paste0("\t\t", cell, ": \n"))
        
        # test data
        if (F) {
          type <- "Body"
          scale_type <- "genomic"
          cell <- "merged"
        }
        
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
        
        # get output
        if (T) {
          file <- paste0("results/2.pic/37.CS_distribution_", type, "_", cell, "_", data_type, "_level_", scale_type, 
                         "_percentage_mRNA_Q", rna_levels, ".pdf")
        }
        
        # ggplot
        if (T) {
          title <- paste0(cell, ": ", data_type, "+", scale_type)
          p <- plot_cell_specifc_CS_distribution_for_mRNA_Q10(input, scale, rms0, title)
          
          ggsave(filename = file, plot = p, width = 10, height = 12)
        }
      }
    }
  }
}

