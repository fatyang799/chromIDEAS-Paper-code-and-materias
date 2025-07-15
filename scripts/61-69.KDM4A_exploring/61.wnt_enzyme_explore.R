# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(stringr)
  library(ggplot2)
  library(reshape2)
  library(qs)
  suppressPackageStartupMessages(library(edgeR))
}

# load the gene data
if (T) {
  file <- "data/saved_data/59.norm_expression_mat_thp1_cd34.qs"
  x <- qread(file, nthreads = 6)
}

# get cpm expression mat
if (T) {
  lcpm <- cpm(x, log=F)
  gene2id <- x$genes
}

# get eraser and writer related genes ID
if (T) {
  writer_k9me3 <- gene2id[grepl("PRDM2|KMT8A", gene2id$name, ignore.case = T) & gene2id$gene_type == "protein_coding", ]
  eraser_k9me3 <- gene2id[grepl("PHF8|KDM7B|JMJD2A|KDM4A|JMJD2B|KDM4B|JMJD2C|KDM4C|JMJD2D|KDM4D|KDM4E", gene2id$name, ignore.case = T) & gene2id$gene_type == "protein_coding", ]
  
  writer_k27me3 <- gene2id[grepl("EZH2", gene2id$name, ignore.case = T) & gene2id$gene_type == "protein_coding", ]
  eraser_k27me3 <- gene2id[grepl("UTX|KDM6A|JMJD3|KDM6B", gene2id$name, ignore.case = T) & gene2id$gene_type == "protein_coding", ]
}

# get deg
if (T) {
  file <- "results/1.tab/59.thp1_vs_cd34_mRNA_deg.csv"
  deg <- read.table(file, header = T, sep = ",", fill = T, comment.char = "")
}

# merge the info
if (T) {
  enzyme_dat <- data.frame(lcpm[c(writer_k9me3$id, eraser_k9me3$id, writer_k27me3$id, eraser_k27me3$id), ])
  
  enzyme_dat$id <- c(writer_k9me3$id, eraser_k9me3$id, writer_k27me3$id, eraser_k27me3$id)
  enzyme_dat$name <- c(writer_k9me3$name, eraser_k9me3$name, writer_k27me3$name, eraser_k27me3$name)
  enzyme_dat$state <- deg$state[match(enzyme_dat$id, deg$id)]
  enzyme_dat$pvalue <- deg$P.Value[match(enzyme_dat$id, deg$id)]
  enzyme_dat$log2fc <- deg$logFC[match(enzyme_dat$id, deg$id)]
  
  enzyme_dat$enzyme <- c(rep("writer", nrow(writer_k9me3)), rep("eraser", nrow(eraser_k9me3)), rep("writer", nrow(writer_k27me3)), rep("eraser", nrow(eraser_k27me3)))
  enzyme_dat$histone <- c(rep("H3K9me3", nrow(writer_k9me3)), rep("H3K9me3", nrow(eraser_k9me3)), rep("H3K27me3", nrow(writer_k27me3)), rep("H3K27me3", nrow(eraser_k27me3)))
  
  enzyme_dat <- enzyme_dat[order(enzyme_dat$enzyme),]
  
  rm(writer_k9me3, eraser_k9me3, writer_k27me3, eraser_k27me3, lcpm, x, deg)
}

# format the ggdat
if (T) {
  enzyme_dat$cd34 <- sapply(1:nrow(enzyme_dat), function(x) {
    sum(enzyme_dat[x, 1:2])/2
  })
  enzyme_dat$thp1 <- sapply(1:nrow(enzyme_dat), function(x) {
    sum(enzyme_dat[x, 3:4])/2
  })
  
  ggdat <- melt(enzyme_dat, id.vars = c("name",  "enzyme", "histone"), measure.vars = c("cd34", "thp1", "cd34_rep1", "cd34_rep2", "thp1_rep4", "thp1_rep5"), 
                variable.name = "Rep", value.name = "log2CPM", factorsAsStrings = F)
  head(ggdat)
  
  ggdat$name <- factor(ggdat$name, levels = enzyme_dat$name)
  ggdat$Cell <- str_split(ggdat$Rep, "_", simplify = T)[, 1]
}

# color setting
if (T) {
  library(RColorBrewer)
  display.brewer.all()
  colors <- brewer.pal(9, "Set1") 
  scales::show_col(colors, labels=T)
  
  col_cell <- colors[1:2]
  names(col_cell) <- c("cd34", "thp1")
}

# ggplot expression level
if (T) {
  head(ggdat)
  
  ggdat_mean <- ggdat[ggdat$Rep %in% c("cd34", "thp1"), ]
  ggdat_rep <- ggdat[!ggdat$Rep %in% c("cd34", "thp1"), ]
  
  p <- ggplot() +
    geom_bar(data=ggdat_mean, aes(x=name, y=log2CPM, fill=Rep), stat="identity", position="dodge") +
    geom_point(data=ggdat_rep, aes(x=name, y=log2CPM, fill=Cell), position=position_jitterdodge(jitter.width=0.2, dodge.width=1), size=2) +
    geom_hline(yintercept = 0, linetype=2) +
    facet_grid(histone ~ enzyme, scales = "free_x", space = "free") +
    ylab("Gene Expression (CPM)") +
    scale_fill_manual(values=col_cell) +
    scale_color_manual(name = "Rep", values = c("cd34" = "blue", "thp1" = "red")) +
    xlab(NULL) +
    cowplot::theme_cowplot() +
    theme(axis.title = element_text(size = rel(1.2)), 
          axis.text = element_text(size = rel(1.1)), 
          strip.text = element_text(size = rel(1.2)), 
          panel.border = element_rect(color="black"), 
          strip.background = element_rect(fill=NA, color=NA))
  
  file <- "results/2.pic/61.k9_27_me3_eraser_writer_enzyme_expression_level_logcpm.pdf"
  ggsave(filename = file, plot = p, width = 10, height = 8)
}

