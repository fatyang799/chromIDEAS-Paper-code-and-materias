# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  suppressPackageStartupMessages(library(circlize))
  suppressPackageStartupMessages(library(ComplexHeatmap))
  library(RColorBrewer)
  library(ggplot2)
}

input <- "data/raw_data/1.emission/chromIDEAS.emission.txt"

# get the emission table
if (T) {
  dat <- read.table(input, header = T, sep = "\t", fill = T, comment.char = "")
}

# format the data
if (T) {
  # modify the rownames
  rownames(dat) <- dat[, 1]
  
  # 热图右侧注释标注每种染色质状态所占基因组比例，由于state0比例过高，故将其调整为基因组占比第2高
  if (T) {
    percentage <- dat[, 1:2]
    percentage$Percentage[percentage$Percentage>sort(percentage$Percentage, decreasing = T)[2]] <- sort(percentage$Percentage, decreasing = T)[2]
  }
  
  # format the matrix
  if (T) {
    dat <- dat[, -c(1:2)]
    dat <- as.matrix(dat)
  }
}

# color setting
if (T) {
  display.brewer.all()
  colors <- brewer.pal(9, "Set1") 
  scales::show_col(colors, labels=T)
}

# ordering
if (T) {
  ord_row <- hclust(d=dist(dat, method = "euclidean"), method = "ward.D2")$labels
  ord_col <- c("ATAC", "H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3", "H3K9me3", "H3K79me2")
  
  dat <- dat[ord_row, ]
  dat <- dat[, ord_col]
  
  qsave(list(ord_row, ord_col), file = "data/saved_data/2.order_of_emission_table.qs", nthreads = 6)
}

# plot complexheatmap (raw)
if (T) {
  col_fun = colorRamp2(c(0, ceiling(quantile(c(dat), 0.95))),
                       c("#FFFFFF", "#3953A4"))
  
  p <- Heatmap(dat, name = "Probablity", col = col_fun, 
               border_gp = gpar(col = "black"), rect_gp = gpar(col = "black"),
               show_heatmap_legend = TRUE, heatmap_legend_param = list(border = T, title_position = "lefttop", legend_direction = "horizontal"),
               row_names_max_width = unit(6, "cm"), row_names_gp = gpar(fontsize = 10), row_names_rot = 0,
               column_names_max_height = unit(6, "cm"), column_names_gp = gpar(fontsize = 10), column_names_rot = 60,
               # right_annotation = rowAnnotation(percentage=anno_barplot(percentage$Percentage, 
               #                                                          gp = gpar(fill = "#3953A4"),
               #                                                          bar_width = 0.8,
               #                                                          ylim = c(0,sort(percentage$Percentage, decreasing = T)[2]),
               #                                                          height = unit(2, "cm"))), 
               cluster_rows = F, cluster_columns = F, show_row_dend = F, show_column_dend = F)
  
  file <- "results/2.pic/2.emission_heatmap_raw.pdf"
  p <- ggplotify::as.ggplot(p)
  ggsave(filename = file, plot = p, width = 4.8, height = 10)
}

