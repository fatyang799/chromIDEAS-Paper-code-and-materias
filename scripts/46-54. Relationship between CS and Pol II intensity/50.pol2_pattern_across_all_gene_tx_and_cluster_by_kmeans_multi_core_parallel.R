library(stringr)
library(reshape2)
library(qs)
library(doParallel)
library(parallel)
library(foreach)
cl <- makeCluster(20, outfile="log50.txt")
registerDoParallel(cl)

cell1 <- "thp1"
cell2 <- "cd34"
bin_size <- 200
up_bin_num <- 5
down_bin_num <- 5
body_bin_num <- 10
rna_levels <- 10

# define the funtion
if (T) {
  # common value
  if (T) {
    name <- c(
      paste0("U", up_bin_num:1), 
      "TSS", 
      paste0("B", 1:body_bin_num), 
      "TES", 
      paste0("D", 1:down_bin_num)
    )
    
    labels <- c(paste0("U", up_bin_num), 
                "TSS", 
                "TES", 
                paste0("D", down_bin_num))
  }
  
  # k-means (modify based on script40)
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
}

# data prepare: cell specific tss/body matrix for each gene
kmeans_single_row_dat <- function(k) {
  # load the environment
  if (T) {
    library(stringr)
    library(reshape2)
    library(qs)
  }
  
  # common value
  if (T) {
    cell1 <- "thp1"
    cell2 <- "cd34"
    bin_size <- 200
    up_bin_num <- 5
    down_bin_num <- 5
    body_bin_num <- 10
    rna_levels <- 10
  }
  
  # define the funtion
  if (T) {
    # common value
    if (T) {
      name <- c(
        paste0("U", up_bin_num:1), 
        "TSS", 
        paste0("B", 1:body_bin_num), 
        "TES", 
        paste0("D", 1:down_bin_num)
      )
      
      labels <- c(paste0("U", up_bin_num), 
                  "TSS", 
                  "TES", 
                  paste0("D", down_bin_num))
    }
    
    # k-means (modify based on script40)
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
  }
  
  # loop body
  if (T) {
    dat <- cluster_kmeans(mat[, name], k, scale="row")
    file <- paste0(dir1, "/50.", data_type, ".", type, ".", cell, ".", k, ".row.qs")
    qsave(dat, file)
  }
}
kmeans_merged_row_dat <- function(k) {
  # load the environment
  if (T) {
    library(stringr)
    library(reshape2)
    library(qs)
  }
  
  # common value
  if (T) {
    cell1 <- "thp1"
    cell2 <- "cd34"
    bin_size <- 200
    up_bin_num <- 5
    down_bin_num <- 5
    body_bin_num <- 10
    rna_levels <- 10
  }
  
  # define the funtion
  if (T) {
    # common value
    if (T) {
      name <- c(
        paste0("U", up_bin_num:1), 
        "TSS", 
        paste0("B", 1:body_bin_num), 
        "TES", 
        paste0("D", 1:down_bin_num)
      )
      
      labels <- c(paste0("U", up_bin_num), 
                  "TSS", 
                  "TES", 
                  paste0("D", down_bin_num))
    }
    
    # k-means (modify based on script40)
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
  }
  
  # loop body
  if (T) {
    dat <- cluster_kmeans(mat[, name], k, scale="row")
    file <- paste0(dir1, "/50.", data_type, ".", type, ".merged.", k, ".row.qs")
    qsave(dat, file)
  }
}
kmeans_single_col_dat <- function(k) {
  # load the environment
  if (T) {
    library(stringr)
    library(reshape2)
    library(qs)
  }
  
  # common value
  if (T) {
    cell1 <- "thp1"
    cell2 <- "cd34"
    bin_size <- 200
    up_bin_num <- 5
    down_bin_num <- 5
    body_bin_num <- 10
    rna_levels <- 10
  }
  
  # define the funtion
  if (T) {
    # common value
    if (T) {
      name <- c(
        paste0("U", up_bin_num:1), 
        "TSS", 
        paste0("B", 1:body_bin_num), 
        "TES", 
        paste0("D", 1:down_bin_num)
      )
      
      labels <- c(paste0("U", up_bin_num), 
                  "TSS", 
                  "TES", 
                  paste0("D", down_bin_num))
    }
    
    # k-means (modify based on script40)
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
  }
  
  # loop body
  if (T) {
    dat <- cluster_kmeans(mat[, name], k, scale="col")
    file <- paste0(dir1, "/50.", data_type, ".", type, ".", cell, ".", k, ".col.qs")
    qsave(dat, file)
  }
}
kmeans_merged_col_dat <- function(k) {
  # load the environment
  if (T) {
    library(stringr)
    library(reshape2)
    library(qs)
  }
  
  # common value
  if (T) {
    cell1 <- "thp1"
    cell2 <- "cd34"
    bin_size <- 200
    up_bin_num <- 5
    down_bin_num <- 5
    body_bin_num <- 10
    rna_levels <- 10
  }
  
  # define the funtion
  if (T) {
    # common value
    if (T) {
      name <- c(
        paste0("U", up_bin_num:1), 
        "TSS", 
        paste0("B", 1:body_bin_num), 
        "TES", 
        paste0("D", 1:down_bin_num)
      )
      
      labels <- c(paste0("U", up_bin_num), 
                  "TSS", 
                  "TES", 
                  paste0("D", down_bin_num))
    }
    
    # k-means (modify based on script40)
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
  }
  
  # loop body
  if (T) {
    dat <- cluster_kmeans(mat[, name], k, scale="col")
    file <- paste0(dir1, "/50.", data_type, ".", type, ".merged.", k, ".col.qs")
    qsave(dat, file)
  }
}
if (T) {
  for (data_type in c("tx")) {
    cat(paste0(data_type, ": \n"))
    
    # mkdir
    if (T) {
      dir1 <- paste0("data/saved_data/50.pol2_pattern/", data_type)
      if (! dir.exists(dir1)) {
        dir.create(dir1, showWarnings = F, recursive = T)
      }
      
      dir2 <- paste0("results/2.pic/50.pol2_pattern/", data_type)
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
          
          # get pol2 matrix: mat
          if (T) {
            cell_dat <- function(cell, type) {
              input <- paste0("data/saved_data/48.mk_signal_distribution/tx/48.", data_type, "_Body_", up_bin_num, "UP_", 
                              down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".pol2.mean_mat.", cell, ".", type, ".qs")
              
              mat <- qread(input, nthreads = 6)
              rownames(mat) <- paste0(mat$gene_id, "@", cell)
              
              return(mat)
            }
            mat <- cell_dat(cell, type)
          }
          
          # get kmeans clusters: dat
          foreach(k = 1:15, .export =ls()) %dopar% kmeans_single_row_dat(k)
          foreach(k = 1:15, .export =ls()) %dopar% kmeans_single_col_dat(k)
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
        
        # get pol2 matrix: mat
        if (T) {
          cell_dat <- function(cell, type) {
            input <- paste0("data/saved_data/48.mk_signal_distribution/tx/48.", data_type, "_Body_", up_bin_num, "UP_", 
                            down_bin_num, "DW_", body_bin_num, "Body.bs", bin_size, ".pol2.mean_mat.", cell, ".", type, ".qs")
            
            mat <- qread(input, nthreads = 6)
            rownames(mat) <- paste0(mat$gene_id, "@", cell)
            
            return(mat)
          }
          dat1 <- cell_dat(cell1, type="single")
          dat2 <- cell_dat(cell2, type="single")
          
          mat <- rbind(dat1, dat2)
        }
        
        # get kmeans clusters: dat
        foreach(k = 1:15, .export =ls()) %dopar% kmeans_merged_row_dat(k)
        foreach(k = 1:15, .export =ls()) %dopar% kmeans_merged_col_dat(k)
      }
    }
  }
}

stopCluster(cl)