# load the environment
if (T) {
  rm(list = ls())
  options(stringAsFactors = F)
  library(stringr)
  library(qs)
}

# annotation check
if (F) {
  file <- "data/raw_data/5.RNA/cd34/gene.info.txt.gz"
  cd34 <- data.table::fread(file, header = T, sep = "\t", data.table = F)
  file <- "data/raw_data/5.RNA/thp1/gene.info.txt.gz"
  thp1 <- data.table::fread(file, header = T, sep = "\t", data.table = F)
  identical(cd34, thp1)
  
  file <- "data/raw_data/5.RNA/cd34/tx.info.txt.gz"
  cd34 <- data.table::fread(file, header = T, sep = "\t", data.table = F)
  file <- "data/raw_data/5.RNA/thp1/tx.info.txt.gz"
  thp1 <- data.table::fread(file, header = T, sep = "\t", data.table = F)
  identical(cd34, thp1)
}

# prepare the matrix: TX
if (T) {
  # count
  if (file.exists("data/saved_data/4.tx_count_matrix_thp1_cd34.qs") + 
      file.exists("data/saved_data/4.tx_geneid_thp1_cd34.qs") < 2) {
    # get all files
    if (T) {
      file <- "data/raw_data/5.RNA/cd34/mRNA_matrix_RSEM.tx.count.txt"
      cd34 <- data.table::fread(file, header = T, sep = "\t", data.table = F)
      file <- "data/raw_data/5.RNA/thp1/mRNA_matrix_RSEM.tx.count.txt"
      thp1 <- data.table::fread(file, header = T, sep = "\t", data.table = F)
      if (identical(cd34$ID, thp1$ID)) {
        dat <- cbind(cd34, thp1[, -1])
        rm(cd34, thp1)
      }
      
      file <- "data/raw_data/5.RNA/cd34/tx.info.txt.gz"
      geneid <- data.table::fread(file, header = T, sep = "\t", data.table = F)
    }
    
    # format the data
    if (T) {
      # rownames
      if (T) {
        # check the duplication for rownames
        if (T) {
          gene_id <- dat[,1]
          message <- ifelse(sum(duplicated(gene_id))==0, 
                            "No duplication for 1st column id with dot",
                            paste0("Have ",sum(sum(duplicated(gene_id)))," duplication for 1st column id with dot"))
          print(message)
          gene_id <- str_split(gene_id, "[.]", simplify = T)[,1]
          message <- ifelse(sum(duplicated(gene_id))==0, 
                            "No duplication for 1st column id with dot",
                            paste0("Have ",sum(sum(duplicated(gene_id)))," duplication for 1st column id with dot"))
          print(message)
        }
        
        # change the rownames
        if (T) {
          print("Using the 1st column name as rownames of matrix")
          rownames(dat) <- dat[,1]
        }
        
        dat <- dat[, -1]
      }
      
      # colnames
      if (T) {
        colnames(dat)
        colnames(dat) <- gsub("_mRNA|_PE", "",colnames(dat))
        colnames(dat)
      }
      
      # rounding
      if (T) {
        dat <- round(dat)
      }
      
      # make sure the order
      if (T) {
        all(geneid$tx_id %in% rownames(dat))
        all(rownames(dat) %in% geneid$tx_id)
        geneid <- geneid[geneid$tx_id %in% rownames(dat), ]
        
        dat <- dat[match(geneid$tx_id, rownames(dat)), ]
        identical(rownames(dat), geneid$tx_id)
      }
    }
    
    # save data
    if (T) {
      file <- "data/saved_data/4.tx_count_matrix_thp1_cd34.qs"
      qsave(dat, file = file, nthreads = 6)
      
      file <- "data/saved_data/4.tx_geneid_thp1_cd34.qs"
      qsave(geneid, file = file, nthreads = 6)
    }
  }
  
  # fpkm
  if (! file.exists("data/saved_data/4.tx_fpkm_matrix_thp1_cd34.qs")) {
    # get all files
    if (T) {
      file <- "data/raw_data/5.RNA/cd34/mRNA_matrix_RSEM.tx.fpkm.txt"
      cd34 <- data.table::fread(file, header = T, sep = "\t", data.table = F)
      file <- "data/raw_data/5.RNA/thp1/mRNA_matrix_RSEM.tx.fpkm.txt"
      thp1 <- data.table::fread(file, header = T, sep = "\t", data.table = F)
      if (identical(cd34$ID, thp1$ID)) {
        dat <- cbind(cd34, thp1[, -1])
        rm(cd34, thp1)
      }
      
      file <- "data/saved_data/4.tx_geneid_thp1_cd34.qs"
      geneid <- qread(file)
    }
    
    # format the data
    if (T) {
      # rownames
      if (T) {
        # check the duplication for rownames
        if (T) {
          gene_id <- dat[,1]
          message <- ifelse(sum(duplicated(gene_id))==0, 
                            "No duplication for 1st column id with dot",
                            paste0("Have ",sum(sum(duplicated(gene_id)))," duplication for 1st column id with dot"))
          print(message)
          gene_id <- str_split(gene_id, "[.]", simplify = T)[,1]
          message <- ifelse(sum(duplicated(gene_id))==0, 
                            "No duplication for 1st column id with dot",
                            paste0("Have ",sum(sum(duplicated(gene_id)))," duplication for 1st column id with dot"))
          print(message)
        }
        
        # change the rownames
        if (T) {
          print("Using the 1st column name as rownames of matrix")
          rownames(dat) <- dat[,1]
        }
        
        dat <- dat[, -1]
      }
      
      # colnames
      if (T) {
        colnames(dat)
        colnames(dat) <- gsub("_mRNA|_PE", "",colnames(dat))
        colnames(dat)
      }
      
      # make sure the order
      if (T) {
        all(geneid$tx_id %in% rownames(dat))
        all(rownames(dat) %in% geneid$tx_id)
        identical(rownames(dat), geneid$tx_id)
        
        dat <- dat[match(geneid$tx_id, rownames(dat)), ]
        identical(rownames(dat), geneid$tx_id)
      }
    }
    
    # save data
    if (T) {
      file <- "data/saved_data/4.tx_fpkm_matrix_thp1_cd34.qs"
      qsave(dat, file = file, nthreads = 6)
    }
  }
}

# prepare the matrix: Gene
if (T) {
  # count
  if (file.exists("data/saved_data/4.gene_count_matrix_thp1_cd34.qs") + 
      file.exists("data/saved_data/4.gene_geneid_thp1_cd34.qs") < 2) {
    # get all files
    if (T) {
      file <- "data/raw_data/5.RNA/cd34/mRNA_matrix_RSEM.gene.count.txt"
      cd34 <- data.table::fread(file, header = T, sep = "\t", data.table = F)
      file <- "data/raw_data/5.RNA/thp1/mRNA_matrix_RSEM.gene.count.txt"
      thp1 <- data.table::fread(file, header = T, sep = "\t", data.table = F)
      if (identical(cd34$ID, thp1$ID)) {
        dat <- cbind(cd34, thp1[, -1])
        rm(cd34, thp1)
      }
      
      file <- "data/raw_data/5.RNA/cd34/gene.info.txt.gz"
      geneid <- data.table::fread(file, header = T, sep = "\t", data.table = F)
    }
    
    # format the data
    if (T) {
      # rownames
      if (T) {
        # check the duplication for rownames
        if (T) {
          gene_id <- dat[,1]
          message <- ifelse(sum(duplicated(gene_id))==0, 
                            "No duplication for 1st column id with dot",
                            paste0("Have ",sum(sum(duplicated(gene_id)))," duplication for 1st column id with dot"))
          print(message)
          gene_id <- str_split(gene_id, "[.]", simplify = T)[,1]
          message <- ifelse(sum(duplicated(gene_id))==0, 
                            "No duplication for 1st column id with dot",
                            paste0("Have ",sum(sum(duplicated(gene_id)))," duplication for 1st column id with dot"))
          print(message)
        }
        
        # change the rownames
        if (T) {
          print("Using the 1st column name as rownames of matrix")
          rownames(dat) <- dat[,1]
        }
        
        dat <- dat[, -1]
      }
      
      # colnames
      if (T) {
        colnames(dat)
        colnames(dat) <- gsub("_mRNA|_PE", "",colnames(dat))
        colnames(dat)
      }
      
      # rounding
      if (T) {
        dat <- round(dat)
      }
      
      # make sure the order
      if (T) {
        all(geneid$id %in% rownames(dat))
        all(rownames(dat) %in% geneid$id)
        geneid <- geneid[geneid$id %in% rownames(dat), ]
        
        dat <- dat[match(geneid$id, rownames(dat)), ]
        identical(geneid$id, rownames(dat))
      }
    }
    
    # save data
    if (T) {
      file <- "data/saved_data/4.gene_count_matrix_thp1_cd34.qs"
      qsave(dat, file = file, nthreads = 6)
      
      file <- "data/saved_data/4.gene_geneid_thp1_cd34.qs"
      qsave(geneid, file = file, nthreads = 6)
    }
  }
  
  # fpkm
  if (! file.exists("data/saved_data/4.gene_fpkm_matrix_thp1_cd34.qs")) {
    # get all files
    if (T) {
      file <- "data/raw_data/5.RNA/cd34/mRNA_matrix_RSEM.gene.fpkm.txt"
      cd34 <- data.table::fread(file, header = T, sep = "\t", data.table = F)
      file <- "data/raw_data/5.RNA/thp1/mRNA_matrix_RSEM.gene.fpkm.txt"
      thp1 <- data.table::fread(file, header = T, sep = "\t", data.table = F)
      if (identical(cd34$ID, thp1$ID)) {
        dat <- cbind(cd34, thp1[, -1])
        rm(cd34, thp1)
      }
      
      file <- "data/saved_data/4.gene_geneid_thp1_cd34.qs"
      geneid <- qread(file, nthreads = 6)
    }
    
    # format the data
    if (T) {
      # rownames
      if (T) {
        # check the duplication for rownames
        if (T) {
          gene_id <- dat[,1]
          message <- ifelse(sum(duplicated(gene_id))==0, 
                            "No duplication for 1st column id with dot",
                            paste0("Have ",sum(sum(duplicated(gene_id)))," duplication for 1st column id with dot"))
          print(message)
          gene_id <- str_split(gene_id, "[.]", simplify = T)[,1]
          message <- ifelse(sum(duplicated(gene_id))==0, 
                            "No duplication for 1st column id with dot",
                            paste0("Have ",sum(sum(duplicated(gene_id)))," duplication for 1st column id with dot"))
          print(message)
        }
        
        # change the rownames
        if (T) {
          print("Using the 1st column name as rownames of matrix")
          rownames(dat) <- dat[,1]
        }
        
        dat <- dat[, -1]
      }
      
      # colnames
      if (T) {
        colnames(dat)
        colnames(dat) <- gsub("_mRNA|_PE", "",colnames(dat))
        colnames(dat)
      }
      
      # make sure the order
      if (T) {
        all(geneid$id %in% rownames(dat))
        all(rownames(dat) %in% geneid$id)
        identical(geneid$id, rownames(dat))
        
        dat <- dat[match(geneid$id, rownames(dat)), ]
        identical(geneid$id, rownames(dat))
      }
    }
    
    # save data
    if (T) {
      file <- "data/saved_data/4.gene_fpkm_matrix_thp1_cd34.qs"
      qsave(dat, file = file, nthreads = 6)
    }
  }
}

