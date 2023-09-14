#!/usr/bin/env Rscript

############################### Load R packages ################################
print('Loading R packages...')
shhh <- suppressPackageStartupMessages
shhh(library(Seurat))
shhh(library(MAST))
shhh(library(SingleCellExperiment))
shhh(library(dreamlet))
shhh(library(zenith))
shhh(library(scater))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(ggplot2))
shhh(library(RColorBrewer))

######################## Functions used in pseudobulkDEA_limmadream.R and scDEA_MAST_glmer.R ##################
# 1. Modify sc-metadata
modify_metadata_sc <- function(so, covs, covs_df){
  # Grab metadata
  metadata <- so@meta.data
  
  # Report
  cells_all <- nrow(metadata)
  print(paste0('Number of initial cells: ', cells_all))
  
  # Remove possible NAs
  for (i in covs){
    na.boolean <- any(is.na(metadata[[i]]))
    none.boolean <- any(metadata[[i]]==0)
    if(na.boolean | none.boolean){
      print(paste0('Removing NAs/none in the phenotype: ', i))
      na.idx <- which(is.na(metadata[[i]]))
      none.idx <- which(metadata[[i]]==0)
      na_none.idx <- c(na.idx, none.idx)
      metadata <- metadata[-na_none.idx,]
      cells_notna <- nrow(metadata)
      print(paste0('Number of annotated cells: ', cells_notna))
    }
  }
  
  # Create random factors variables in the metadata
  metadata <- metadata %>% separate(Donor_Pool, c('Donor', 'Pool'), ';;', remove = FALSE)
  
  # Phenotype modifications
  if('SEX'%in%covs){
    print(paste0('Relevel SEX...'))
    metadata[['SEX']] <-  ifelse(metadata[['SEX']]==1,'M','F')
    Group_order <- c('M','F')
    metadata[['SEX']] <- factor(metadata[['SEX']],
                                levels = Group_order)
  }
  
  if('age_cat'%in%covs){
    print(paste0('Creating a new metadata variable: age_cat...'))
    metadata[['age_cat']] <- ifelse(metadata$age<=40, 'Y',
                                    ifelse(metadata$age>=60, 'O', 'M'))
    metadata <- metadata[metadata[['age_cat']]%in%c('Y','O'),] #only keep Y and O samples (remove if we want to consider all samples)
    Group_order <- c('Y','M','O')
    metadata[['age_cat']] <- factor(metadata[['age_cat']],
                                    levels = Group_order)
  }
  
  if('age_cat_all'%in%covs){
    print(paste0('Creating a new metadata variable:  age_cat_all...'))
    metadata[['age_cat_all']] <- ifelse(metadata$age<=40, 'Y', 'O')
    Group_order <- c('Y','O')
    metadata[['age_cat_all']] <- factor(metadata[['age_cat_all']],
                                        levels = Group_order)
  }
  
  if('age_squared'%in%covs){
    print(paste0('Creating a new metadata variable: age_squared...'))
    metadata[['age_squared']] <- metadata$age^2
  }
  
  # Declare factors
  factor_vars <- covs_df[covs_df$class=='factor',]$covariate
  metadata %>% 
    mutate_at(all_of(factor_vars), as.factor) %>%
    as.data.frame() -> metadata
  metadata %>% 
    mutate_at(all_of(factor_vars), droplevels) %>%
    as.data.frame() -> metadata
  
  # Add new metadata
  cells_final <- nrow(metadata)
  print(paste0('Number of final cells: ', cells_final))
  cells_kept <- rownames(metadata)
  so <- so[,cells_kept]
  so@meta.data <- metadata
  so$Barcode <- rownames(so@meta.data)
  rownames(so@meta.data) <- so@meta.data$Barcode
  
  return(so)
}

# 2. Filter donors (Donor_Pool) based on genotypes and QC filtering
filter_donors <- function(so, smf_fn, qtlInput_Pcs_fn){
  cell_md <- so@meta.data
  ncells <- nrow(cell_md)
  ndonors <- length(unique(cell_md$Donor_Pool))
  qtlInput_Pcs <- read.delim(qtlInput_Pcs_fn, check.names = FALSE)
  smf <- read.delim(smf_fn, check.names = FALSE)
  rownames(qtlInput_Pcs) <- qtlInput_Pcs[,1]
  qtlInput_Pcs <- qtlInput_Pcs[,-1]
  smf.vec <- unique(smf$phenotype_id)
  qtlInput_Pcs.vec <- rownames(qtlInput_Pcs)
  donor_pool.vec <- intersect(smf.vec, qtlInput_Pcs.vec)
  cell_md.filtered <- droplevels(cell_md[cell_md$Donor_Pool%in%donor_pool.vec,])
  ncells.filtered <- nrow(cell_md.filtered)
  ndonors.filtered <- length(unique(cell_md.filtered$Donor_Pool))
  print(paste0('# of initial Donor_Pool: ', ndonors, ' (nCells = ', ncells, ')'))
  print(paste0('# of filtered Donor_Pool: ', ndonors.filtered, ' (nCells = ', ncells.filtered, ')'))
  cells_kept <- rownames(cell_md.filtered)
  so <- so[,cells_kept]
  so@meta.data <- cell_md.filtered
  rownames(so@meta.data) <- so@meta.data$Barcode
  return(so)
}

# 3. Filter replicates: pick only one sample (Donor_Pool) per donor (Donor)
filter_replicates <- function(donor, cell_md){
  cell_md.donor <- droplevels(cell_md[cell_md$Donor==donor,])
  Donor_Pool.tab <- table(cell_md.donor$Donor_Pool)
  Donor_Pool.vec <- as.vector(Donor_Pool.tab)
  names(Donor_Pool.vec) <- names(Donor_Pool.tab)
  ncells_max <- max(Donor_Pool.vec)
  Donor_Pool.max <- Donor_Pool.vec[Donor_Pool.vec==ncells_max]
  if(length(Donor_Pool.max)>1){
    Donor_Pool.max <- Donor_Pool.max[1]
  }
  Donor_Pool.kept <- names(Donor_Pool.max)
  return(Donor_Pool.kept)
}