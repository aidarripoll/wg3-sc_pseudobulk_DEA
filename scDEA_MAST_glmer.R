#!/usr/bin/env Rscript

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("-l", "--cell_level"), action="store", default="L1", type='character',
              help="Cell type levels: Low resolution (L1) or high resolution (L2)"),
  make_option(c("-c", "--cell_type"), action="store", default=NULL, type='character',
              help="Cell types in the low cell level resolution (l1) or in the high cell level resolution (l2)"),
  make_option(c("-v", "--phenotype"), action="store", default=NULL, type='character',
              help="Variable/Phenotype (SEX/age/age_cat/age_cat_all/age_squared)"),
  make_option(c("--covariates"), action="store", default='scDEA.covariates.tab', type='character',
              help="Covariates file."),
  make_option(c("-r","--residuals"), action="store", default=FALSE, type='logical',
              help="Collect the residuals."),
  make_option(c("--combined_data"), action="store", default=FALSE, type='logical',
              help="If the data from more than one study has been previously merged in one object."),
  make_option(c("--genoPC1"), action="store", default=FALSE, type='logical',
              help="If the data contains more than one genetic ancestry."),
  make_option(c("-i", "--in_dir"), action="store", default='input', type='character',
              help="Input directory: WG2 output directory."),
  make_option(c("-o","--out_dir"), action="store", default="scDEA_MAST_glmer", type='character',
              help="Output directory.")
)
opt = parse_args(OptionParser(option_list=option_list))

################################# Set Variables ################################
 # Set working directory
cwd <- getwd()
setwd(cwd)

# Loading functions
functions.fn <- 'scripts/scDEA_functions.R'
print(paste0('Loading functions from: ', functions.fn))
source(functions.fn)

common_functions.fn <- 'scripts/common_functions.R'
print(paste0('Loading common functions from: ', common_functions.fn))
source(common_functions.fn)

# # MAST version should be >= 1.16.0 --> with R version 4.0.3 (2020-10-10), we got MAST version 1.17.3
# mast_version.user <- as.character(packageVersion("MAST"))
# mast_version.required <- '1.16.0'
# mast_cv <- compareVersion(mast_version.user, mast_version.required)
# if(mast_cv!=1){
#   err_message <- print(paste0('Your MAST version (', mast_version.user, ') is out-dated. It should be >= ', mast_version.required, 
#                               '. Please up-grade it with: devtools::install_github("RGLab/MAST"). You can check the package version with: packageVersion("MAST").')) 
#   stop(err_message)
# }

# Input/Output directories
in.dir <- paste0(opt$in_dir, '/', opt$cell_level, '/')
out.dir <- paste0(opt$out_dir, '/', opt$cell_level, '/', opt$phenotype, '/', opt$cell_type, '/')
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

# Covariates 
covs.fn <- opt$covariates

# Seurat object
so.fn <- paste0(in.dir, opt$cell_type, '.Qced.Normalized.SCs.Rds')

# smf file
smf.fn <- paste0(opt$in_dir, '/smf.txt')

# qtlInput Pcs file
qtlInput_Pcs.fn <- paste0(in.dir, opt$cell_type, '.qtlInput.Pcs.txt')

# Print report
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Cell type: ', opt$cell_type))
print(paste0('Phenotype: ', opt$phenotype))
print(paste0('Covariates file: ', opt$covariates))
print(paste0('Residuals: ', opt$residuals))
print(paste0('Input directory: ', opt$in_dir))
print(paste0('Input file: ', so.fn))
print(paste0('Output directory: ', out.dir))
cat('\n')

################################ Read input data ###############################
# Read covariates file
print(paste0('Reading covariates file in: ',covs.fn))
covs.df <- read.table(covs.fn, header = T)
covariates <- covs.df$covariate
phenotypes <- covs.df[covs.df$type=='fixed',]$covariate

# Read seurat object file
print(paste0('Reading PBMC seurat object file in: ',so.fn))
system.time(pbmc_i <- readRDS(so.fn))
DefaultAssay(pbmc_i) <- "RNA"

# Testing
# pbmc_raw <- pbmc_i
# set.seed(123)
# cells <- sample(colnames(pbmc_i), 5000)
# pbmc_i <- pbmc_i[,cells]

########################## Pre-processing (sc/psedobulk) ###########################
# Modify metadata
print('Modifying metadata...')
pbmc <- modify_metadata_sc(pbmc_i, covariates, covs.df)

# Filter donors (Donor_Pool) based on genotypes and QC filtering
print('Filter donors (Donor_Pool) based on genotypes and QC filtering...')
pbmc <- filter_donors(pbmc, smf.fn, qtlInput_Pcs.fn)

# Filter replicates: picking only one sample (Donor_Pool) per donor (Donor)
cell.md <- pbmc@meta.data
donor.md <- unique(cell.md[,c('Donor', 'Donor_Pool')])
if(any(table(donor.md$Donor)>1)){
  print('Filter replicates: picking only one sample (Donor_Pool) per donor (Donor)...')
  # Check
  donors <- unique(donor.md$Donor)
  ndonors <- length(donors)
  samples <- unique(donor.md$Donor_Pool)
  nsamples <- length(samples)
  donor_replicates.vec <- table(donor.md$Donor)[table(donor.md$Donor)>1]
  donor_replicates <- names(donor_replicates.vec)
  donor_replicates.n <- length(donor_replicates)
  print(paste0('# of initial Donor: ', ndonors))
  print(paste0('# of initial Donor_Pool: ', nsamples))
  print(paste0('# of initial Donor with replicates: ', donor_replicates.n))
  
  # Filtering samples
  samples_kept.list <- sapply(donors, function(i) filter_replicates(i, cell.md), simplify = FALSE)
  samples_kept <- unname(unlist(samples_kept.list))
  print(paste0('# of filtered Donor_Pool: ', length(samples_kept)))
  cell_md.filtered <- droplevels(cell.md[cell.md$Donor_Pool%in%samples_kept,])
  ncells.filtered <- nrow(cell_md.filtered)
  ndonors.filtered <- length(unique(cell_md.filtered$Donor))
  print(paste0('# of filtered Donor: ', ndonors.filtered, ' (nCells = ', ncells.filtered, ')'))
  cells_kept <- rownames(cell_md.filtered)
  
  # Filtering seurat object
  pbmc <- pbmc[,cells_kept]
  pbmc@meta.data <- cell_md.filtered
  rownames(pbmc@meta.data) <- pbmc@meta.data$Barcode
}

# (optional) Add a new variable for the previously merged studies
if(opt$combined_data){
  # Add a new variable in the seurat object
  pbmc$Dataset <- str_split_fixed(rownames(pbmc@meta.data), '_', 2)[,1]
  pbmc$Dataset <- as.factor(pbmc$Dataset)
  nDataset <- length(unique(pbmc$Dataset))
  print('The previous datasets with nCells were merged:')
  print(table(pbmc$Dataset))
  
  # Add a new covariate in the covariates file
  Dataset.df <- data.frame(covariate = 'Dataset',
                           type = 'random',
                           class = 'factor')
  covs.df <- rbind(covs.df, Dataset.df)
}

# (optional) Add a new variable for genoPC1
if(opt$genoPC1){
  print('The data contains more than one genetic ancestry...')
  
  # Read qtlInput_Pcs
  qtlInput_Pcs <- read.delim(qtlInput_Pcs.fn, check.names = FALSE)
  rownames(qtlInput_Pcs) <- qtlInput_Pcs[,1]
  qtlInput_Pcs <- qtlInput_Pcs[,-1]
  
  # Select only the filtered Donor_Pool samples
  Donor_Pool.vec <- unique(pbmc$Donor_Pool)
  qtlInput_Pcs.filt <- qtlInput_Pcs[rownames(qtlInput_Pcs)%in%Donor_Pool.vec,]
  qtlInput_Pcs.filt.df <- data.frame(Donor_Pool = rownames(qtlInput_Pcs.filt),
                                     genoPC1 = qtlInput_Pcs.filt$genoPC1)
  
  # Add a new variable in the seurat object
  pbmc@meta.data <- merge(pbmc@meta.data, qtlInput_Pcs.filt.df, by = 'Donor_Pool')
  rownames(pbmc@meta.data) <- pbmc@meta.data$Barcode
  
  # Add a new covariate in the covariates file
  genoPC1.df <- data.frame(covariate = 'genoPC1',
                           type = 'fixed',
                           class = 'integer')
  covs.df <- rbind(covs.df, genoPC1.df)
}

# Check covariates
print('Checking covariates...')
cov_vars <- covs.df$covariate
metadata_vars <- colnames(pbmc@meta.data)
if(!all(cov_vars%in%metadata_vars)){
  err_message <- paste0('Your metadata variables in: ', opt$covariates, ' should be in the seurat object metadata: ', so.fn)
  stop(err_message)
}

# Check phenotype
print('Checking phenotype...')
contrast_vars <- covs.df[covs.df$type=='fixed',]$covariate
contrast_vars.add <- c(contrast_vars, c('age_cat', 'age_cat_all', 'age_squared'))
if(!opt$phenotype%in%contrast_vars.add){
  err_message <- paste0('You should provide one of the following values in --phenotype argument: ', paste(contrast_vars.add, collapse = ', '))
}
cat('\n')

################################ Pre-processing (sc-specific) ################################
# Convert Seurat object to SCE object; and then from SCE object to SCA object
print('Converting Seurat object to Single Cell Assay object...')
system.time(sca_raw <- so_to_sca(pbmc))
sca_raw.genes <- nrow(sca_raw)

# Preprocessing: Filter out lowly variable genes --> CDR (=cngeneson; cellular detection rate) --> Filter out lowly expressed genes
print('Pre-processing Single Cell Assay object...')
system.time(sca <- preprocess_sca(sca_raw))
sca.genes <- nrow(sca)

# Save Seurat object (raw) and SCA objects (raw/filtered)
# Seurat object
pbmc_so.fn <- paste0(out.dir,'pbmc_so.rds')
print(paste0('Saving Seurat object (raw) in: ',pbmc_so.fn))
system.time(saveRDS(pbmc, pbmc_so.fn))

# SCA objects
## raw
sca_raw.fn <- paste0(out.dir,'pbmc_sca_raw.rds')
print(paste0('Saving Single Cell Assay object (raw) object in: ',sca_raw.fn))
system.time(saveRDS(sca_raw, sca_raw.fn))

## filtered
sca.fn <- paste0(out.dir,'pbmc_sca.rds')
print(paste0('Saving Single Cell Assay object (filtered) object in: ',sca.fn))
system.time(saveRDS(sca, sca.fn))

# Report for scDEA
genes_tested.prop <- round(sca.genes/sca_raw.genes,2)
print(paste0('Number of genes expressed: ', sca.genes, ' out of ', sca_raw.genes, ' genes'))
print(paste0('Proportion of genes expressed: ', genes_tested.prop))
ncells_tested <- ncol(sca)
print(paste0('Number of cells: ', ncells_tested))
cat('\n')

########################### sc-DEA with MAST GLMER #############################
# Set Fixed and Random effect variables
print('Defining the fixed and random effect variables...')

## Fixed effects
fixed_effects.covs <- covs.df[covs.df$type=='fixed',]$covariate
fixed_effects.vars <- c('cngeneson', fixed_effects.covs)
contrast.var_out <- opt$phenotype
if(opt$phenotype%in%c('age_cat', 'age_cat_all')){
  contrast.var_out <- 'age'
}
fixed_effects.vars <- fixed_effects.vars[!fixed_effects.vars%in%contrast.var_out]
fixed_effects.keep <- unlist(sapply(fixed_effects.vars, function(i) drop_terms(i, sca), simplify=FALSE))
fixed_effects.vars <- fixed_effects.vars[fixed_effects.keep]
  
## Random effects
random_effects.vars <- covs.df[covs.df$type=='random',]$covariate
random_effects.keep <- unlist(sapply(random_effects.vars, function(i) drop_terms(i, sca), simplify=FALSE))
random_effects.vars <- random_effects.vars[random_effects.keep]

# Run sc-DEA
# Testing
# sca_raw <- sca
# set.seed(123)
# genes <- sample(rownames(sca), 20)
# sca <- sca[genes,]

print('Running the sc-DEA with MAST glmer...')
system.time(de_glmer.res <- de_glmer.func(sca_object = sca,
                                          contrast = opt$phenotype,
                                          fixed_effects = fixed_effects.vars,
                                          random_effects = random_effects.vars,
                                          res = opt$residuals,
                                          out_dir = out.dir))
# Saving sc-DEA results
de_glmer.fn <- paste0(out.dir, 'de_glmer_nagq0.rds')
print(paste0('Saving sc-DEA results with MAST glmer in: ', de_glmer.fn))
saveRDS(de_glmer.res, de_glmer.fn)
cat('\n')

################################## Get DEGs ####################################
# Get sc-DEGs from MAST glmer output
print('Getting the sc-DEGs...')
degs_df <- get_degs(de_glmer.res, opt$phenotype)
print(degs_df)

# Saving sc-DEGs
get_degs.fn <- paste0(out.dir, 'de_glmer_nagq0.degs.rds')
print(paste0('Saving sc-DEGs results with MAST glmer in: ', get_degs.fn))
saveRDS(degs_df, get_degs.fn)
