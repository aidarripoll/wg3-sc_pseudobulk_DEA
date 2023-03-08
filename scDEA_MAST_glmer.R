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
  make_option(c("--covariates"), action="store", default='scDEA_covariates.tab', type='character',
              help="Covariates file."),
  make_option(c("-r","--residuals"), action="store", default=FALSE, type='logical',
              help="Collect the residuals."),
  make_option(c("-i", "--in_dir"), action="store", default='inputs', type='character',
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

# MAST version should be >= 1.16.0 --> with R version 4.0.3 (2020-10-10), we got MAST version 1.17.3
mast_version.user <- as.character(packageVersion("MAST"))
mast_version.required <- '1.16.0'
mast_cv <- compareVersion(mast_version.user, mast_version.required)
if(mast_cv!=1){
  err_message <- print(paste0('Your MAST version (', mast_version.user, ') is out-dated. It should be >= ', mast_version.required, 
                              '. Please up-grade it with: devtools::install_github("RGLab/MAST"). You can check the package version with: packageVersion("MAST").')) 
  stop(err_message)
}

# Input/Output directories
in.dir <- paste0(opt$in_dir, '/', opt$cell_level, '/')
out.dir <- paste0(opt$out_dir, '/', opt$cell_level, '/', opt$phenotype, '/', opt$cell_type, '/')
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

# Covariates 
covs.fn <- opt$covariates

# Seurat object
so.fn <- paste0(in.dir, opt$cell_type, '.Qced.Normalized.SCs.Rds')

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

# Read seurat object file
print(paste0('Reading PBMC seurat object file in: ',so.fn))
system.time(pbmc_i <- readRDS(so.fn))
DefaultAssay(pbmc_i) <- "RNA"

# Modify metadata
print('Modifying metadata...')
pbmc <- modify_metadata_sc(pbmc_i, opt$phenotype, covs.df)

# Subset seurat object (if needed --> if opt$phenotype==age_cat, keep only Y and O samples)
cells_initial <- ncol(pbmc)
cells_final <- cells_initial
print(paste0('Number of inital cells: ', cells_initial))
if(opt$phenotype=='age_cat'){
  print('Subsetting seurat object: keeping only Y and O samples...')
  pbmc <- subset(x = pbmc, subset = age_cat %in% c("Y","O"))
  cells_final <- ncol(pbmc)
}
pbmc@meta.data <- droplevels(pbmc@meta.data)
print(paste0('Number of final cells: ', cells_final))

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

################################ Pre-processing ################################
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

## Random effects
random_effects.vars <- covs.df[covs.df$type=='random',]$covariate

# Run sc-DEA
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