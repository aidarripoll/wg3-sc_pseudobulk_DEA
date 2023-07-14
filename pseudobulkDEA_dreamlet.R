#!/usr/bin/env Rscript

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("-l", "--cell_level"), action="store", default="L1", type='character',
              help="Cell type levels: Low resolution (L1) or high resolution (L2)"),
  make_option(c("-c", "--cell_type"), action="store", default=NULL, type='character',
              help="Cell types in the low cell level resolution (l1) or in the high cell level resolution (l2)"),
  make_option(c("--covariates"), action="store", default='pseudobulkDEA_dreamlet.covariates.tab', type='character',
              help="Covariates file."),
  make_option(c("--vp_reduced"), action="store", default=FALSE, type='character',
              help="Not estimate the effect of the batch factor in the VariancePartition analysis."),
  make_option(c("-i", "--in_dir"), action="store", default='inputs', type='character',
              help="Input directory: WG2 output directory."),
  make_option(c("-o","--out_dir"), action="store", default="pseudobulkDEA_dreamlet", type='character',
              help="Output directory.")
)
opt = parse_args(OptionParser(option_list=option_list))

################################# Set Variables ################################
 # Set working directory
cwd <- getwd()
setwd(cwd)

# Loading functions
functions.fn <- '/tools/wg3-sc_pseudobulk_DEA/scripts/pseudobulkDEA_dreamlet_functions.R'
print(paste0('Loading functions from: ', functions.fn))
source(functions.fn)

# Input/Output directories
in.dir <- paste0(opt$in_dir, '/', opt$cell_level, '/')
out.dir <- paste0(opt$out_dir, '/', opt$cell_level, '/', opt$cell_type, '/')
if(opt$vp_reduced){out.dir <- paste0(out.dir, 'vp_reduced/')}
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

# Covariates 
covs.fn <- opt$covariates

# Seurat object
so.fn <- paste0(in.dir, opt$cell_type, '.Qced.Normalized.SCs.Rds')

# Print report
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Cell type: ', opt$cell_type))
print(paste0('Covariates file: ', opt$covariates))
print(paste0('Input directory: ', opt$in_dir))
print(paste0('Input file: ', so.fn))
print(paste0('Output directory: ', out.dir))
cat('\n')

################################ Read input data ###############################
# Read covariates file
print(paste0('Reading covariates file in: ',covs.fn))
covs.df <- read.table(covs.fn, header = T)
phenotypes <- covs.df[covs.df$type=='fixed',]$covariate

# Read Seurat object file
print(paste0('Reading PBMC seurat object file in: ',so.fn))
system.time(pbmc_i <- readRDS(so.fn))
DefaultAssay(pbmc_i) <- "RNA"

# Modify metadata
print('Modifying metadata...')
pbmc <- modify_metadata_sc(pbmc_i, phenotypes, covs.df)

# Check covariates
print('Checking covariates...')
cov_vars <- covs.df$covariate
metadata_vars <- colnames(pbmc@meta.data)
if(!all(cov_vars%in%metadata_vars)){
  err_message <- paste0('Your metadata variables in: ', opt$covariates, ' should be in the seurat object metadata: ', so.fn)
  stop(err_message)
}

# Get contrasts according to the phenotypes
contrast_coefName.list <- sapply(phenotypes, function(i) get_coefName(i, pbmc), simplify = FALSE)

# Convert Seurat object to SCE object
print('Convert Seurat to SCE object...')
sce <- as.SingleCellExperiment(pbmc)
cat('\n')

################################ Aggregate to pseudobulk ################################
## Dreamlet, like muscat, performs analysis at the pseudobulk-level by summing raw counts across cells for a given sample and cell type. 
## aggregateToPseudoBulk is substantially faster for large on-disk datasets than muscat::aggregateData.
print('Performing pseudobulk...')
colData(sce)$cell_type <- opt$cell_type
colData(sce)$cell_type <- as.factor(colData(sce)$cell_type)
colData(sce)$Donor_Pool <- as.factor(colData(sce)$Donor_Pool)
system.time(pb <- aggregateToPseudoBulk(sce,
                                        assay = "counts",     
                                        cluster_id = "cell_type", 
                                        sample_id = "Donor_Pool",
                                        verbose = FALSE))
pb_raw <- pb
cat('\n')

########################### DEA and VariancePartition with dreamlet #############################
print('Running dreamlet...')
system.time(dreamlet.res <- dreamlet.func(ge_dge = pb,
                                          covariates = covs.df, 
                                          contrast_list = contrast_coefName.list, 
                                          vp_reduced = opt$vp_reduced,
                                          out_dir = out.dir))
