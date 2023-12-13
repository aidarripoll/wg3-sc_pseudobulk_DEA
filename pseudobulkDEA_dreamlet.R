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
  make_option(c("--min_prop"), action="store", default=0.2, type='double',
              help="In dreamlet::processAssays(), minimum proportion of retained samples with non-zero counts for a gene to be retained."),
  make_option(c("--vp_reduced"), action="store", default=FALSE, type='logical',
              help="Not estimate the effect of the batch factor in the VariancePartition analysis."),
  make_option(c("--combined_data"), action="store", default=FALSE, type='logical',
              help="If the data from more than one study has been previously merged in one object."),
  make_option(c("--genoPC1"), action="store", default=FALSE, type='logical',
              help="If the data contains more than one genetic ancestry."),
  make_option(c("-i", "--in_dir"), action="store", default='input', type='character',
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
functions.fn <- 'scripts/pseudobulkDEA_dreamlet_functions.R'
print(paste0('Loading functions from: ', functions.fn))
source(functions.fn)

common_functions.fn <- 'scripts/common_functions.R'
print(paste0('Loading common functions from: ', common_functions.fn))
source(common_functions.fn)

# Input/Output directories
in.dir <- paste0(opt$in_dir, '/', opt$cell_level, '/')
out.dir <- paste0(opt$out_dir, '/', opt$cell_level, '/', opt$cell_type, '/')
if(opt$vp_reduced){out.dir <- paste0(out.dir, 'vp_reduced/')}
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
print(paste0('Covariates file: ', opt$covariates))
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

# Read Seurat object file
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

################################ Aggregate to pseudobulk ################################
# Get contrasts according to the phenotypes
contrast_coefName.list <- sapply(phenotypes, function(i) get_coefName(i, pbmc), simplify = FALSE)
contrast_coefName.list <- Filter(Negate(is.null), contrast_coefName.list)
                                 
# Convert Seurat object to SCE object
print('Convert Seurat to SCE object...')
sce <- as.SingleCellExperiment(pbmc)
cat('\n')

# Aggregate to pseudobulk: Dreamlet, like muscat, performs analysis at the pseudobulk-level by summing raw counts across cells for a given sample and cell type. 
# aggregateToPseudoBulk is substantially faster for large on-disk datasets than muscat::aggregateData.
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
# Testing
# pb_raw <- pb
# set.seed(123)
# genes <- sample(rownames(pb), 5000)
# pb <- pb[genes,]

print('Running dreamlet...')
system.time(dreamlet.res <- dreamlet.func(ge_dge = pb,
                                          covariates = covs.df, 
                                          min_prop = opt$min_prop,
                                          contrast_list = contrast_coefName.list, 
                                          vp_reduced = opt$vp_reduced,
                                          out_dir = out.dir))
