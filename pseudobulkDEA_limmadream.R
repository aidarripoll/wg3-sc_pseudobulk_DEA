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
  make_option(c("--expr"), action="store", default=NULL, type='double',
              help="Expression threshold (50% expression)."),
  make_option(c("--cv"), action="store", default=NULL, type='double',
              help="Top two quartiles based on their squared coefficient of variation (CV^2 = variance / mean^2) calculated across all cells of each different cell-type."),
  make_option(c("--span"), action="store", default=0.5, type='double',
              help="If dream() fails due to genes with high weights (in some cases, when using sum to aggregate single-cell to pseudobulk). If we want to keep all the genes, we can reduce the span parameter to make sure we won't have this issue anymore. This change won't have an important impact on other genes as it is only affecting the right tail (where there are very few points)."),
  make_option(c("--weights"), action="store", default=NULL, type='double',
              help="If dream() fails due to genes with high weights (in some cases, when using sum to aggregate single-cell to pseudobulk). If we want to remove the cases with high weights to avoid errors when using dream(), we can detect these few cases by picking the outliers (high weights) from the weights distribution (voomWithDreamWeights() output), remove them from the gene expression matrix, and run dream()"),
  make_option(c("--covariates"), action="store", default='pseudobulkDEA_covariates.tab', type='character',
              help="Covariates file."),
  make_option(c("--aggregates"), action="store", default='pseudobulkDEA_aggregates.tab', type='character',
              help="Aggregates file."),
  make_option(c("--eBayes"), action="store", default=TRUE, type='logical',
              help="Run eBayes on dream() fit."),
  make_option(c("-i", "--in_dir"), action="store", default='inputs', type='character',
              help="Input directory: WG2 output directory."),
  make_option(c("-o","--out_dir"), action="store", default="pseudobulkDEA_limmadream", type='character',
              help="Output directory.")
)
opt = parse_args(OptionParser(option_list=option_list))

################################# Set Variables ################################
 # Set working directory
cwd <- getwd()
setwd(cwd)

# Loading functions
functions.fn <- 'scripts/pseudobulkDEA_functions.R'
print(paste0('Loading functions from: ', functions.fn))
source(functions.fn)

# Input/Output directories
in.dir <- paste0(opt$in_dir, '/', opt$cell_level, '/')
out.dir <- paste0(opt$out_dir, '/', opt$cell_level, '/', opt$phenotype, '/', opt$cell_type, '/')
if(opt$eBayes){out.dir <- paste0(out.dir, 'eBayes/')}
if(opt$span!=0.5){out.dir <- paste0(out.dir, 'span_', as.character(opt$span), '/')}
if(!is.null(opt$weights)){out.dir <- paste0(out.dir, 'weights_', as.character(opt$weights), '/')}
if(!is.null(opt$expr)){out.dir <- paste0(out.dir, 'expr_', as.character(opt$expr), '/')}
if(!is.null(opt$cv)){out.dir <- paste0(out.dir, 'cv_', as.character(opt$cv), '/')}
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

# Covariates
covs.fn <- opt$covariates

# Aggregates
aggregates.fn <- opt$aggregates

# Metadata
metadata.fn <- paste0(in.dir, opt$cell_type, '.covariates.txt')

# Gene expression
mean_mat.fn <- paste0(in.dir, '/', opt$cell_type, '.Exp.txt')

# Print report
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Cell type: ', opt$cell_type))
print(paste0('Span: ', opt$span))
print(paste0('Weights: ', opt$weights))
print(paste0('Expression threshold: ', opt$expr))
print(paste0('Squared CV threshold: ', opt$cv))
print(paste0('Phenotype: ', opt$phenotype))
print(paste0('Covariates file: ', covs.fn))
print(paste0('Input directory: ', opt$in_dir))
print(paste0('Input file: ', mean_mat.fn))
print(paste0('Output directory: ', out.dir))
cat('\n')

################################ Read input data ###############################
# Read covariates file
print(paste0('Reading covariates file in: ',covs.fn))
covs.df <- read.table(covs.fn, header = T)

# Read aggregates file
print(paste0('Reading aggregates file in: ',covs.fn))
aggregates.df <- read.table(aggregates.fn, header = T)

# Read metadata file
print(paste0('Reading metadata file in: ',metadata.fn))
metadata_df <- read.delim(metadata.fn, check.names = FALSE)

# Read gene expression matrix (mean+PFlogPF)
print(paste0('Reading gene expression (mean+PFlogPF) file in: ',mean_mat.fn))
system.time(mean_pseudo.mat <- read.delim(mean_mat.fn, check.names = FALSE))

# Check phenotype
print('Checking phenotype...')
contrast_vars <- covs.df[covs.df$type=='fixed',]$covariate
contrast_vars.add <- c(contrast_vars, c('age_cat', 'age_cat_all', 'age_squared'))
if(!opt$phenotype%in%contrast_vars.add){
  err_message <- paste0('You should provide one of the following values in --phenotype argument: ', paste(contrast_vars.add, collapse = ', '))
}
cat('\n')

################################ Pre-processing ################################
# Modify metadata
print('Modifying metadata...')
system.time(metadata_list <- modify_metadata_pseudo(metadata_df, opt$phenotype, covs.df, aggregates.df))
metadata <- metadata_list$metadata
aggregates <- metadata_list$aggregates

# Check covariates
print('Checking covariates...')
cov_vars <- covs.df$covariate
metadata_vars <- colnames(metadata)
if(!all(cov_vars%in%metadata_vars)){
  err_message <- paste0('Your metadata variables in: ', opt$covariates, ' should be in the metadata: ', metadata.fn)
  stop(err_message)
}

# Modify gene expression
print('Modifying gene expression matrix...')
system.time(mean_pseudo <- modify_ge(mean_pseudo.mat, aggregates))

# Subset gene expression matrix (if needed --> if opt$phenotype==age_cat, keep only Y and O samples)
ndonors <- ncol(mean_pseudo)
print(paste0('Number of inital donors: ', ndonors))
if(opt$phenotype=='age_cat'){
  print('Subsetting gene expression matrix: keeping only Y and O samples...')
  mean_pseudo <- subset_ge(mean_pseudo, metadata)
}
ndonors_final <- ncol(mean_pseudo)
print(paste0('Number of final donors: ', ndonors_final))

mean_pseudo.expr <- mean_pseudo
if(!is.null(opt$expr)){
  # Pick highly expressed genes 
  print(paste0('Filter gene expression matrix by expression frequency across donors (', opt$expr, ')'))
  system.time(mean_pseudo.expr <- highly_expressed(mean_pseudo, opt$expr))
}

mean_pseudo.expr.hvgs <- mean_pseudo.expr
if(!is.null(opt$cv)){
  # Pick highly variable genes
  print(paste0('Filter gene expression matrix by highly variable genes (top ', as.character(opt$cv), ' quartiles based on their CV-squared)'))
  system.time(mean_pseudo.expr.hvgs <- highly_variable(mean_pseudo.expr, opt$cv))
}
geneExpr <- mean_pseudo.expr.hvgs

# Checking donor ids order (metadata and gene exression matrix)
cat('\n')
print('Checking if the donor ids from the metadata and gene expression matrix are in the same order...')
system.time(aggregate_metadata <- check_donor_ids(geneExpr, metadata))
  
############################## Limma dream #####################################
# Set Fixed and Random effect variables
print('Defining the fixed and random effect variables...')

## Fixed effects
fixed_effects.vars <- covs.df[covs.df$type=='fixed',]$covariate
contrast.var_out <- opt$phenotype
if(opt$phenotype%in%c('age_cat', 'age_cat_all', 'age_squared')){
  contrast.var_out <- 'age'
}
fixed_effects.vars <- fixed_effects.vars[!fixed_effects.vars%in%contrast.var_out]

## Random effects
random_effects.vars <- covs.df[covs.df$type=='random',]$covariate

# Run limma dream
print('Running limma dream()...')
system.time(dreamer.res <- dreamer(ge_dge = geneExpr,
                                   contrast = opt$phenotype,
                                   fixed_effects = fixed_effects.vars,
                                   random_effects = random_effects.vars,
                                   metadata = aggregate_metadata,
                                   eBayes_var = opt$eBayes,
                                   span_var = opt$span,
                                   weights = opt$weights,
                                   out_dir = out.dir))

################################## Get DEGs ####################################
# Get pseudobulk-DEGs from limma dream output
print('Getting the pseudobulk-DEGs...')
degs_list <- sapply(names(dreamer.res), function(i) get_degs(i, dreamer.res), simplify = FALSE)
print(degs_list)

# Saving pseudobulk-DEGs
get_degs.fn <- paste0(out.dir, opt$phenotype, '.combinations.degs.rds')
print(paste0('Saving pseudobulk-DEGs results with limma dream in: ', get_degs.fn))
saveRDS(degs_list, get_degs.fn)