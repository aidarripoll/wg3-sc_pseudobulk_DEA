#!/usr/bin/env Rscript

############################### Load R packages ################################
print('Loading R packages...')
shhh <- suppressPackageStartupMessages
shhh(library(MAST)) #MAST version should be >= 1.16.0 --> with R version 4.0.3 (2020-10-10), we got MAST version 1.17.3
shhh(library(Seurat))
shhh(library(SeuratObject))
shhh(library(data.table))
shhh(library(lme4))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(tidyr))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(sparseMatrixStats))

######################## Functions used in scDEA_MAST_glmer.R ##################
# 1. Convert Seurat object to SCE object; and then from SCE object to SCA object
so_to_sca <- function(so){
  sca <- as(as.SingleCellExperiment(so), 'SingleCellAssay')
  return(sca)
}

# 2. Calculate the frequency of expression (i.e., the proportion of non-zero values in sc) --> mymic MAST::freq() but using sparseMatrixStats package which is faster
freq_sparse <- function(sc, na.rm = TRUE){
  stopifnot(is(sc, "SingleCellAssay"))
  sc_mat <- t(assay(sc, withDimnames = FALSE)) > 0
  out <- sparseMatrixStats::colMeans2(sc_mat, rows = NULL, cols = NULL, na.rm = na.rm)
  names(out) <- colnames(sc_mat)
  return(out)
  }
  
# 3. Filter out lowly variable genes --> CDR (=cngeneson; cellular detection rate) --> Filter out lowly expressed genes
preprocess_sca <- function(sca, freq_expr = 0.1){
  # Initial # of genes
  ngenes_inital <- nrow(sca)
  print(paste0('Initial # of genes: ', ngenes_inital))
  
  # Filter out lowly variable genes
  print('Filtering out lowly variable genes...')
  system.time(freq_sca <- freq_sparse(sca)) #MAST::freq() adapted --> faster
  # system.time(freq_sca <- freq(sca)) #MAST::freq() adapted --> too slow (used before)
  
  # table(freq_sca>0)
  system.time(sca <- sca[freq_sca>0,]) #returns the frequency of expression, i.e., the proportion of non-zero values in sc
  freq_sca <- freq_sca[names(freq_sca)%in%rownames(sca)]
  ngenes_variable <- nrow(sca)
  print(paste0('Highly variable # of genes: ', ngenes_variable))
  
  # Modify metadata (colData and rowData)
  colData(sca)$wellKey <- colnames(sca)
  rowData(sca)$primerid <- rownames(rowData(sca))
  
  # Calculating CDR (=cngeneson)
  cdr2 <- colSums(assay(sca, withDimnames = FALSE)>0) #assay(sca) is working on assays(sca)$logcounts
  colData(sca)$cngeneson <- scale(cdr2)
  
  # Filter out lowly expressed genes
  print(paste0('Filtering out lowly expressed genes: select the genes found in at least ', freq_expr, ' of the cells (proportion of non-zero cells)'))
  expressed_genes <- freq_sca > freq_expr
  sca <- sca[expressed_genes,]
  ngenes_expressed <- nrow(sca)
  print(paste0('Highly expressed # of genes: ', ngenes_expressed))
  
  return(sca)
}

# 4. Check model variables and remove terms
drop_terms <- function(x, sca_object){
  print(x)
  sca_md <- as.data.frame(colData(sca_object))
  var.boolean <- TRUE
  if(class(sca_md[[x]])%in%c('numeric', 'integer')){
    print('Checking variance>0...')
    x_var <- stats::var(sca_md[[x]])
    if(x_var==0){var.boolean <- FALSE}
  }else{
    print('Checking nLevels>1...')
    x_nlevels <- length(unique(sca_md[[x]]))
    if(x_nlevels==1){var.boolean <- FALSE}
  }
  print(var.boolean)
  cat('\n')
  return(var.boolean)
}

# 5. scDEA with MAST glmer
de_glmer.func <- function(sca_object, contrast, fixed_effects, random_effects, res, out_dir){
  # sc-DEA with MAST glmer (zlm)
  print(paste0('Testing: ', contrast))
  random_effects.fmla <- paste(paste0('(1|',random_effects,')'),collapse='+') #try other nomenclature (=option 1, default)
  contrast_fixed.fmla <- paste(c(contrast,fixed_effects),collapse='+')
  zlm_vars <- paste0('~',paste(c(contrast_fixed.fmla, random_effects.fmla), collapse='+'))
  zlm_formula <- as.formula(zlm_vars)
  
  contrast_LRT <- contrast
  if(contrast%in%c('SEX','age_cat', 'age_cat_all')){
    contrast_LRT <- paste0(contrast, levels(colData(sca_object)[[contrast]])[2])
  }
  print(paste0('Fitting glmer: ',zlm_vars))
  print(paste0('Trying to mitigate model failing to converge for some genes by passing nAGQ=0 to the fitting function zlm (fitArgsD=list(nAGQ=0))'))
  zlmCond <- zlm(zlm_formula, 
                 sca_object,
                 method='glmer', ebayes=FALSE, fitArgsD=list(nAGQ=0))
  summaryCond <- summary(zlmCond, doLRT=contrast_LRT, fitArgsD=list(nAGQ=0))
  summaryDt <- summaryCond$datatable
  
  # Collect residuals (with glm, not working with glmer)
  if(res){
    print('Collecting residuals with MAST glm...')
    # Collecting residuals
    fixed_effects.fmla <- paste0('~',paste(fixed_effects,collapse='+'))
    print(paste0('Collecting residuals (with glm): ',fixed_effects.fmla))
    fixed_effects.formula <- as.formula(fixed_effects.fmla)
    
    window <- function(x1) lapply(assays(x1), function(x2) x2[, 1:2])
    
    ## discrete residuals
    z1 <- zlm(fixed_effects.formula, sca_object, hook=discrete_residuals_hook)
    z1_residuals <- collectResiduals(z1, sca_object)
    # window(z1_residuals)
    
    ## continuous residuals
    z2 <- zlm(fixed_effects.formula, sca_object, hook=continuous_residuals_hook)
    z2_residuals <- collectResiduals(z2, sca_object)
    # window(z2_residuals)
    
    ## combined residuals
    z3 <- zlm(fixed_effects.formula, sca_object, hook=combined_residuals_hook)
    z3_residuals <- collectResiduals(z3, sca_object)
    # window(z3_residuals)
    
    residualsList <- list(discrete = assays(z1_residuals)$Residuals,
                          continuous = assays(z2_residuals)$Residuals,
                          combined = assays(z3_residuals)$Residuals)
    
    summary_bygene <- lapply(residualsList, function(x) summary(t(x))) #check summary of each of the residuals (discrete, continuous and combined) per gene
    
    # Saving residuals
    residuals.fn <- paste0(out_dir, 'residuals_glm.rds')
    print(paste0('Saving residuals with MAST glm in: ', residuals.fn))
    saveRDS(residualsList, residuals.fn)
  }
  return(summaryDt)
}

# 6. Get sc-DEGs from MAST glmer output
get_degs <- function(summaryDt, phenotype){
  contrast_LRT <- phenotype
  if(phenotype%in%c('SEX','age_cat', 'age_cat_all')){
    contrast_LRT <- grep(phenotype, levels(summaryDt$contrast), value = TRUE)
  }
  fcHurdle <- merge(summaryDt[contrast==contrast_LRT & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast==contrast_LRT & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], #logFC coefficients
                    by = 'primerid')
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdle <- as.data.frame(fcHurdle)
  fcHurdle <- fcHurdle[order(fcHurdle$fdr),]
  fcHurdle$direction <- ifelse(is.na(sign(fcHurdle$coef)), NA, 
                               ifelse(sign(fcHurdle$coef)==1, 'up', 'down'))
  fcHurdle$ss <- ifelse(is.na(sign(fcHurdle$fdr)), NA, 
                        ifelse(fcHurdle$fdr<=0.05, 'ss', 'ns'))
  fcHurdle <- fcHurdle[order(fcHurdle$coef, fcHurdle$fdr),]
  fcHurdleSig <- droplevels(fcHurdle[!is.na(fcHurdle$fdr) & fcHurdle$fdr<=0.05,])
  fcHurdleSig <- fcHurdleSig[order(fcHurdleSig$coef, fcHurdleSig$fdr),]
  
  return(fcHurdleSig)
}
