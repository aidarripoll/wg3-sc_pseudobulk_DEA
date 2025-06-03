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
shhh(library(DescTools))
shhh(library(lme4))
shhh(library(lmerTest))
shhh(library(broom))
shhh(library(broom.mixed))
shhh(library(parallel))
shhh(library(caret))
shhh(library(tidyverse))

######################## Functions used in pseudobulkDEA_limmadream.R ##################
# 1. Filtering: processAssays()
## accessory - Get contrasts
get_coefName <- function(phe, so){
  md <- so@meta.data
  contrast_var <- phe
  if(is.factor(md[[phe]])){
    contrast_var.levels <- levels(md[[phe]])
    if(length(contrast_var.levels)>1){
      contrast_var <- paste0(phe, contrast_var.levels[[2]])
    }else{
      print(paste0(phe, ' only has 1 level: ', contrast_var.levels, '. NOT consider this phenotype.'))
      contrast_var <- NULL
    }
  }
  return(contrast_var)
}

### accessory - Define formula (VP or DEA)
define_form <- function(gt, df, vp){
  # forms
  print(gt)
  fixed_var_dea <- df[df$type=='fixed',]$covariate
  random_var_dea <- df[df$type=='random',]$covariate
  if(vp){
    if(gt=='VP'){
      print('Not considering the batch effect in the VariancePartition...')
      random_var_dea.idx <- which(df$covariate==random_var_dea)
      df <- df[-random_var_dea.idx,]
    }
  }
  model_vars <- df$covariate
  
  if(gt=='VP'){
    df[df$type=='fixed' & df$class=='factor',]$type <- 'random' 
  }
  
  fixed_var <- df[df$type=='fixed',]$covariate
  fixed.fmla <- paste(fixed_var,collapse='+')
  random_var <- df[df$type=='random',]$covariate
  random.fmla <- NULL
  if(length(random_var)>0){
    random.fmla <- paste(paste0('(1|',random_var,')'),collapse='+')
  }
  form_vars <- paste(c(fixed.fmla,random.fmla), collapse='+')
  form_vars <- paste0('~',form_vars)
  print(paste0('Fitting lmer: ',form_vars))
  form <- as.formula(form_vars)
  
  # specificy colors
  sex.hex <- brewer.pal(9, 'Greens')[7]
  age.hex <- brewer.pal(9, 'Blues')[7]
  random.hex <- brewer.pal(9, 'Greys')[7]
  residuals.hex <- brewer.pal(9, 'Greys')[3]
  cols_vars <- c(sex.hex, age.hex, random.hex, residuals.hex)
  names(cols_vars) <- c(fixed_var_dea[c(1,2)], random_var_dea[1], 'Residuals')
  if(length(random_var_dea)>1 | length(fixed_var_dea)>1){
    if(length(random_var_dea)>1){
      random_added.hex <- brewer.pal(9, 'Greys')[5]
      names(random_added.hex) <- random_var_dea[2]
      cols_vars <- c(cols_vars, random_added.hex)
    }
    if(length(fixed_var_dea)>1){
      fixed_added.hex <- brewer.pal(9, 'Oranges')[5]
      names(fixed_added.hex) <- fixed_var_dea[3]
      cols_vars <- c(cols_vars, fixed_added.hex)
    }
  }
  model_vars_in <- c(model_vars, 'Residuals')
  cols_vars <- cols_vars[names(cols_vars)%in%model_vars_in]
  cols_vars <- cols_vars[match(model_vars_in, names(cols_vars))]
  
  # output
  out <- list(form = form,
              cols = cols_vars)
  
  return(out)
}

## main - processAssays()
# ge_dge = pb
# covariates = covs.df
# min_prop = opt$min_prop
# cell_counts = opt$nCells_per_donor
# vp_reduced = opt$vp_reduced
# out_dir = out.dir
# gene_test = c('VP','DEA')
process_data <- function(ge_dge, covariates, min_prop, cell_counts, vp_reduced, out_dir, gene_test = c('VP','DEA')){
  ### Defining the VP/DEA formulas ###
  print('Defining the VP/DEA formulas...')
  gene_test.forms <- sapply(gene_test, function(i) define_form(i, covariates, vp_reduced), simplify = FALSE)
  
  #### Normalize and apply voom/voomWithDreamWeights ####
  # Run processAssays()
  form <- gene_test.forms$DEA$form
  print('Normalizing the pseudobulk-data...')
  system.time(res.proc <- processAssays(ge_dge, form, 
                                        min.cells = 5, 
                                        min.count = 5, 
                                        min.samples = 4, 
                                        min.prop = min_prop))
  
  # View details of dropping samples
  details(res.proc)
  
  # Check nSamples and nGenes tested
  genes_all <- rownames(ge_dge)
  genes_tested <- rownames(as.data.frame(res.proc))
  genes_all.n <- nrow(ge_dge)
  genes_tested.n <- nrow(as.data.frame(res.proc))
  genes_tested.prop <- round(genes_tested.n/genes_all.n,3)
  samples_all <- colnames(ge_dge)
  samples_tested <- colnames(as.data.frame(res.proc))
  samples_all.n <- ncol(ge_dge)
  samples_tested.n <- ncol(as.data.frame(res.proc))
  samples_tested.prop <- round(samples_tested.n/samples_all.n,3)
  print(paste0('# Genes tested: ', genes_tested.n, ', out of ', genes_all.n, ' (', genes_tested.prop, ')'))
  print(paste0('# Samples tested: ', samples_tested.n, ', out of ', samples_all.n, ' (', samples_tested.prop, ')'))
  
  # Show voom plot for each cell clusters
  ## Here the mean-variance trend from voom is shown for each cell type. Cell types with sufficient number of cells and reads show a clear mean-variance trend. While in rare cell types like megakaryocytes, fewer genes have sufficient reads and the trend is less apparent.
  plotVoom.p <- plotVoom(res.proc)
  plotVoom.fn <- paste0(out_dir, 'plotVoom.png')
  ggsave(plotVoom.fn, plotVoom.p)
  
  # Filter pb 
  ## Keep samples & genes after processAssays()
  res_ct.proc <- res.proc[[1]]
  voomWithDreamWeights.mat <- res_ct.proc$E
  samples_kept <- colnames(voomWithDreamWeights.mat)
  genes_kept <- rownames(voomWithDreamWeights.mat)
  
  ## Get formula
  de_form <- res_ct.proc$formula
  
  ## Filter pb
  pb <- ge_dge
  pb.ge <- as.matrix(assays(pb)[[1]])
  pb_filt.ge <- pb.ge[rownames(pb.ge)%in%genes_kept, colnames(pb.ge)%in%samples_kept]
  
  # Get metadata
  pb.md <- as.data.frame(colData(pb))
  pb_filt.md <- pb.md[rownames(pb.md)%in%samples_kept,]
  
  # Add nCells per donor
  if(cell_counts){
    # Get cell counts per donor
    cellCounts_df <- as.data.frame(cellCounts(pb))
    
    # Add cell counts to metadata
    cellCounts_df <- cellCounts_df %>% rownames_to_column("rownames")
    colnames(cellCounts_df) <- c('rownames', 'nCells')
    pb_filt_tmp.md <- pb_filt.md
    pb_filt_tmp.md <- pb_filt_tmp.md %>% rownames_to_column("rownames")
    pb_filt.md <- pb_filt_tmp.md %>% left_join(cellCounts_df, by = "rownames")
    pb_filt.md <- pb_filt.md %>% column_to_rownames("rownames")
    
    # Add cell counts to formula
    de_form.deparse_ncells <- paste0(deparse(de_form), ' + nCells')
    de_form <- as.formula(de_form.deparse_ncells)
  }
  
  #edgeR::cpm() transformation
  log2cpm.mat <- edgeR::cpm(pb_filt.ge, log = TRUE)
  log2cpm.mat <- log2cpm.mat[,match(rownames(pb_filt.md), colnames(log2cpm.mat))]
  identical(colnames(log2cpm.mat), rownames(pb_filt.md)) #check
  
  # variancePartition::voomWithDreamWeights() transformation
  voomWithDreamWeights.mat <- voomWithDreamWeights.mat[,match(rownames(pb_filt.md), colnames(voomWithDreamWeights.mat))]
  identical(colnames(voomWithDreamWeights.mat), rownames(pb_filt.md)) #check

  # Output
  ## check
  dim(pb_filt.ge)
  dim(voomWithDreamWeights.mat)
  dim(log2cpm.mat)
  pb_filt.ge[1:3,1:3]
  voomWithDreamWeights.mat[1:3,1:3]
  log2cpm.mat[1:3,1:3]
  
  ## return
  out <- list(expr_counts = pb_filt.ge, #raw
              expr_voomWithDreamWeights = voomWithDreamWeights.mat, #voomWithDreamWeights
              expr_log2cpm = log2cpm.mat, #log2cpm
              donor_metadata = pb_filt.md,
              de_form = de_form)
  return(out)
}
  
# 2. LMER
## accessory - Inverse-normal rank transformation (INRT)
# x <- gene_expr
rankTransform <- function(x){
  x_norm <- x
  print(paste0('nSamples: ', length(x)))
  notNA <- which(!is.na(x))
  print(paste0('nSamples (not NA): ', length(x[notNA])))
  percentile <- rank(x[notNA], ties.method='random', na.last = NA)/(length(x)+1)
  # percentile <- rank(x[notNA], ties.method='random', na.last = NA)/(length(x[notNA])+1) #check with Maxime
  mean_level <- mean(Winsorize(x[notNA]))
  sd_level <- sd(Winsorize(x[notNA]))
  x[notNA] <- qnorm(percentile, mean_level, sd_level)
  
  # check normal distribution
  ## If the p-value > 0.05: Fail to reject the null hypothesis, meaning the data is likely normal.
  ## If the p-value ≤ 0.05: Reject the null hypothesis, meaning the data is not normal.
  # x_norm_inrt <- x
  # shapiro.test(x_norm)
  # shapiro.test(x_norm_inrt)
  # ks.test(x_norm, "pnorm", mean = mean(x_norm), sd = sd(x_norm))
  # ks.test(x_norm_inrt, "pnorm", mean = mean(x_norm_inrt), sd = sd(x_norm_inrt))
  
  return(x)
}

## main - LMM using lme4::lmer()
# g <- genes_expressed[1]
# norm <- norm_methods[2]
# process_data_list = process_data.list
# phenotype = opt$phenotype
# freqCut_nzv = 95/5
# uniqueCut_nzv = 10
lmer_func <- function(g, norm, process_data_list, phenotype, freqCut_nzv = 95/5, uniqueCut_nzv = 10){
  print(paste0('# Normalize: ', norm))
  print(paste0('# Gene: ', g))
  
  ### prepare data ###
  # get expression data
  ## normalized
  norm_i <- paste0('expr_', norm)
  norm_data <- process_data_list[[norm_i]]
  gene_vec <- norm_data[g,]
  
  ## inverse-normal rank transformation
  donors <- names(gene_vec)
  gene_expr <- unname(gene_vec)
  gene_expr_inrt <- rankTransform(gene_expr)
  names(gene_expr_inrt) <- donors
  
  # get donor metadata
  donor_md <- process_data_list$donor_metadata
  donor_md$assignment <- rownames(donor_md)
  
  # dataframe
  gene_df <- data.frame(assignment=names(gene_expr_inrt), value=unname(gene_expr_inrt))
  df_i <- merge(gene_df, donor_md, by = 'assignment')
  rownames(df_i) <- df_i$assignment
  df_i <- df_i[,-1]
  df_i <- df_i %>% mutate_if(is.character, as.factor)

  # check covariates variance --> nearZeroVar()
  metrics <- 'value'
  nearZeroVar.metrics <- nearZeroVar(df_i[,metrics], freqCut = freqCut_nzv, uniqueCut = uniqueCut_nzv, names=TRUE)
  nearZeroVar.df <- nearZeroVar(df_i[,metrics], freqCut = freqCut_nzv, uniqueCut = uniqueCut_nzv, saveMetrics=TRUE)
  nearZeroVar.df$nObs <- nrow(df_i)
  if(length(nearZeroVar.metrics)>0){
    print(paste0('nearZeroVar: ', nearZeroVar.metrics))
  }

  # formula (testing)
  # covs <- c('Age', 'Gender', '(1|date)')
  # fmla_covs <- paste(covs, collapse = '+')
  
  form_tmp <- process_data_list$de_form
  fmla_covs <- deparse(form_tmp)
  fmla <- paste0('value',fmla_covs)
  form <- as.formula(fmla)
  print(paste0('Fitting lmer: ',fmla))

  # check variance of the fitted value
  if(var(df_i$value)!=0){
    # lmer
    mod <-  lmerTest::lmer(form, data = df_i)
    tidy_mod <- broom.mixed::tidy(mod, conf.int = TRUE, effects = "fixed")
    tidy_mod <- as.data.frame(tidy_mod)
  }else{
    print(paste0('There is no variance of the fitted value. All values equal to: ', as.character(unique(df_i$value))))
    tidy_mod <- data.frame(effect = rep('fixed',3),
                           term = c('(Intercept)', 'Age', 'GenderF'))
    tidy_mod.na <- as.data.frame(matrix(NA, nrow = 3, ncol = 7))
    colnames(tidy_mod.na) <- c('estimate', 'std.error', 'statistic',
                               'df', 'p.value', 'conf.low', 'conf.high')
    tidy_mod <- cbind(tidy_mod, tidy_mod.na)
  }
  out <- list(model = tidy_mod,
              nearZeroVar = nearZeroVar.df)

  cat('\n')
  return(out)
}

