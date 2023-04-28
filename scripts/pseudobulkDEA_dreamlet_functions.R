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
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(ggplot2))
shhh(library(RColorBrewer))

######################## Functions used in pseudobulkDEA_limmadream.R ##################
# 1. Modify sc-metadata
modify_metadata_sc <- function(so, phes, covs_df){
  # Grab metadata
  metadata <- so@meta.data
  cells_all <- nrow(metadata)
  print(paste0('Number of initial cells: ', cells_all))
  
  # Remove possible NAs
  for (i in phes){
    na.boolean <- any(is.na(metadata[[i]]))
    if(na.boolean){
      print('Removing NAs in the phenotype...')
      metadata <- metadata[-which(is.na(metadata[[i]])),]
      cells_notna <- nrow(metadata)
      print(paste0('Number of not NA cells: ', cells_notna))
    }
  }
  
  # Create random factors variables in the metadata
  metadata[c('Donor','Pool')] <- str_split_fixed(metadata$Donor_Pool, ';;', 2)
  
  # Phenotype modifications
  if('SEX'%in%phes){
    print(paste0('Relevel SEX...'))
    metadata[['SEX']] <-  ifelse(metadata[['SEX']]==1,'M','F')
    Group_order <- c('M','F')
    metadata[['SEX']] <- factor(metadata[['SEX']],
                                levels = Group_order)
  }
  
  if('age_cat'%in%phes){
    print(paste0('Creating a new metadata variable: age_cat...'))
    metadata[['age_cat']] <- ifelse(metadata$age<=40, 'Y',
                                    ifelse(metadata$age>=60, 'O', 'M'))
    metadata <- metadata[metadata[['age_cat']]%in%c('Y','O'),] #only keep Y and O samples (remove if we want to consider all samples)
    Group_order <- c('Y','M','O')
    metadata[['age_cat']] <- factor(metadata[['age_cat']],
                                    levels = Group_order)
  }
  
  if('age_cat_all'%in%phes){
    print(paste0('Creating a new metadata variable:  age_cat_all...'))
    metadata[['age_cat_all']] <- ifelse(metadata$age<=40, 'Y', 'O')
    Group_order <- c('Y','O')
    metadata[['age_cat_all']] <- factor(metadata[['age_cat_all']],
                                        levels = Group_order)
  }
  
  if('age_squared'%in%phes){
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
  rownames(so@meta.data) <- so@meta.data$Barcode
  
  return(so)
}

# 2. Get contrasts
get_coefName <- function(phe, so){
  md <- so@meta.data
  contrast_var <- phe
  if(is.factor(md[[phe]])){
    contrast_var.levels <- levels(md[[phe]])
    contrast_var <- paste0(phe, contrast_var.levels[[2]])
  }
  return(contrast_var)
}

# 3. Define formula (VP or DEA)
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
  random.fmla <- paste(paste0('(1|',random_var,')'),collapse='+')
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
  names(cols_vars) <- c(fixed_var_dea, random_var_dea, 'Residuals')
  model_vars_in <- c(model_vars, 'Residuals')
  cols_vars <- cols_vars[names(cols_vars)%in%model_vars_in]
  
  # output
  out <- list(form = form,
              cols = cols_vars)
  
  return(out)
}

# 4. DEA extract and plots
extract_plots <- function(i, dea_res, contrast_var, vp_res, cols, o_dir){
  # Extract results (topTable --> DEGs)
  ### Each entry in res.dl stores a model fit by dream(), and results can be extracted using topTable() as in limma by specifying the coefficient of interest. 
  ### The results shows the gene name, log fold change, average expression, t-statistic, p-value, FDR (i.e. adj.P.Val).
  genes <- rownames(dea_res[[1]]$residuals)
  topTable.res <- topTable(dea_res, 
                           coef = contrast_var,
                           number = length(genes))
  degs <- topTable.res[topTable.res$adj.P.Val<=0.05,]$ID
  
  ### DEA plots (for all genes) ###
  # Volcano plots
  ## The volcano plot can indicate the strength of the differential expression signal with each cell type. Red points indicate FDR < 0.05.
  plotVolcano.p <- plotVolcano(dea_res, coef = contrast_var)
  plotVolcano.fn <- paste0(o_dir, 'plotVolcano.png')
  print(paste0('Saving plotVolcano in: ', plotVolcano.fn))
  ggsave(plotVolcano.fn, plotVolcano.p)
  
  ### DEA and VP plots (only if DEGs) ###
  if(length(degs)>0){
    ### DEA plots ###
    # Gene-level heatmap
    ## For each cell type and specified gene, show z-statistic from dreamlet analysis. 
    ## Grey indicates that insufficient reads were observed to include the gene in the analysis.
    plotGeneHeatmap.p <- plotGeneHeatmap(dea_res, coef=contrast_var, genes=degs)
    plotGeneHeatmap.fn <- paste0(o_dir, 'plotGeneHeatmap.png')
    print(paste0('Saving plotGeneHeatmap in: ', plotGeneHeatmap.fn))
    ggsave(plotGeneHeatmap.fn, plotGeneHeatmap.p)
    
    ## Forest plot
    ## A forest plot shows the log fold change and standard error of a given gene across all cell types. The color indicates the FDR.
    # os_dir <- paste0(o_dir, '/plotForest/')
    # if(!dir.exists(os_dir)){dir.create(os_dir, recursive = T)}
    # plotForest.save <- lapply(degs, function(i){
    #   plotForest.p <- plotForest(dea_res, coef = contrast_var, gene = i)
    #   plotForest.fn <- paste0(os_dir, i, '.png')
    #   print(paste0('Saving plotForest in: ', plotForest.fn))
    #   ggsave(plotForest.fn, plotForest.p)
    #   return(NULL)
    # })
    
    ### VP plots ###
    # Pick only DEGs in the VP results
    vp_res.degs <- vp_res[vp_res$gene%in%degs,]
    cnames <- c('assay','gene', names(cols))
    vp_res.degs <- vp_res.degs[,match(cnames, colnames(vp_res.degs))]
    vp_res.degs <- vp_res.degs[match(degs, vp_res.degs$gene),]
    
    # plotPercentBars --> some genes show differences when estimating only Gender/Age or Gender/Age/date
    plotPercentBars.p <- plotPercentBars(vp_res.degs, cols)
    plotPercentBars.fn <- paste0(o_dir, 'plotPercentBars.png')
    print(paste0('Saving plotPercentBars in: ', plotPercentBars.fn))
    ggsave(plotPercentBars.fn, plotPercentBars.p)
    
    # plotVarPart
    plotVarPart.p <- plotVarPart(vp_res.degs, cols, label.angle=60) 
    plotVarPart.fn <- paste0(o_dir, 'plotVarPart.png')
    print(paste0('Saving plotVarPart in: ', plotVarPart.fn))
    ggsave(plotVarPart.fn, plotVarPart.p)
  }
  cat('\n')
  res <- topTable.res
  return(res)
} 

# 5. DEA + VP extract and plots (by phenotype)
extract_plots_by_phe <- function(phe, dea_res, vp_res, c_list, cols, o_dir){
  print(phe)
  
  # create output dir
  out_sdir <- paste0(o_dir, '/', phe, '/')
  if(!dir.exists(out_sdir)){dir.create(out_sdir, recursive = T)}
  
  # pick contrast variable
  contrast_coefName <- c_list[[phe]]
  
  # extract and plots
  res <- extract_plots(i = phe,
                       dea_res = dea_res,
                       contrast_var = contrast_coefName,
                       vp_res = vp_res,
                       cols = cols,
                       o_dir = out_sdir)
  return(res)
}

# 6. dreamer
dreamlet.func <- function(ge_dge, covariates, contrast_list, vp_reduced, gene_test = c('VP','DEA'), out_dir = out.dir){
  ### Defining the VP/DEA formulas ###
  print('Defining the VP/DEA formulas...')
  gene_test.forms <- sapply(gene_test, function(i) define_form(i, covariates, vp_reduced), simplify = FALSE)
  
  #### Normalize and apply voom/voomWithDreamWeights ####
  # Run processAssays()
  form <- gene_test.forms$DEA$form
  print('Normalizing the pseudobulk-data...')
  system.time(res.proc <- processAssays(ge_dge, form, min.count=5))
  
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
  
  ### Differential expression ###
  ## Since the normalized expression data and metadata are stored within res.proc, only the regression formula remains to be specified.
  ## Here we only included the stimulus status, but analyses of larger datasets can include covariates and random effects.
  ## With formula ~ StimStatus, an intercept is fit and coefficient StimStatusstim log fold change between simulated and controls.
  ## Differential expression analysis within each assay, evaluated on the voom normalized data
  print('Running DEA...')
  system.time(res.dl <- dreamlet(res.proc, form))
  
  ### Variance partitioning ###
  ## The variancePartition package uses linear and linear mixed models to quanify the contribution of multiple sources of expression variation at the gene-level.
  ## For each gene it fits a linear (mixed) model and evalutes the fraction of expression variation explained by each variable.
  ## Variance fractions can be visualized at the gene-level for each cell type using a bar plot, or genome-wide using a violin plot.
  # Mymic DEA model: https://github.com/GabrielHoffman/dreamlet/issues/4#issuecomment-1507767030
  # vp.lst = fitVarPart(res.proc, form) # not working --> 2 alternatives (use option 1):
  #### 1. Use ~(1|Sex)+Age+(1|Batch). variancePartition works best when categorical variables are modeled as a random effects. It't not an issue that this formula isn't identical the the differential expression formula.
  #### 2. We can try to regress out the Batch variable, and do fitVarPart() on the residuals --> You could do that, but you'd have to use variancePartition::fitExtractVarPartModel() directly.
  print('Running VariancePartition...')
  system.time(vp.lst <- fitVarPart(res.proc, gene_test.forms$VP$form))
  cols_vars <- gene_test.forms$VP$cols
  
  ### Extract results and plots (DEA and VP) ###
  extract_plots_by_phe.res <- sapply(names(contrast_list),
                                     function(i) extract_plots_by_phe(phe = i,
                                                                      dea_res = res.dl,
                                                                      vp_res = vp.lst,
                                                                      c_list = contrast_list,
                                                                      cols = cols_vars,
                                                                      o_dir = out_dir), simplify = FALSE)
  
  ### Save outputs ###
  res <- list(dea = res.dl,
              vp = vp.lst,
              topTable = extract_plots_by_phe.res)
  res_fn <- paste0(out_dir, 'dea_vp_topTable.rds')
  print(paste0('Saving dreamlet results: ', res_fn))
  saveRDS(res, res_fn)
  
  return(res)
}


