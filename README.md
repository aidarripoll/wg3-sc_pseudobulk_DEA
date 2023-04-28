# sc-eQTLGen WG3 pipeline (II): sc- and pseudobulk-differential expression analysis 

We provide **three main scripts** to peform **differential expression analysis (DEA)** with different conditions, such as human phenotypes (e.g., sex or age) or stimulation conditions, using **single-cell RNA-seq data (scRNA-seq) (i.e., 10x Genomics)** at two different levels:

1. **[single-cell level (sc-DEA)](/scDEA_MAST_glmer.R)**: using the [MAST](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5) glmer implementation.
2. **[pseudobulk level (pseudobulk-DEA updated)](/pseudobulkDEA_dreamlet.R)**: using the [dreamlet](https://www.biorxiv.org/content/10.1101/2023.03.17.533005v1) glmer implementation.
3. **[pseudobulk level (pseudobulk-DEA outdated)](/pseudobulkDEA_limmadream.R)**: using the [limma dream](https://academic.oup.com/bioinformatics/article/37/2/192/5878955) glmer implementation. This script was replaced by the **[previous script (2)](/pseudobulkDEA_dreamlet.R)** using the [dreamlet](https://www.biorxiv.org/content/10.1101/2023.03.17.533005v1) glmer implementation. 

**Of note**: 
* You will only need to run the two following scripts: [**sc-DEA**](/scDEA_MAST_glmer.R) (script 1) and [**pseudobulk-DEA updated**](/pseudobulkDEA_dreamlet.R) (script 2). 

* These analyses are meant to be run on scRNA-seq data composed by **only one sample per donor**. In case you have more than one sample per donor (e.g., stimulated vs. non-stimulated samples from the same donor, etc...) you should modify some configuration files ([scDEA_covariates.tab](/scDEA_covariates.tab), [pseudobulkDEA_dreamlet.covariates.tab](/pseudobulkDEA_dreamlet.covariates.tab), [pseudobulkDEA_limmadream.covariates.tab](/pseudobulkDEA_limmadream.covariates.tab), and [pseudobulkDEA_limmadream.aggregates.tab](/pseudobulkDEA_limmadream.aggregates.tab))

* To run these scripts you should have **successfully run** the following sc-eQTLGen consortium pipelines: **WG1**, **WG2** and **WG3 (I)** 

## Contact
If you have any questions or issues, feel free to open an issue or directly email Aida Ripoll-Cladellas (aida.ripoll@bsc.es)

-------

## Required Software
**R** >=4.1.2 version: You need to install the packages loaded in the:
* Three main DEA scripts: [sc-DEA](/scDEA_MAST_glmer.R), [pseudobulk-DEA updated](/pseudobulkDEA_dreamlet.R), and [pseudobulk-DEA outdated](/pseudobulkDEA_limmadream.R).
* Additional scripts in the [scripts](/scripts/) directory (which contain the functions called in the three main DEA scripts).

-------

## Required Input
This section explains the input data and it’s structure to run the three main scripts: [sc-DEA](/scDEA_MAST_glmer.R), [pseudobulk-DEA updated](/pseudobulkDEA_dreamlet.R), and [pseudobulk-DEA outdated](/pseudobulkDEA_limmadream.R). 

**Of note**: 
* Remember that you only need to run the two following scripts: [**sc-DEA**](/scDEA_MAST_glmer.R) and [**pseudobulk-DEA updated**](/pseudobulkDEA_dreamlet.R).
* To follow better the explanations in the **Required Input** section, you can clone this repository and change your current working directory. 

```
git clone https://github.com/aidarripoll/wg3-sc_pseudobulk_DEA.git
cd wg3-sc_pseudobulk_DEA
```

### Test Data
We have provided some **testing inputs** in the **[inputs directory](inputs/)** that contains the B cells outputs (Azimuth's level 1) from WG3 (I) pipeline. 

**Of note**: 

* These files have been anonymized and they are significantly down-sized and sub-sampled versions of the whole B cells outputs from WG3 (I). The total number of cells is 663 from 40 donors, and the number of genes is 50 and 500 for the sc- and pseudobulk-DEA, respectively.

* Here is the structure of the [testing input directory](/inputs/). This input directory (*inputs/*) should have the same structure as the WG3 (I) pipeline output directory. We will need only the following files since the other ones will be used for the eQTL calling pipeline in the WG3 (II):

**inputs/**    
```bash
|-- L1
|   |-- B.Exp.txt
|   |-- B.Qced.Normalized.SCs.Rds
|   |-- B.covariates.txt 
|-- donor_pool_stim.txt

````

### Required Data
**wg3-sc_pseudobulk_DEA/**  
```bash
|-- inputs
|   |-- L1
|   |   |-- B.Exp.txt
|   |   |-- B.Qced.Normalized.SCs.Rds
|   |   |-- B.covariates.txt 
|   |-- donor_pool_stim.txt
|-- scDEA_covariates.tab 
|-- pseudobulkDEA_dreamlet.covariates.tab 
|-- pseudobulkDEA_limmadream.covariates.tab 
|-- pseudobulkDEA_limmadream.aggregates.tab
```

#### 1. Main inputs: [inputs/](/inputs/)
It contains a directory for each Azimuth's level (**[L1](/inputs/L1/)** or L2) with the main outputs per cell type from WG3 (I): 
* **${cell_type}.Exp.txt:** Pseudobulk gene expression matrix per donor-pool combination (PFlogPF normalization + mean on QC-filtered single-cell gene expression matrix)
* **${cell_type}.Qced.Normalized.SCs.Rds:** QC-filtered single-cell gene expression matrix
* **${cell_type}.covariates.txt:** Sample metadata (from the psam file in WG1 pipeline)

#### 2. DEA covariates files: sc-DEA ([scDEA_covariates.tab](/scDEA_covariates.tab)), pseudobulk-DEA updated ([pseudobulkDEA_dreamlet.covariates.tab](/pseudobulkDEA_dreamlet.covariates.tab)), and pseudobulk-DEA outdated ([pseudobulkDEA_limmadream.covariates.tab](/pseudobulkDEA_limmadream.covariates.tab))
A priori, these files should not be modified. Each tsv file that has in the:
* 1st column (covariate): Covariates included in the model
* 2nd column (type): Fixed/random effect
* 3rd column (class): Categorical (factor) or quantitative (integer, double)

*Of note*:
* Tab separated.
* The values of the 1st column (covariate) should be the same as the ones in the [pseudobulk-](/inputs/L1/B.Exp.txt) and [sc-](/inputs/L1/Qced.Normalized.SCs.Rds)gene expression files, and in the [sample metadata file](/inputs/L1/B.covariates.txt).  
* This file must have this header. 
* The covariates files provided for testing have the following structure:

**2.1. sc-DEA** ([scDEA_covariates.tab](/scDEA_covariates.tab)):

| covariate  | type | class  | 
| ------------- | ------------- | ------------- | 
| SEX  | fixed  | factor  | 
| age  | fixed  | integer  | 
| Donor_Pool  | random  | factor  | 
| Pool  | random  | factor  | 

**2.2. pseudobulk-DEA updated** ([pseudobulkDEA_dreamlet.covariates.tab](/pseudobulkDEA_dreamlet.covariates.tab)):

| covariate  | type | class  | 
| ------------- | ------------- | ------------- | 
| SEX  | fixed  | factor  | 
| age  | fixed  | integer  | 
| Pool  | random  | factor  |

**2.3. pseudobulk-DEA outdated** ([pseudobulkDEA_limmadream.covariates.tab](/pseudobulkDEA_limmadream.covariates.tab)):

| covariate  | type | class  | 
| ------------- | ------------- | ------------- | 
| SEX  | fixed  | factor  | 
| age  | fixed  | integer  | 
| CellCount  | fixed  | integer  | 
| Pool  | random  | factor  |

#### 3. Pseudobulk aggregation file (only used in the [pseudobulk-DEA outdated](/pseudobulkDEA_limmadream.R)): [pseudobulkDEA_limmadream.aggregates.tab](/pseudobulkDEA_limmadream.aggregates.tab)
A priori, this filee should not be modified. Each tsv file that has in the:
* 1st column (simplified): Simplified aggregation variable.
* 2nd column (complete): Complete aggregation variable.

*Of note*:
* Tab separated.
* All the values should appear as column names of the [sample metadata file](/inputs/L1/B.covariates.txt). 
* The values of the 2nd column (complete) correspond to the column names of the [pseudobulk-gene expressio nfile](/inputs/L1/B.Exp.txt).
* This file must have this header. 
* The pseudobulk aggregation file files provided for testing have the following structure:

| simplified  | complete |       
| ------------- | ------------- |     
| Donor  | Donor;Pool  |     


-------

## Running the sc/pseudobulk-DEA
*Of note*: 

* If you have not done it yet, the first step would be to *clone** this repository and change your current working directory. 

```
git clone https://github.com/aidarripoll/wg3-sc_pseudobulk_DEA.git
cd wg3-sc_pseudobulk_DEA
```

* These analyses are meant to be run on scRNA-seq data composed by **only one sample per donor**. In case you have more than one sample per donor (e.g., stimulated vs. non-stimulated samples from the same donor, etc...) you should modify some configuration files ([scDEA_covariates.tab](/scDEA_covariates.tab), [pseudobulkDEA_dreamlet.covariates.tab](/pseudobulkDEA_dreamlet.covariates.tab), [pseudobulkDEA_limmadream.covariates.tab](/pseudobulkDEA_limmadream.covariates.tab), and [pseudobulkDEA_limmadream.aggregates.tab](/pseudobulkDEA_limmadream.aggregates.tab)). For this analysis, we are only using european individuals and non-stimulated samples, which has been previously selected in WG3 (I) pipeline.

* The **functions** called in the **sc/pseudobulk-DEA scripts** ([sc-DEA](/scDEA_MAST_glmer.R), [pseudobulk-DEA updated](/pseudobulkDEA_dreamlet.R).[pseudobulk-DEA outdated](/pseudobulkDEA_limmadream.R)) are defined in the [additional scripts](/scripts/).


### Running the sc-DEA script

**1.** Set common environmental variables:  
```
cell_level=L1 #default
cell_type=B
phenotype=SEX #or age
input_directory=inputs #default
output_directory=scDEA_MAST_glmer #default
```

**2.** Running the **[sc-DEA](/scDEA_MAST_glmer.R)**:

```
Rscript scDEA_MAST_glmer.R -l $cell_level -c $cell_type -v $phenotype -i $input_directory -o $output_directory
```

The output directory (**[scDEA_MAST_glmer/](/scDEA_MAST_glmer/))** has the following structure:
```bash
L1
|-- SEX
|   |-- B
|       |-- de_glmer_nagq0.degs.rds
|       |-- de_glmer_nagq0.rds
|       |-- pbmc_sca.rds
|       |-- pbmc_sca_raw.rds
|       |-- pbmc_so.rds
|-- age
    |-- B
        |-- de_glmer_nagq0.degs.rds
        |-- de_glmer_nagq0.rds
        |-- pbmc_sca.rds
        |-- pbmc_sca_raw.rds
        |-- pbmc_so.rds

```

*Of note:*
* You should **run** the [sc-DEA](/scDEA_MAST_glmer.R) per **each combination** of Azimuth's cell level - cell type - phenotype:
    - Azimuth's cell level (`-l`): L1 or L2.
    - Cell type (`-c`): L1 and L2 cell types.
    - Phenotype (`-v`): SEX and age.

### Running the pseudobulk-DEA script (updated)

**1.** Set common environmental variables:  
```
cell_level=L1 #default
cell_type=B
input_directory=inputs #default
output_directory=pseudobulkDEA_dreamlet #default
```

**2.** Running the **[pseudobulk-DEA updated](/pseudobulkDEA_dreamlet.R)**:

```
Rscript pseudobulkDEA_dreamlet.R -l $cell_level -c $cell_type -i $input_directory -o $output_directory
```

The output directory (**[pseudobulkDEA_dreamlet/](/pseudobulkDEA_dreamlet/)**) has the following structure:
```bash
L1
|   |-- B
|        |-- SEX
|            |-- plotVolcano.png
|            |-- plotGeneHeatmap.png
|            |-- plotPercentBars.png
|            |-- plotVarPart.png
|        |-- age
|            |-- plotVolcano.png
|            |-- plotGeneHeatmap.png
|            |-- plotPercentBars.png
|            |-- plotVarPart.png
|    |-- plotVoom.png
|    |-- dea_vp_topTable.rds

```

*Of note:*
* You should **run** the [pseudobulk-DEA updated](/pseudobulkDEA_dreamlet.R) per **each combination** of Azimuth's cell level - cell type:
    - Azimuth's cell level (`-l`): L1 or L2.
    - Cell type (`-c`): L1 and L2 cell types.

* This script is running both a **pseudobulk-DEA** and a **VariancePartition analysis**. However, some of the outputs will only be generated if we find differentially expressed genes (DEGs) (FDR<=0.05) with each condition (plotGeneHeatmap.png, plotPercentBars.png, and plotVarPart.png).

* It runs the **pseudobulk-DEA** on the 'fixed' variables in [pseudobulkDEA_dreamlet.covariates.tab](/pseudobulkDEA_dreamlet.covariates.tab)

**3.** Additional parameters of [pseudobulk-DEA updated](/pseudobulkDEA_dreamlet.R):

3.1. If you do not want to estimate the variance of your batch variable (Pool), you can set the `--vp_reduced TRUE`. The structure of the output directory will be the same as the previous one but it will have a subdirectory called 'vp_reduced'. 

### Running the pseudobulk-DEA script (outdated)

**1.** Set common environmental variables:  
```
cell_level=L1 #default
cell_type=B
phenotype=SEX #or age
input_directory=inputs #default
output_directory=pseudobulkDEA_limmadream #default
```

**2.** Running the **[pseudobulk-DEA outdated](/pseudobulkDEA_limmadream.R)**:

```
Rscript pseudobulkDEA_limmadream.R -l $cell_level -c $cell_type -v $phenotype -i $input_directory -o $output_directory
```

The output directory (**[pseudobulkDEA_limmadream/](/pseudobulkDEA_limmadream/)**) has the following structure:
```bash
L1
|-- SEX
|   |-- B
|       |-- eBayes
|           |-- SEX.combinations.degs.rds
|           |-- SEX.combinations.rds
|           |-- SEXM_SEXF.rds
|           |-- SEXM_SEXF.tsv
|           |-- SEXM_SEXF.vars.rds
|-- age
    |-- B
        |-- eBayes
            |-- age.combinations.degs.rds
            |-- age.combinations.rds
            |-- age.rds
            |-- age.tsv
            |-- age.vars.rds
```

*Of note:*
* You should **run** the [pseudobulk-DEA outdated](/pseudobulkDEA_limmadream.R) per **each combination** of Azimuth's cell level - cell type - phenotype:
    - Azimuth's cell level (`-l`): L1 or L2.
    - Cell type (`-c`): L1 and L2 cell types.
    - Phenotype (`-v`): SEX and age.

**3.** Additional parameters of [pseudobulk-DEA outdated](/pseudobulkDEA_limmadream.R):

3.1. It can be that `variancePartition::dream()` fails due to genes with high weights. It has been observed in some specific cases, specially when using sum to aggregate single-cell to pseudobulk. These are some potential solutions proposed by the `variancePartition` developer (See the [github issue](https://github.com/GabrielHoffman/variancePartition/issues/66) for further details):

* `span`: If we want to keep all the genes, we can reduce the `span` parameter (default = 0.5) to make sure we will not have this issue anymore. This change will not have an important impact on other genes as it is only affecting the right tail (where there are very few points). You can try to set it to 0.1.

* `weights`: If we want to remove the genes with high weights to avoid errors, we can detect these few cases by picking the outliers (high weights) from the weights distribution (`variancePartition::voomWithDreamWeights()` output), remove them from the pseudobulk-gene expression matrix, and run again `variancePartition::dream()`. You can try to set it to 0.001 or 0.005. 

3.2. At this moment, we are not setting any gene expression threshold. We could change this by changing the following parameters:

* `expr`: Select genes with a minimum of frequency across donors. You can try to set it to 0.5.
* `cv`: Select genes with a minimum of variation across donors. You can try to set it to 0.5, which will be equivalent to the top two quartiles based on their squared coefficient of variation (CV^2 = variance / mean^2) calculated across all cells of each different cell-type.

-------

## Running time and memory requirements
* [sc-DEA](/scDEA_MAST_glmer.R), [pseudobulk-DEA updated](/pseudobulkDEA_dreamlet.R), and [pseudobulk-DEA outdated](/pseudobulkDEA_limmadream.R): To speed up the running time and improve the memory requirements of these three main scripts, we recommend to submit each of the commands of the sections *'Running the sc-DEA script'* and *'Running the pseudobulk-DEA script updated'* as an independent job on your HPC infrastructure (i.e., run each job as an element of a job array). In the case of the [**sc-DEA**](/scDEA_MAST_glmer.R), each job will be defined by the combination of: `Azimuth's cell level - cell type - phenotype`. In the case of the [**pseudobulk-DEA updated**](/pseudobulkDEA_dreamlet.R), each job will be defined by the combination of: `Azimuth's cell level - cell type`.

* The HPC resources in for the major Azimuth's l1 cell type (CD4 T cells) in the V2 unstimulated data from [**Oelen et al, 2022**](https://www.nature.com/articles/s41467-022-30893-), which is composed by 72 individuals and 33,505 cells, using the defined SLURM's `sbatch` parameters (`--nodes=1`, `--ntasks=1`, `--cpus-per-task=48`, `--tasks-per-node=1`) in [MareNostrum 4](https://www.bsc.es/supportkc/docs/MareNostrum4/intro/), where each node has 48 cores of 1.880GB/core, were:
1. [**sc-DEA**](/scDEA_MAST_glmer.R): 
* testing SEX: Elapsed (05:17:54), MaxRSS (19690752K) and CPUTime (10-14:19:12)
* testing age: Elapsed (05:17:42), MaxRSS (19191280K) and CPUTime (10-14:09:36)

2. [**pseudobulk-DEA updated**](/pseudobulkDEA_dreamlet.R): 
* testing SEX and age (in the same script): 





