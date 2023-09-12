# sc-eQTLGen WG3 pipeline (II): sc- and pseudobulk-differential expression analysis 

We provide **two main scripts** to perform **differential expression analysis (DEA)** with different conditions, such as human phenotypes (e.g., sex or age) or stimulation conditions, using **single-cell RNA-seq data (scRNA-seq) (i.e., 10x Genomics)** at two different levels:

1. **[single-cell level (sc-DEA)](/scDEA_MAST_glmer.R)**: using the [MAST](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5) glmer implementation.
2. **[pseudobulk level (pseudobulk-DEA)](/pseudobulkDEA_dreamlet.R)**: using the [dreamlet](https://www.biorxiv.org/content/10.1101/2023.03.17.533005v1) glmer implementation.

**Of note**: 
* These analyses are meant to be run on scRNA-seq data composed by **only one sample per donor**. Some considerations:
- If you have biological replicates (e.g., stimulated vs. non-stimulated samples from the same donor, etc...), you will have to modify some configuration files ([scDEA.covariates.tab](/scDEA.covariates.tab), [pseudobulkDEA_dreamlet.covariates.tab](/pseudobulkDEA_dreamlet.covariates.tab)).
- If you have technical replicates, we will only select the sample with the largest number of cells.

* In these analyses, we are only using European individuals and non-stimulated samples, which have been previously selected in WG3 (I) pipeline.
  
* To run these scripts you should have **successfully run** the following sc-eQTLGen consortium pipelines: **WG1**, **WG2** and **WG3 (I)** 

## Contact
If you have any questions or issues, feel free to open an issue or directly email Aida Ripoll-Cladellas (aida.ripoll@bsc.es)

-------

## Required Software
**R** >=4.1.2 version: You need to install the packages loaded in the:
* Two main DEA scripts: [sc-DEA](/scDEA_MAST_glmer.R) and [pseudobulk-DEA](/pseudobulkDEA_dreamlet.R).
* Additional scripts in the [scripts](/scripts/) directory (which contain the functions called in the three main DEA scripts).

-------

## Required Input
This section explains the input data and its structure to run the two main scripts: [sc-DEA](/scDEA_MAST_glmer.R), [pseudobulk-DEA](/pseudobulkDEA_dreamlet.R).

**Of note**: 
* To follow better the explanations in the **Required Input** section, you can clone this repository and change your current working directory. 

```
git clone https://github.com/aidarripoll/wg3-sc_pseudobulk_DEA.git
cd wg3-sc_pseudobulk_DEA
```

### Test Data
We have provided some **testing inputs** in the **[input directory](input/)** that contains the B cells outputs (Azimuth's level 1) from WG3 (I) pipeline. 

**Of note**: 

* These files have been anonymized and they are significantly down-sized and sub-sampled versions of the whole B cells outputs from WG3 (I). The total number of cells is 663 from 40 donors, and the number of genes is 100.

* Here is the structure of the [testing input directory](/input/). This input directory (*input/*) should have the same structure as the WG3 (I) pipeline output directory. We will need only the following files since the other ones will be used for the eQTL calling pipeline in the WG3 (II):

**input/**    
```bash
|-- L1
|   |-- B.Qced.Normalized.SCs.Rds 
|   |-- B.qtlInput.Pcs.txt
|-- smf.txt

````

### Required Data
**wg3-sc_pseudobulk_DEA/**  
```bash
|-- input
|   |-- L1
|   |   |-- B.Qced.Normalized.SCs.Rds 
|   |   |-- B.qtlInput.Pcs.txt
|   |-- smf.txt
|-- scDEA.covariates.tab 
|-- pseudobulkDEA_dreamlet.covariates.tab 
```

#### 1. Main inputs: [input/](/input/)
1. Sampling matching file (**[smf.txt](/input/smf.txt)**). 

2. **[L1](/input/L1/)** (or L2) directory: A directory for each Azimuth's level with the main outputs per cell type from WG3 (I): 
* **${cell_type}.Qced.Normalized.SCs.Rds:** QC-filtered Seurat object
* **${cell_type}.qtlInput.Pcs.txt:** Sample's PCs

#### 2. DEA covariates files: sc-DEA ([scDEA.covariates.tab](/scDEA.covariates.tab)) and pseudobulk-DEA ([pseudobulkDEA_dreamlet.covariates.tab](/pseudobulkDEA_dreamlet.covariates.tab))
A priori, these files should not be modified. Each tsv file has:
* 1st column (covariate): Covariates included in the model
* 2nd column (type): Fixed/random effect
* 3rd column (class): Categorical (factor) or quantitative (integer, double)

*Of note*:
* Tab separated.
* The values of the 1st column (covariate) should be the same as the ones in the metadata from the [QC-filtered Seurat object](/input/L1/B.Qced.Normalized.SCs.Rds)
* This file must have this header. 
* The covariates files provided for testing have the following structure:

**2.1. sc-DEA** ([scDEA.covariates.tab](/scDEA.covariates.tab)):

| covariate  | type | class  | 
| ------------- | ------------- | ------------- | 
| SEX  | fixed  | factor  | 
| age  | fixed  | integer  | 
| Donor  | random  | factor  | 
| Pool  | random  | factor  | 

**2.2. pseudobulk-DEA** ([pseudobulkDEA_dreamlet.covariates.tab](/pseudobulkDEA_dreamlet.covariates.tab)):

| covariate  | type | class  | 
| ------------- | ------------- | ------------- | 
| SEX  | fixed  | factor  | 
| age  | fixed  | integer  | 
| Pool  | random  | factor  |

-------

## Running the sc/pseudobulk-DEA
*Of note*: 

* If you have not done it yet, the first step would be to *clone** this repository and change your current working directory. 

```
git clone https://github.com/aidarripoll/wg3-sc_pseudobulk_DEA.git
cd wg3-sc_pseudobulk_DEA
```

* These analyses are meant to be run on scRNA-seq data composed by **only one sample per donor**:
1. If you have biological replicates (e.g., stimulated vs. non-stimulated samples from the same donor, etc...), you will have to modify some configuration files ([scDEA.covariates.tab](/scDEA.covariates.tab), [pseudobulkDEA_dreamlet.covariates.tab](/pseudobulkDEA_dreamlet.covariates.tab))
2. If you have technical replicates, we will only select the sample with the largest number of cells.

* In these analyses, we are only using European individuals and non-stimulated samples, which have been previously selected in WG3 (I) pipeline.

* The **functions** called in the **sc/pseudobulk-DEA scripts** ([sc-DEA](/scDEA_MAST_glmer.R), [pseudobulk-DEA](/pseudobulkDEA_dreamlet.R) are defined in the [additional scripts](/scripts/).


### Running the sc-DEA script
**1.** Set common environment variables:  
```
cell_level=L1 #default
cell_type=B
phenotype=SEX #or age
input_directory=input #default
output_directory=scDEA_MAST_glmer #default
```

**2.** Running the **[sc-DEA](/scDEA_MAST_glmer.R)**:

```
Rscript scDEA_MAST_glmer.R -l $cell_level -c $cell_type -v $phenotype -i $input_directory -o $output_directory
```

*Of note:*
* If you are using the **WG3 singularity image**, you should set another common environment variable: `covariates_file=/tools/wg3-sc_pseudobulk_DEA/scDEA.covariates.tab` and run **[sc-DEA](/scDEA_MAST_glmer.R)** as:

 ```
Rscript scDEA_MAST_glmer.R -l $cell_level -c $cell_type -v $phenotype --covariates $covariates_file -i $input_directory -o $output_directory
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

**3.** Additional parameters: 
* If your dataset contains data from more than one study that has been previously merged into one [QC-filtered Seurat object](/input/L1/B.Qced.Normalized.SCs.Rds), you should set the argument `--combined_data TRUE` to control for the differences between studies.

* If your dataset contains data from more than one genetic ancestry, you should set the argument `--genoPC1 TRUE` to control for the differences between ancestries.
  
### Running the pseudobulk-DEA script

**1.** Set common environment variables:  
```
cell_level=L1 #default
cell_type=B
input_directory=input #default
output_directory=pseudobulkDEA_dreamlet #default
```

**2.** Running the **[pseudobulk-DEA](/pseudobulkDEA_dreamlet.R)**:

```
Rscript pseudobulkDEA_dreamlet.R -l $cell_level -c $cell_type -i $input_directory -o $output_directory
```

*Of note:*
* If you are using the **WG3 singularity image**, you should set another common environment variable: `covariates_file=/tools/wg3-sc_pseudobulk_DEA/pseudobulkDEA_dreamlet.covariates.tab` and run **[pseudobulk-DEA](/pseudobulkDEA_dreamlet.R)** as:

 ```
Rscript pseudobulkDEA_dreamlet.R -l $cell_level -c $cell_type --covariates $covariates_file -i $input_directory -o $output_directory
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
* You should **run** the [pseudobulk-DEA](/pseudobulkDEA_dreamlet.R) per **each combination** of Azimuth's cell level - cell type:
    - Azimuth's cell level (`-l`): L1 or L2.
    - Cell type (`-c`): L1 and L2 cell types.

* This script is running both a **pseudobulk-DEA** and a **VariancePartition analysis**. However, some of the outputs will only be generated if we find differentially expressed genes (DEGs) (FDR<=0.05) with each condition (plotGeneHeatmap.png, plotPercentBars.png, and plotVarPart.png).

* It runs the **pseudobulk-DEA** on the 'fixed' variables in [pseudobulkDEA_dreamlet.covariates.tab](/pseudobulkDEA_dreamlet.covariates.tab)

**3.** Additional parameters: 
* If your dataset contains data from more than one study that has been previously merged into one [QC-filtered Seurat object](/input/L1/B.Qced.Normalized.SCs.Rds), you should set the argument `--combined_data TRUE` to control for the differences between studies.

* If your dataset contains data from more than one genetic ancestry, you should set the argument `--genoPC1 TRUE` to control for the differences between ancestries.
  
* If you do not want to estimate the variance of your batch variable (Pool), you can set the argument `--vp_reduced TRUE`. The structure of the output directory will be the same as the previous one but it will have a subdirectory called 'vp_reduced'. 

-------

## Running time and memory requirements
* [sc-DEA](/scDEA_MAST_glmer.R), [pseudobulk-DEA](/pseudobulkDEA_dreamlet.R): To speed up the running time and improve the memory requirements of these three main scripts, we recommend to submit each of the commands of the sections *'Running the sc-DEA script'* and *'Running the pseudobulk-DEA script'* as an independent job on your HPC infrastructure (i.e., run each job as an element of a job array). In the case of the [**sc-DEA**](/scDEA_MAST_glmer.R), each job will be defined by the combination of: `Azimuth's cell level - cell type - phenotype`. In the case of the [**pseudobulk-DEA**](/pseudobulkDEA_dreamlet.R), each job will be defined by the combination of: `Azimuth's cell level - cell type`.

* The HPC resources for the major Azimuth's l1 cell type (CD4 T cells) in the V2 unstimulated data from [**Oelen et al, 2022**](https://www.nature.com/articles/s41467-022-30893-), which is composed by 72 individuals and 33,505 cells, using the defined SLURM's `sbatch` parameters (`--nodes=1`, `--ntasks=1`, `--cpus-per-task=48`, `--tasks-per-node=1`) in [MareNostrum 4](https://www.bsc.es/supportkc/docs/MareNostrum4/intro/), where each node has 48 cores of 1.880GB/core, were:
1. [**sc-DEA**](/scDEA_MAST_glmer.R): 
* testing SEX: Elapsed (05:17:54), MaxRSS (19690752K) and CPUTime (10-14:19:12)
* testing age: Elapsed (05:17:42), MaxRSS (19191280K) and CPUTime (10-14:09:36)

2. [**pseudobulk-DEA**](/pseudobulkDEA_dreamlet.R): 
* testing SEX and age (in the same script): Elapsed (00:19:35), MaxRSS (3467004K) and CPUTime (15:40:00)





