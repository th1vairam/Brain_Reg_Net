---
title: "Covariate analysis of MSBB microarray data"
author: "Thanneer Perumal"
date: "Mon Jan  4 14:13:06 2016"
output: html_document
---

```r
library(knit2synapse)
library(synapseClient)

synapseLogin()

knit2synapse::knitToFolderEntity(file = "./MSSM_Array.Rmd",
                                 parentId ="syn5570248",
                                 entityName = 'MSBB Microarray')
```




Following brain regions were analysed for differential expression between case and control

| Brain Region                        | Synapse ID |
| ----------------------------------- | ---------- |
| middle temporal                     | syn3191101 |
| superior temporal                   | syn3191103 |
| inferior temporal                   | syn3191099 |
| precuneus superior parietal lobule  | syn3191119 |
| caudal anterior cingulate           | syn3191107 |
| frontal pole                        | syn3191095 |
| para hippocampal gyrus              | syn3191109 |


### Extract expression data and adjust covariates 
Working on Middle Temporal

Running PCA and calculating correlations for:
Scaled expression  data in PCA; PVE >= 1%; pearson correlations 
${image?fileName=covariate%2Eanalysis%2D1%2Epng&align=none&scale=100}
Significant covariates are: 
Fitting linear model with the following coefficients: DxAD,DxControl

Running PCA and calculating correlations for:
Scaled residual expression  data in PCA; PVE >= 1%; pearson correlations 
Working on Superior Temporal

Running PCA and calculating correlations for:
Scaled expression  data in PCA; PVE >= 1%; pearson correlations 
${image?fileName=covariate%2Eanalysis%2D2%2Epng&align=none&scale=100}
Significant covariates are: 
Fitting linear model with the following coefficients: DxAD,DxControl

Running PCA and calculating correlations for:
Scaled residual expression  data in PCA; PVE >= 1%; pearson correlations 
Working on Inferior Temporal

Running PCA and calculating correlations for:
Scaled expression  data in PCA; PVE >= 1%; pearson correlations 
${image?fileName=covariate%2Eanalysis%2D3%2Epng&align=none&scale=100}
Significant covariates are: NTrSum
Fitting linear model with the following coefficients: DxAD,DxControl

Running PCA and calculating correlations for:
Scaled residual expression  data in PCA; PVE >= 1%; pearson correlations 
Following coefficients have to be included in the model: NTrSum
Working on Superior Parietal Lobule

Running PCA and calculating correlations for:
Scaled expression  data in PCA; PVE >= 1%; pearson correlations 
Significant covariates are: 
Fitting linear model with the following coefficients: DxAD,DxControl

Running PCA and calculating correlations for:
Scaled residual expression  data in PCA; PVE >= 1%; pearson correlations 
Working on Caudal Anterior Cingulate

Running PCA and calculating correlations for:
Scaled expression  data in PCA; PVE >= 1%; pearson correlations 
Significant covariates are: 
Fitting linear model with the following coefficients: DxAD,DxControl

Running PCA and calculating correlations for:
Scaled residual expression  data in PCA; PVE >= 1%; pearson correlations 
Working on Parahippocampal Gyrus

Running PCA and calculating correlations for:
Scaled expression  data in PCA; PVE >= 1%; pearson correlations 
Significant covariates are: 
Fitting linear model with the following coefficients: DxAD,DxControl

Running PCA and calculating correlations for:
Scaled residual expression  data in PCA; PVE >= 1%; pearson correlations 
Working on Frontal Pole

Running PCA and calculating correlations for:
Scaled expression  data in PCA; PVE >= 1%; pearson correlations 
Significant covariates are: 
Fitting linear model with the following coefficients: DxAD,DxControl

Running PCA and calculating correlations for:
Scaled residual expression  data in PCA; PVE >= 1%; pearson correlations 
### Store files in synapse

