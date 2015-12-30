# Script to migrate CMC normalised data files for network inference
## It is assumed your working directory is where this file is

### Clear R console screen output
cat("\014")  

# Load required libraries
library(synapseClient)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(plyr)
library(dplyr)
library(psych)
library(tidyr)

library(limma)

library(knitr)
library(rGithubClient) ## Needs the dev branch

synapseLogin()

# Source utility files from ../R/lib folder
file.sources = list.files('../R/lib',pattern="*.R", full.names = T)
tmp = sapply(file.sources,source,.GlobalEnv)

makeCopy <- function(synId, parentId, FName = NULL, executed = list()){
  OLD_OBJ = synGet(synId)
  NEW_OBJ = File(OLD_OBJ@filePath, name = FName, parentId = parentId)
  NEW_OBJ = synStore(NEW_OBJ, used = OLD_OBJ, executed = executed)
}

# Make copies of CMC data
## CMC - DLPFC
EXP = makeCopy('syn3493960', 
               'syn5570049', 
               FName = 'Adjusted Voom Normalised Weighted Residuals', 
               executed = )

COV_OBJ = synGet('syn3493927')
COV = fread(COV_OBJ@filePath, data.table=F, header=T)

## CMC - ACC

## HBCC - DLPFC

## HBCC - ACC

## HBCC - DLPFC (ARRAY)