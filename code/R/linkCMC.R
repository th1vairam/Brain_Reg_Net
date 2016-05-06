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

# Utility functions
makeCopy <- function(synId, parentId, FName = NULL, executed = list(), activityName = NULL){
  OLD_OBJ = synGet(synId)
  NEW_OBJ = File(OLD_OBJ@filePath, name = FName, parentId = parentId)
  NEW_OBJ = synStore(NEW_OBJ, used = OLD_OBJ, executed = executed, activityName = activityName)
}

ActivityName <- 'Make copy of data'

thisFileName <- 'linkCMC.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/Brain_Reg_Net", 
                    ref="branch", 
                    refName='CMC')

thisFile <- getPermlink(repository = thisRepo,
                        repositoryPath=paste0('code/R/', thisFileName))


# Make copies of CMC data
## CMC - DLPFC
EXP = makeCopy('syn3493960', 
               'syn5570049', 
               FName = 'Adjusted Voom Normalised Weighted Residuals', 
               executed = thisFile,
               activityName = ActivityName)

COV = makeCopy('syn3493927', 
               'syn5570049', 
               FName = 'Covariates', 
               executed = thisFile,
               activityName = ActivityName)

## CMC - ACC
EXP = makeCopy('syn4985422', 
               'syn5570054', 
               FName = 'Adjusted Voom Normalised Weighted Residuals', 
               executed = thisFile,
               activityName = ActivityName)

COV = makeCopy('syn4985413', 
               'syn5570054', 
               FName = 'Covariates', 
               executed = thisFile,
               activityName = ActivityName)

## HBCC - DLPFC (ARRAY)
EXP = makeCopy('syn4941606', 
               'syn5570055', 
               FName = 'Adjusted Voom Normalised Weighted Residuals (Array)', 
               executed = thisFile,
               activityName = ActivityName)

COV = makeCopy('syn4941610', 
               'syn5570055', 
               FName = 'Covariates (Array)', 
               executed = thisFile,
               activityName = ActivityName)