#### Function to extract gene sets from ENRICHR  source files ####
## It is assumed your working directory is where this file

### Clear R console screen output
cat("\014") 

## Load required libraries
library(CovariateAnalysis)
library(data.table)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)

library(synapseClient)
library(knitr)
library(githubr)

library(parallel)
library(doParallel)
library(foreach)

cl = makeCluster(detectCores() - 2)
registerDoParallel(cl)

synapseLogin()

#### Synapse specific parameters ####
parentId = 'syn5570248';
activityName = 'Extract Genesets';
activityDescription = 'Extract genesets from enrichr';

# Github link
thisRepo <- getRepo(repository = "th1vairam/Brain_Reg_Net", ref="branch", refName='geneSetCuration')
thisFile <- getPermlink(repository = thisRepo, repositoryPath='code/R/extractGeneSets.R')

#### Query Synapse for source files ####
all.files = synQuery('select id,name from file where parentId == "syn4598359"')

GeneSets = plyr::dlply(all.files, .(file.id), .fun = function(y){
  print(paste('Started',y$file.name))
  tmp = read.table(synGet(y$file.id)@filePath, fill = NA, row.names = NULL, header = FALSE)
  tmp1 = plyr::dlply(tmp, .(V1), .fun = function(x){
    as.character(unlist(x[-(1)])) %>% unique() %>% setdiff('')
  }, 
  .parallel = T,
  .paropts = list(.packages = c('dplyr')))
  print(paste('Completed',y$file.name))
  return(tmp1)
 }, .progress = 'text')
names(GeneSets) = gsub('.txt','',all.files$file.name)
save(list = 'GeneSets', file = 'allEnrichrGeneSets.RData')
stopCluster(cl)

#### Store in synapse ####
obj = File('allEnrichrGeneSets.RData', name = 'Gene Sets in RList Format ', parentId = syn4867780)
obj = synStore(obj, used = all.files$file.id, activityName = activityName, 
               executed = thisFile, activityDescription = activityDescription)