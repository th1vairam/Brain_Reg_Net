## It is assumed your working directory is where this file is

# Clear R console screen output
cat("\014")  

# Load required libraries
library(CovariateAnalysis)
library(data.table)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)

library(synapseClient)
library(knitr)
library(githubr)

library(biomaRt)
library(ComplexHeatmap)

synapseLogin()

# Get the latest commit of the code from git
thisFileName <- 'collateDiffExp.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/Brain_Reg_Net", ref="branch", refName='AMPAD')
thisFile <- getPermlink(repository = thisRepo, repositoryPath=paste0('code/R/',thisFileName))

# Get protein coding genes from biomaRt
ensembl=useMart("ensembl")
ensemblHSapiens = useDataset("hsapiens_gene_ensembl",mart=ensembl)
humanProteinCodingGenes = getBM(attributes=c("ensembl_gene_id","external_gene_name","description"), 
                                filters='biotype', 
                                values=c('protein_coding'), 
                                mart=ensemblHSapiens)

#### Get all covariates
covariates.id = c('ROSMAP' = 'syn6114444', 
                  'MSBB' = 'syn5573110',
                  'Mayo' = 'syn6088521')
covariates = lapply(covariates.id, downloadFile)

covariates$ROSMAP = covariates$ROSMAP %>%
  dplyr::select(Sampleid, msex, cogdx) %>%
  dplyr::rename(ID = Sampleid, Gender = msex, Status = cogdx) %>%
  dplyr::mutate(Gender = factor(Gender, labels = c('1' = 'male','2'='female')),
                Status = factor(Status, labels = c('1' = 'NCI', '2' = 'MCI', '3' = 'MCI.CI', 
                                                   '4' = 'AD', '5' = 'AD.CI', '6' = 'OD')),
                BrainRegion = 'DLPFC') %>%
  dplyr::filter(Status %in% c('NCI', 'AD'))

covariates$MSBB = covariates$MSBB %>%
  dplyr::select(SampleId, SEX, BrainRegion.Dx) %>%
  tidyr::separate(BrainRegion.Dx, c('BrainRegion', 'Dx'), sep = '\\.') %>%
  dplyr::rename(ID = SampleId, Gender = SEX, Status = Dx) %>%
  dplyr::mutate(Gender = factor(Gender, labels = c('M' = 'male','F'='female')),
                Status = factor(Status)) %>%
  dplyr::filter(Status %in% c('ND', 'SD'))

covariates = rbindlist(covariate, use.names = T, fill = T, idcol = 'Study')

#### Get raw expression values
logcpm.id = c('ROSMAP' = 'syn6114447', 'MSBB' = 'syn6117295')
logcpm = lapply(logcpm.id, downloadFile)

# Get all differential expresison results from synapse
diffexp.id = c('ROSMAP' = 'syn6114453',
               'MSBB' = 'syn5609009',
               'Mayo' = 'syn6088525')

diffexp = lapply(diffexp.id, downloadFile) %>%
  rbindlist(idcol = 'DataSetName', fill=T, use.names=T) %>%
  filter(ensembl_gene_id %in% humanProteinCodingGenes$ensembl_gene_id)

# Plot ROSMAP
filter(diffexp, adj.P.Val <= 0.05, Comparison == '')
  

# # Store differential expression results in synapse
# write.table(diffExp, file = 'differentialExpression.tsv', sep = '\t', row.names=F, quote=F)
# DEXP_OBJ = File('differentialExpression.tsv', 
#                 name = 'Differential expression results of all the RNASeq data (protein coding only)', 
#                 parentId = 'syn5569102')
# DEXP_OBJ = synStore(DEXP_OBJ, 
#                     used = as.character(diffExp.ids), 
#                     activityName = 'Collate all the differential expression analysis results', 
#                    executed = thisFile)
# 
# logFC = diffExp %>%
#   filter(adj.P.Val <= 0.05, abs(logFC) <= 1) %>%
#   dplyr::select(DataSetName, Comparison, ensembl_gene_id, logFC) %>%
#   unite(DataSet.Comparison, DataSetName, Comparison, sep = ' ') %>%
#   spread(DataSet.Comparison, logFC)