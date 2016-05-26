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
thisRepo <- getRepo(repository = "th1vairam/Brain_Reg_Net", ref="branch", refName='metaAnal')
thisFile <- getPermlink(repository = thisRepo, repositoryPath=paste0('code/R/',thisFileName))

# Get protein coding genes from biomaRt
# ensembl=useMart(biomart = "ensembl", host="www.ensembl.org", dataset = "hsapiens_gene_ensembl")
ensembl=useMart(biomart = "ENSEMBL_MART_ENSEMBL", host="dec2015.archive.ensembl.org", dataset = "hsapiens_gene_ensembl") # Use old version if ensembl is down
humanProteinCodingGenes = getBM(attributes=c("ensembl_gene_id","external_gene_name","description"), 
                                filters='biotype', 
                                values=c('protein_coding'), 
                                mart=ensembl)

#### Get all covariates
covariates.id = c('ROSMAP' = 'syn6114444', 
                  'MSBB' = 'syn5573110',
                  'MAYO' = 'syn6088521')
covariates = lapply(covariates.id, downloadFile)

covariates$ROSMAP = covariates$ROSMAP %>%
  dplyr::select(Sampleid, msex, cogdx) %>%
  dplyr::rename(ID = Sampleid, Gender = msex, Status = cogdx) %>%
  dplyr::mutate(Gender = factor(Gender, labels = c('1' = 'MALE','2'='FEMALE')),
                Status = factor(Status, labels = c('1' = 'NCI', '2' = 'MCI', '3' = 'MCI.CI', 
                                                   '4' = 'AD', '5' = 'AD.CI', '6' = 'OD')),
                BrainRegion = 'DLPFC') %>%
  dplyr::filter(Status %in% c('NCI', 'AD'))

covariates$MSBB = covariates$MSBB %>%
  dplyr::select(SampleId, SEX, BrainRegion.Dx) %>%
  tidyr::separate(BrainRegion.Dx, c('BrainRegion', 'Dx'), sep = '\\.') %>%
  dplyr::rename(ID = SampleId, Gender = SEX, Status = Dx) %>%
  dplyr::mutate(Gender = factor(Gender, labels = c('M' = 'MALE','F'='FEMALE')),
                Status = factor(Status),
                BrainRegion = factor(BrainRegion, labels = c('BM_10' = 'FP', 'BM_22' = 'STG', 
                                                             'BM_36' = 'PHG', 'BM_44' = 'IFG'))) %>%
  dplyr::filter(Status %in% c('ND', 'SD'))

covariates$MAYO = covariates$MAYO %>%
  dplyr::select(ID, Gender, BrainRegion.Diagnosis) %>%
  tidyr::separate(BrainRegion.Diagnosis, c('BrainRegion', 'Dx'), sep = '\\.') %>%
  dplyr::rename(Status = Dx) %>%
  dplyr::mutate(Gender = factor(Gender, labels = c('M' = 'MALE','F'='FEMALE')),
                Status = factor(Status)) %>%
  dplyr::filter(Status %in% c('Control', 'AD'))

covariates = rbindlist(covariates, use.names = T, fill = T, idcol = 'Study')

#### Get raw expression values
logcpm.id = c('ROSMAP' = 'syn6114447', 'MSBB' = 'syn6117295', 'MAYO' = 'syn6121652')
logcpm = lapply(logcpm.id, downloadFile) %>%
  lapply(function(x){
    x = dplyr::filter(x, ensembl_gene_id %in% humanProteinCodingGenes$ensembl_gene_id) %>%
      dplyr::select(-hgnc_symbol)
  }) %>%
  join_all(type = 'full')

# Get all differential expresison results from synapse
diffexp.id = c('ROSMAP' = 'syn6114453',
               'MSBB' = 'syn5609009',
               'Mayo' = 'syn6088525')

diffexp = lapply(diffexp.id, downloadFile) %>%
  rbindlist(idcol = 'Study', fill=T, use.names=T) %>%
  filter(ensembl_gene_id %in% humanProteinCodingGenes$ensembl_gene_id,
         adj.P.Val <= 0.05, Study == 'ROSMAP')
  
# Filter logcpm based on differential expression
counts = filter(logcpm, ensembl_gene_id %in% diffexp$ensembl_gene_id)
rownames(counts) = counts$ensembl_gene_id
counts$ensembl_gene_id = NULL
counts[is.na(counts)] = 0

covariates = data.frame(covariates)
rownames(covariates) = covariates$ID

# Filter and arrange covariates and cpm ocunts
counts = counts[,intersect(rownames(covariates), colnames(counts))]
covariates = covariates[intersect(rownames(covariates), colnames(counts)),]

# Arrange covariates
covariates = covariates %>%
  arrange(Study, BrainRegion, Status, Gender) %>%
  dplyr::select(ID, Study, BrainRegion, Status, Gender) %>%
  data.frame
rownames(covariates) = covariates$ID

# Scale counts
counts = counts[,covariates$ID]
counts = scale(counts)
counts = t(scale(t(counts)))

ha = HeatmapAnnotation(covariates[,-(1)])
Heatmap(counts, top_annotation = ha, name = '', 
        show_row_names = F, show_column_names = F, 
        cluster_columns = F, show_column_dend = F)

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