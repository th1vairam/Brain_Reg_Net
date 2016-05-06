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

# Get all differential expresison results from synapse
diffExp.ids = c('ROSMAP DLPFC' = 'syn5608845',
                'MSBB Frontal Pole' = 'syn6037382',
                'MSBB Superior Temporal Gyrus'	= 'syn6037356',
                'MSBB Para Hippocampal Gyrus' = 'syn6037276',
                'MSBB Inferior Frontal Gyrus' = 'syn6037095',
                'Mayo Cerebellum' = 'syn5609850',
                'Mayo Temporal Cortex' = 'syn5609813')
diffExp = lapply(diffExp.ids, downloadFile) %>%
  rbindlist(idcol = 'DataSetName', fill=T, use.names=T) %>%
  dplyr::select(-genes, -position) %>%
  filter(ensembl_gene_id %in% humanProteinCodingGenes$ensembl_gene_id)

# Store differential expression results in synapse
write.table(diffExp, file = 'differentialExpression.tsv', sep = '\t', row.names=F, quote=F)
DEXP_OBJ = File('differentialExpression.tsv', 
                name = 'Differential expression results of all the RNASeq data (protein coding only)', 
                parentId = 'syn5569102')
DEXP_OBJ = synStore(DEXP_OBJ, 
                    used = as.character(diffExp.ids), 
                    activityName = 'Collate all the differential expression analysis results', 
                   executed = thisFile)

logFC = diffExp %>%
  filter(adj.P.Val <= 0.05, abs(logFC) <= 1) %>%
  dplyr::select(DataSetName, Comparison, ensembl_gene_id, logFC) %>%
  unite(DataSet.Comparison, DataSetName, Comparison, sep = ' ') %>%
  spread(DataSet.Comparison, logFC)

