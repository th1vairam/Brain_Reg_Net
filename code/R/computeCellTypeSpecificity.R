## Compute cell type specificity metrics
## NOTE: Metrics were obtained from A benchmark of gene expression tissue-specificity metrics
##       DOI: 10.1093/bib/bbw008

# Load libraries
library(CovariateAnalysis)
library(data.table)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)
library(matrixStats)

library(synapseClient)
library(knitr)
library(githubr)

synapseLogin()

# Get reference expression and cell type labels from Darmanis., S et al., 2015
ALL_USED_IDs = c('syn8077142', 'syn8077191')
ref.expr = downloadFile('syn8077142')
rownames(ref.expr) = ref.expr$hgnc_symbol
ref.expr$hgnc_symbol = NULL

ref.cellTypes = downloadFile('syn8077191') %>%
  dplyr::select(SampleID, CellType)

# Summarise reference expression interms of cell types
ref.expr.sum = list()
for(celltype in unique(ref.cellTypes$CellType)){
  ids = ref.cellTypes$SampleID[ref.cellTypes$CellType == celltype]
  ref.expr.sum[[celltype]] = log2(apply((2^ref.expr[,ids]+0.5), 1, median, na.rm = T)) %>%
    rownameToFirstColumn('hgnc_symbol') %>%
    plyr::rename(c('DF' = celltype))
}
ref.expr.sum = join_all(ref.expr.sum)

# Compute tau
tmp = data.matrix(ref.expr.sum[,-(1)])
rownames(tmp) = ref.expr.sum$hgnc_symbol
tau = rowSums(1 - (tmp/rowMaxs(tmp)), na.rm = T)
tau = (tau/(dim(tmp)[2]-1)) %>%
  rownameToFirstColumn('hgnc_symbol') %>%
  dplyr::rename(tau = DF)

# Compute TSI
TSI = rowMaxs(tmp)/rowSums(tmp) %>%
  as.data.frame()
colnames(TSI) = 'TSI'
TSI$hgnc_symbol = rownames(TSI)

# Compute Hg
pi = tmp/rowSums(tmp)
Hg = -rowSums(pi*log2(pi)) %>%
  rownameToFirstColumn('hgnc_symbol') %>%
  dplyr::rename(Hg = DF)

# Get github commit
thisRepo <- getRepo(repository = "th1vairam/Brain_Reg_Net", 
                    ref="branch", 
                    refName='celltype')
thisFile <- getPermlink(repository = thisRepo, 
                        repositoryPath='code/R/computeCellTypeSpecificity.R')

# Store results to synapse
plyr::join_all(list(tau, TSI, Hg), type = 'full', match = 'all') %>%
  write.table(file = 'CellSpecificityMetrics.tsv', sep = '\t', row.names = F, quote=F)
obj = File('CellSpecificityMetrics.tsv', parentId = 'syn8077138')
obj = synStore(obj, activityName = 'Compute cell type specificity',
               used = ALL_USED_IDs,
               executed = thisFile)