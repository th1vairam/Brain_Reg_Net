## R code to extract and curate GSE67835 expression data from GEO
library(CovariateAnalysis)
library(data.table)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)

library(synapseClient)
library(knitr)
library(githubr)

library(GEOquery)

synapseLogin()

# Get commit from github
thisRepo <- getRepo(repository = "th1vairam/Brain_Reg_Net", ref="branch", refName='celltype')
thisFile <- getPermlink(repository = thisRepo, repositoryPath=paste0('code/R/','getGSE67835.R'))

# Get raw data from geo
geo.gse = getGEO('GSE67835')
pheno.data = rbind(pData(geo.gse[[2]]),
                   pData(geo.gse[[1]])) %>%
  dplyr::select(title, geo_accession, characteristics_ch1, characteristics_ch1.1, 
                characteristics_ch1.2, characteristics_ch1.3, characteristics_ch1.4) %>%
  tidyr::separate(characteristics_ch1, c('tmp1', 'Tissue'), sep = '\\: ') %>%
  tidyr::separate(characteristics_ch1.1, c('tmp2', 'CellType'), sep = '\\: ') %>%
  tidyr::separate(characteristics_ch1.2, c('tmp3', 'Age'), sep = '\\: ') %>%
  tidyr::separate(characteristics_ch1.3, c('tmp4', 'ChipID'), sep = '\\: ') %>%
  tidyr::separate(characteristics_ch1.4, c('tmp5', 'SampleName'), sep = '\\: ') %>%
  dplyr::select(-starts_with('tmp')) %>%
  filter(CellType %in% c('oligodendrocytes', 'astrocytes', 'OPC', 'microglia', 'neurons', 'endothelial'))

# Get supplementary files from GEO
geo.files = getGEOSuppFiles('GSE67835')
ind = grep('_RAW.tar', rownames(geo.files))
untar(rownames(geo.files)[ind])
all.files = dir(pattern = '*.csv.gz')
counts = lapply(all.files, function(x){
  dat = read.table(x, sep = '\t', header = F)
  x = stringr::str_split(x, '\\_')[[1]][1]
  colnames(dat) = c('Gene.Names',x)
  return(dat)
}) %>%
  plyr::join_all(type = 'full') %>%
  dplyr::select(Gene.Names, one_of(as.character(pheno.data$geo_accession)))

# Modify covariates
covariates = pheno.data[,c('Tissue', 'CellType', 'Age', 'ChipID')]
covariates = as.data.frame(lapply(covariates, factor))
rownames(covariates) = pheno.data$geo_accession

# Convert counts to cpm
expr = counts[,-(1)]
rownames(expr) = counts$Gene.Names
expr = voom(expr, normalize.method = 'quantile', plot = T)$E

ind = intersect(rownames(covariates), colnames(expr))
expr = expr[,ind]
covariates = covariates[ind,]

# Write files to synapse
expr %>%
  rownameToFirstColumn('hgnc_symbol') %>%
  write.table(file = 'LogCPM.tsv', sep = '\t', quote = F, row.names = F)
expr.obj = File('LogCPM.tsv', parentId = 'syn8077138', name = 'Expression (Log CPM)')
expr.obj = synStore(expr.obj, 
                    activityName = 'Curate GEO data', 
                    used = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835',
                    executed = thisFile)

counts %>%
  rownameToFirstColumn('hgnc_symbol') %>%
  write.table(file = 'Counts.tsv', sep = '\t', quote = F, row.names = F)
counts.obj = File('Counts.tsv', parentId = 'syn8077138', name = 'Expression (Raw counts)')
counts.obj = synStore(counts.obj, 
                      activityName = 'Curate GEO data', 
                      used = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835',
                      executed = thisFile)

covariates %>%
  rownameToFirstColumn('SampleID') %>%
  write.table(file = 'covariates.tsv', sep = '\t', quote = F, row.names = F)
covariates.obj = File('covariates.tsv', parentId = 'syn8077138', name = 'Covariates')
covariates.obj = synStore(covariates.obj, 
                          activityName = 'Curate GEO data', 
                          used = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835',
                          executed = thisFile)