library(synapseClient)
library(knitr)
library(githubr)

library(CovariateAnalysis)
library(data.table)
library(tidyr)
library(plyr)
library(dplyr)
library(stringr)

synapseLogin()

diff.exp = lapply(c(MAYO = 'syn8468023',
                    MSSM = 'syn10157628',
                    ROSMAP = 'syn8456721'),
                  downloadFile)

############################################
# AD-CONTROL
gs = list()
gs[['DLPFC.AD-CONTROL.UP']] = diff.exp$ROSMAP %>%
  dplyr::filter(Model == 'Diagnosis', Comparison == 'AD-CONTROL', adj.P.Val <= 0.05, logFC >= 0) %>%
  dplyr::select(hgnc_symbol) %>%
  unlist %>% unique()

gs[['DLPFC.AD-CONTROL.DOWN']] = diff.exp$ROSMAP %>%
  dplyr::filter(Model == 'Diagnosis', Comparison == 'AD-CONTROL', adj.P.Val <= 0.05, logFC <= 0) %>%
  dplyr::select(hgnc_symbol) %>%
  unlist %>% unique()

tmp = diff.exp$MSSM %>%
  dplyr::filter(Model == 'Diagnosis', 
                BrainRegion %in% c('IFG', 'PHG' , 'FP', 'STG'), 
                Comparison == 'CONTROL-AD', adj.P.Val <= 0.05, logFC <= 0) %>%
  plyr::dlply(.(BrainRegion), .fun = function(x){
    dplyr::select(x, hgnc_symbol) %>%
      unlist %>% unique()
  })
names(tmp) = paste(names(tmp), 'AD-CONTROL.UP', sep = '.')
gs = c(gs, tmp)

tmp = diff.exp$MSSM %>%
  dplyr::filter(Model == 'Diagnosis', 
                BrainRegion %in% c('IFG', 'PHG' , 'FP', 'STG'), 
                Comparison == 'CONTROL-AD', adj.P.Val <= 0.05, logFC >= 0) %>%
  plyr::dlply(.(BrainRegion), .fun = function(x){
    dplyr::select(x, hgnc_symbol) %>%
      unlist %>% unique()
  })
names(tmp) = paste(names(tmp), 'AD-CONTROL.DOWN', sep = '.')
gs = c(gs, tmp)

tmp = diff.exp$MAYO %>%
  dplyr::filter(BrainRegion.ref %in% c('CER', 'TCX'),
                Model == 'Diagnosis',
                Comparison == 'Control-AD',
                adj.P.Val <= 0.05, logFC >= 0) %>%
  plyr::dlply(.(BrainRegion.ref), .fun = function(x){
    dplyr::select(x, hgnc_symbol) %>%
      unlist %>% unique()
  })
names(tmp) = paste(names(tmp), 'AD-CONTROL.DOWN', sep = '.')
gs = c(gs, tmp)

tmp = diff.exp$MAYO %>%
  dplyr::filter(BrainRegion.ref %in% c('CER', 'TCX'),
                Model == 'Diagnosis',
                Comparison == 'Control-AD',
                adj.P.Val <= 0.05, logFC <= 0) %>%
  plyr::dlply(.(BrainRegion.ref), .fun = function(x){
    dplyr::select(x, hgnc_symbol) %>%
      unlist %>% unique()
  })
names(tmp) = paste(names(tmp), 'AD-CONTROL.UP', sep = '.')
gs = c(gs, tmp)
############################################

############################################
# AD - CONTROL (Male vs Female)
gs1 = list()
tmp = diff.exp$ROSMAP %>%
  dplyr::filter(Model == 'Diagnosis.Gender', Comparison %in% c("AD-CONTROL.IN.FEMALE", "AD-CONTROL.IN.MALE"), 
                adj.P.Val <= 0.05, logFC >= 0) %>%
  plyr::dlply(.(Comparison), .fun = function(x){
    select(x, hgnc_symbol) %>%
      unlist %>% unique()
  })
names(tmp) = c('DLPFC.AD-CONTROL.FEMALE.UP', 'DLPFC.AD-CONTROL.MALE.UP')
gs1 = c(gs1, tmp)

tmp = diff.exp$ROSMAP %>%
  dplyr::filter(Model == 'Diagnosis.Gender', Comparison %in% c("AD-CONTROL.IN.FEMALE", "AD-CONTROL.IN.MALE"), 
                adj.P.Val <= 0.05, logFC <= 0) %>%
  plyr::dlply(.(Comparison), .fun = function(x){
    select(x, hgnc_symbol) %>%
      unlist %>% unique()
  })
names(tmp) = c('DLPFC.AD-CONTROL.FEMALE.DOWN', 'DLPFC.AD-CONTROL.MALE.DOWN')
gs1 = c(gs1, tmp)

tmp = diff.exp$MSSM %>%
  dplyr::filter(Model == 'Diagnosis.SEX', 
                Comparison %in% c('FP.AD.FEMALE-FP.CONTROL.FEMALE','FP.AD.MALE-FP.CONTROL.MALE',
                                  'IFG.AD.FEMALE-IFG.CONTROL.FEMALE','IFG.AD.MALE-IFG.CONTROL.MALE',
                                  'PHG.AD.FEMALE-PHG.CONTROL.FEMALE','PHG.AD.MALE-PHG.CONTROL.MALE',
                                  'STG.AD.FEMALE-STG.CONTROL.FEMALE','STG.AD.MALE-STG.CONTROL.MALE'), 
                adj.P.Val <= 0.05, logFC >= 0) %>%
  plyr::dlply(.(Comparison), .fun = function(x){
    dplyr::select(x, hgnc_symbol) %>%
      unlist %>% unique()
  })
names(tmp) = c("FP.AD-CONTROL.FEMALE.UP", "FP.AD-CONTROL.MALE.UP", 
               "IFG.AD-CONTROL.MALE.UP", 
               "PHG.AD-CONTROL.FEMALE.UP", "PHG.AD-CONTROL.MALE.UP",
               "STG.AD-CONTROL.FEMALE.UP", "STG.AD-CONTROL.MALE.UP")
gs1 = c(gs1, tmp)

tmp = diff.exp$MSSM %>%
  dplyr::filter(Model == 'Diagnosis.SEX', 
                Comparison %in% c('FP.AD.FEMALE-FP.CONTROL.FEMALE','FP.AD.MALE-FP.CONTROL.MALE',
                                  'IFG.AD.FEMALE-IFG.CONTROL.FEMALE','IFG.AD.MALE-IFG.CONTROL.MALE',
                                  'PHG.AD.FEMALE-PHG.CONTROL.FEMALE','PHG.AD.MALE-PHG.CONTROL.MALE',
                                  'STG.AD.FEMALE-STG.CONTROL.FEMALE','STG.AD.MALE-STG.CONTROL.MALE'), 
                adj.P.Val <= 0.05, logFC <= 0) %>%
  plyr::dlply(.(Comparison), .fun = function(x){
    dplyr::select(x, hgnc_symbol) %>%
      unlist %>% unique()
  })
names(tmp) = c("FP.AD-CONTROL.FEMALE.DOWN", "FP.AD-CONTROL.MALE.DOWN", 
               "IFG.AD-CONTROL.FEMALE.DOWN", "IFG.AD-CONTROL.MALE.DOWN", 
               "PHG.AD-CONTROL.FEMALE.DOWN", "PHG.AD-CONTROL.MALE.DOWN",
               "STG.AD-CONTROL.FEMALE.DOWN", "STG.AD-CONTROL.MALE.DOWN")
gs1 = c(gs1, tmp)

tmp = diff.exp$MAYO %>%
  tidyr::unite(Gender, Gender.ref, Gender.ag, sep = '_') %>%
  dplyr::filter(BrainRegion.ref %in% c('CER', 'TCX'), 
                Model == 'Diagnosis.Gender',
                Comparison == 'AD-Control', 
                Gender %in% c('M_M', 'F_F'),
                adj.P.Val <= 0.05, logFC >= 0) %>%
  plyr::dlply(.(BrainRegion.ref, Gender), .fun = function(x){
    dplyr::select(x, hgnc_symbol) %>%
      unlist %>% unique()
  })
names(tmp) = c("CER.AD-CONTROL.FEMALE.UP", "CER.AD-CONTROL.MALE.UP", 
               "TCX.AD-CONTROL.FEMALE.UP", "TCX.AD-CONTROL.MALE.UP")
gs1 = c(gs1, tmp)

tmp = diff.exp$MAYO %>%
  tidyr::unite(Gender, Gender.ref, Gender.ag, sep = '_') %>%
  dplyr::filter(BrainRegion.ref %in% c('CER', 'TCX'), 
                Model == 'Diagnosis.Gender',
                Comparison == 'AD-Control', 
                Gender %in% c('M_M', 'F_F'),
                adj.P.Val <= 0.05, logFC <= 0) %>%
  plyr::dlply(.(BrainRegion.ref, Gender), .fun = function(x){
    dplyr::select(x, hgnc_symbol) %>%
      unlist %>% unique()
  })
names(tmp) = c("CER.AD-CONTROL.FEMALE.DOWN", "CER.AD-CONTROL.MALE.DOWN", 
               "TCX.AD-CONTROL.FEMALE.DOWN", "TCX.AD-CONTROL.MALE.DOWN")
gs1 = c(gs1, tmp)
############################################

############################################
# ApoE
gs2 = list()
tmp = diff.exp$ROSMAP %>%
  dplyr::filter(Model == 'Apoe', 
                Comparison %in% c("apoe_genotype2-apoe_genotype0", 
                                  "apoe_genotype2-apoe_genotype1",
                                  "apoe_genotype1-apoe_genotype0"), 
                adj.P.Val <= 0.05, logFC >= 0) %>%
  plyr::dlply(.(Comparison), .fun = function(x){
    select(x, hgnc_symbol) %>%
      unlist %>% unique()
  })

tmp = diff.exp$ROSMAP %>%
  dplyr::filter(Model == 'Apoe', 
                Comparison %in% c("apoe_genotype2-apoe_genotype0", 
                                  "apoe_genotype2-apoe_genotype1",
                                  "apoe_genotype1-apoe_genotype0"), 
                adj.P.Val <= 0.05, logFC <= 0) %>%
  plyr::dlply(.(Comparison), .fun = function(x){
    select(x, hgnc_symbol) %>%
      unlist %>% unique()
  })

tmp = diff.exp$MSSM %>%
  dplyr::filter(Model == 'ApoE', 
                BrainRegion %in% c('IFG', 'PHG' , 'FP', 'STG'),
                Comparison %in% c('0-1','0-2','1-2'), 
                adj.P.Val <= 0.05, logFC >= 0) %>%
  plyr::dlply(.(BrainRegion, Comparison), .fun = function(x){
    dplyr::select(x, hgnc_symbol) %>%
      unlist %>% unique()
  })
names(tmp) = c("PHG.ApoE2-ApoE0.DOWN", "PHG.ApoE2-ApoE1.DOWN",
               "STG.ApoE2-ApoE0.DOWN", "STG.ApoE2-ApoE1.DOWN")
gs2 = c(gs2, tmp)

tmp = diff.exp$MSSM %>%
  dplyr::filter(Model == 'ApoE', 
                BrainRegion %in% c('IFG', 'PHG' , 'FP', 'STG'),
                Comparison %in% c('0-1','0-2','1-2'), 
                adj.P.Val <= 0.05, logFC <= 0) %>%
  plyr::dlply(.(BrainRegion, Comparison), .fun = function(x){
    dplyr::select(x, hgnc_symbol) %>%
      unlist %>% unique()
  })
names(tmp) = c("PHG.ApoE2-ApoE0.UP", "PHG.ApoE2-ApoE1.UP",
               "STG.ApoE2-ApoE1.UP")
gs2 = c(gs2, tmp)

tmp = diff.exp$MAYO %>%
  dplyr::filter(Model == 'ApoE', 
                BrainRegion.ref %in% c('CER', 'TCX'),
                Comparison %in% c('ApoE2-ApoE0','ApoE2-ApoE1','ApoE1-ApoE0'), 
                adj.P.Val <= 0.05, logFC >= 0) %>%
  plyr::dlply(.(BrainRegion.ref, Comparison), .fun = function(x){
    dplyr::select(x, hgnc_symbol) %>%
      unlist %>% unique()
  })
names(tmp) = c("CER.ApoE2-ApoE0.UP", 
               "TCX.ApoE1-ApoE0.UP", "TCX.ApoE2-ApoE0.UP")
gs2 = c(gs2, tmp)

tmp = diff.exp$MAYO %>%
  dplyr::filter(Model == 'ApoE', 
                BrainRegion.ref %in% c('CER', 'TCX'),
                Comparison %in% c('ApoE2-ApoE0','ApoE2-ApoE1','ApoE1-ApoE0'), 
                adj.P.Val <= 0.05, logFC <= 0) %>%
  plyr::dlply(.(BrainRegion.ref, Comparison), .fun = function(x){
    dplyr::select(x, hgnc_symbol) %>%
      unlist %>% unique()
  })
names(tmp) = c("CER.ApoE2-ApoE0.DOWN", 
               "TCX.ApoE1-ApoE0.DOWN", "TCX.ApoE2-ApoE0.DOWN", "TCX.ApoE2-ApoE1.DOWN")
gs2 = c(gs2, tmp)
############################################

############################################
# AOD
gs3 = list()
tmp = diff.exp$ROSMAP %>%
  dplyr::filter(Model == 'Diagnosis.AOD', 
                Comparison %in% c("AOD.AD-AOD.CONTROL"), 
                adj.P.Val <= 0.05, logFC >= 0) %>%
  plyr::dlply(.(Comparison), .fun = function(x){
    select(x, hgnc_symbol) %>%
      unlist %>% unique()
  })
names(tmp) = 'DLPFC.AD-CONTROL.AOD.UP'
gs3 = c(gs3, tmp)

tmp = diff.exp$ROSMAP %>%
  dplyr::filter(Model == 'Diagnosis.AOD', 
                Comparison %in% c("AOD.AD-AOD.CONTROL"), 
                adj.P.Val <= 0.05, logFC <= 0) %>%
  plyr::dlply(.(Comparison), .fun = function(x){
    select(x, hgnc_symbol) %>%
      unlist %>% unique()
  })
names(tmp) = 'DLPFC.AD-CONTROL.AOD.DOWN'
gs3 = c(gs3, tmp)

tmp = diff.exp$MSSM %>%
  dplyr::filter(Model == 'Diagnosis.AOD', 
                Comparison %in% c("IFG.AD.AOD-IFG.CONTROL.AOD",
                                  "STG.AD.AOD-STG.CONTROL.AOD",
                                  "PHG.AD.AOD-PHG.CONTROL.AOD",
                                  "FP.AD.AOD-FP.CONTROL.AOD"), 
                adj.P.Val <= 0.05, logFC >= 0) %>%
  plyr::dlply(.(Comparison), .fun = function(x){
    dplyr::select(x, hgnc_symbol) %>%
      unlist %>% unique()
  })
names(tmp) = c("FP.AD-CONTROL.AOD.UP", "IFG.AD-CONTROL.AOD.UP",
               "PHG.AD-CONTROL.AOD.UP", "STG.AD-CONTROL.AOD.UP")
gs3 = c(gs3, tmp)

tmp = diff.exp$MSSM %>%
  dplyr::filter(Model == 'Diagnosis.AOD', 
                Comparison %in% c("IFG.AD.AOD-IFG.CONTROL.AOD",
                                  "STG.AD.AOD-STG.CONTROL.AOD",
                                  "PHG.AD.AOD-PHG.CONTROL.AOD",
                                  "FP.AD.AOD-FP.CONTROL.AOD"), 
                adj.P.Val <= 0.05, logFC <= 0) %>%
  plyr::dlply(.(Comparison), .fun = function(x){
    dplyr::select(x, hgnc_symbol) %>%
      unlist %>% unique()
  })
names(tmp) = c("FP.AD-CONTROL.AOD.DOWN", "IFG.AD-CONTROL.AOD.DOWN",
               "PHG.AD-CONTROL.AOD.DOWN", "STG.AD-CONTROL.AOD.DOWN")
gs3 = c(gs3, tmp)
############################################

# Get github commit
thisFileName <- 'collateDiffExp.R'
thisRepo <- getRepo(repository = "th1vairam/Brain_Reg_Net", ref="branch", refName='AMPAD_reprocessing')
thisFile <- getPermlink(repository = thisRepo, repositoryPath=paste0('code/R/',thisFileName))

# Store results in synapse
amp.ad.de.geneSets = c(gs, gs1, gs2, gs3)
save(list='amp.ad.de.geneSets', file = 'all.diff.exp.gs.RData')
obj = File('all.diff.exp.gs.RData', name = 'All Differential Expression GeneSets (RData format)', parentId = 'syn8672415')
obj = synStore(obj, used = c('syn8468023', 'syn10157628', 'syn8456721'), executed = thisFile, activityName = 'Collate differential expression results')
