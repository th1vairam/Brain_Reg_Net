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
all.used.ids = c('syn8468023', 'syn10157628', 'syn8456721')
############################################
# AD-CONTROL
all.de = rbindlist(list(
  diff.exp$ROSMAP %>%
    dplyr::filter(Model == 'Diagnosis', Comparison == 'AD-CONTROL') %>%
    dplyr::mutate(Study = 'ROSMAP') %>%
    dplyr::select(-Model, -Comparison) %>%
    dplyr::rename(Tissue = Region),
  diff.exp$MSSM %>%
    dplyr::filter(Model == 'Diagnosis', BrainRegion %in% c('IFG', 'PHG' , 'FP', 'STG'), Comparison == 'CONTROL-AD') %>%
    dplyr::mutate(Study = 'MSSM', logFC = -logFC, CI.L = -CI.L, CI.R = -CI.R) %>%
    dplyr::select(-reference, -against, -Model, -Comparison) %>%
    dplyr::rename(Tissue = BrainRegion, CI.L = CI.R, CI.R = CI.L),
  diff.exp$MAYO %>%
    dplyr::filter(BrainRegion.ref %in% c('CER', 'TCX'), Model == 'Diagnosis', Comparison == 'Control-AD') %>%
    dplyr::mutate(Study = 'MAYO', logFC = -logFC, CI.L = -CI.L, CI.R = -CI.R) %>%
    dplyr::select(-Model, -Comparison, -BrainRegion.ag, -Gender.ref, -Gender.ag, -AOD.ref, -AOD.ag) %>%
    dplyr::rename(Tissue = BrainRegion.ref, CI.L = CI.R, CI.R = CI.L)),
  use.names = T, fill = T) %>%
  dplyr::mutate(Model = 'Diagnosis', Gender = 'ALL', Comparison = 'AD-CONTROL')
all.de$Direction[all.de$logFC > 0] = 'UP'
all.de$Direction[all.de$logFC < 0] = 'DOWN'
############################################

############################################
# AD - CONTROL (Male vs Female)
gender.de = rbindlist(list(
  diff.exp$ROSMAP %>%
    dplyr::filter(Model == 'Diagnosis.Gender', Comparison %in% c("AD-CONTROL.IN.FEMALE", "AD-CONTROL.IN.MALE")) %>%
    dplyr::mutate(Study = 'ROSMAP', Comparison = gsub('AD-CONTROL.IN.','',Comparison)) %>%
    dplyr::select(-Model) %>%
    dplyr::rename(Tissue = Region, Gender = Comparison),
  diff.exp$MSSM %>%
    dplyr::filter(Model == 'Diagnosis.SEX',
                  Comparison %in% c('FP.AD.FEMALE-FP.CONTROL.FEMALE','FP.AD.MALE-FP.CONTROL.MALE',
                                    'IFG.AD.FEMALE-IFG.CONTROL.FEMALE','IFG.AD.MALE-IFG.CONTROL.MALE',
                                    'PHG.AD.FEMALE-PHG.CONTROL.FEMALE','PHG.AD.MALE-PHG.CONTROL.MALE',
                                    'STG.AD.FEMALE-STG.CONTROL.FEMALE','STG.AD.MALE-STG.CONTROL.MALE')) %>%
    dplyr::mutate(Study = 'MSSM') %>%
    tidyr::separate(Comparison, c('ref', 'ag'), sep = '\\-') %>%
    tidyr::separate(ref, c('Tissue', 'a', 'Gender'), sep = '\\.') %>%
    dplyr::select(-reference, -against, -Model, -ag, -a, -BrainRegion),
  diff.exp$MAYO %>%
    tidyr::unite(Gender, Gender.ref, Gender.ag, sep = '_') %>%
    dplyr::filter(BrainRegion.ref %in% c('CER', 'TCX'), Model == 'Diagnosis.Gender', Comparison == 'AD-Control', Gender %in% c('M_M', 'F_F')) %>%
    dplyr::mutate(Study = 'MAYO', Gender = gsub('M_M', 'MALE', Gender), Gender = gsub('F_F', 'FEMALE', Gender)) %>%
    dplyr::select(-Model, -Comparison, -BrainRegion.ag, -AOD.ref, -AOD.ag) %>%
    dplyr::rename(Tissue = BrainRegion.ref)),
  use.names = T, fill = T) %>%
  dplyr::mutate(Model = 'Diagnosis.Gender', Comparison = 'AD-CONTROL')
gender.de$Direction[gender.de$logFC > 0] = 'UP'
gender.de$Direction[gender.de$logFC < 0] = 'DOWN'
############################################

############################################
# AD - CONTROL (AOD)
aod.de = rbindlist(list(
  diff.exp$ROSMAP %>%
    dplyr::filter(Model == 'Diagnosis.AOD', Comparison %in% c("AOD.AD-AOD.CONTROL")) %>%
    dplyr::mutate(Study = 'ROSMAP', Comparison = gsub('AOD.','',Comparison)) %>%
    dplyr::select(-Model) %>%
    dplyr::rename(Tissue = Region),
  diff.exp$MSSM %>%
    dplyr::filter(Model == 'Diagnosis.AOD',
                  Comparison %in% c("FP.AD.AOD-FP.CONTROL.AOD",
                                    "IFG.AD.AOD-IFG.CONTROL.AOD",
                                    "PHG.AD.AOD-PHG.CONTROL.AOD",
                                    "STG.AD.AOD-STG.CONTROL.AOD")) %>%
    dplyr::mutate(Study = 'MSSM') %>%
    tidyr::separate(Comparison, c('ref', 'ag'), sep = '\\-') %>%
    tidyr::separate(ref, c('Tissue', 'a', 'b'), sep = '\\.') %>%
    dplyr::select(-reference, -against, -Model, -ag, -a, -BrainRegion, -b),
  diff.exp$MAYO %>%
    dplyr::filter(BrainRegion.ref %in% c('CER', 'TCX'), Model == 'Diagnosis.AOD', Comparison == 'AD-Control') %>%
    dplyr::mutate(Study = 'MAYO') %>%
    dplyr::select(-Model, -Comparison, -BrainRegion.ag, -AOD.ref, -AOD.ag, -Gender.ref, -Gender.ag) %>%
    dplyr::rename(Tissue = BrainRegion.ref)),
  use.names = T, fill = T) %>%
  dplyr::mutate(Model = 'Diagnosis.AOD', Gender = 'ALL', Comparison = 'AD-CONTROL')
aod.de$Direction[aod.de$logFC > 0] = 'UP'
aod.de$Direction[aod.de$logFC < 0] = 'DOWN'
############################################

############################################
# ApoE4 allele
apoe.de = rbindlist(list(
  diff.exp$ROSMAP %>%
    dplyr::filter(Model == 'Apoe', Comparison %in% c("apoe_genotype2-apoe_genotype0", "apoe_genotype2-apoe_genotype1", "apoe_genotype1-apoe_genotype0")) %>%
    dplyr::mutate(Study = 'ROSMAP', Comparison = gsub('apoe_genotype','',Comparison)) %>%
    dplyr::select(-Model) %>%
    dplyr::rename(Tissue = Region),
  diff.exp$MSSM %>%
    dplyr::filter(Model == 'ApoE',
                  Comparison %in% c("0-1", "0-2", "1-2")) %>%
    dplyr::mutate(Study = 'MSSM', logFC = -logFC, CI.L = -CI.L, CI.R = -CI.R) %>%
    tidyr::separate(Comparison, c('ref', 'ag'), sep = '\\-') %>%
    tidyr::unite(Comparison, ag, ref, sep = '-') %>%
    dplyr::select(-reference, -against) %>%
    dplyr::rename(Tissue = BrainRegion, CI.L = CI.R, CI.R = CI.L),
  diff.exp$MAYO %>%
    dplyr::filter(BrainRegion.ref %in% c('CER', 'TCX'), Model == 'ApoE', Comparison %in% c("ApoE2-ApoE0", "ApoE2-ApoE1", "ApoE1-ApoE0")) %>%
    dplyr::mutate(Study = 'MAYO', Comparison = gsub('ApoE','', Comparison)) %>%
    dplyr::select(-BrainRegion.ag, -AOD.ref, -AOD.ag, -Gender.ref, -Gender.ag) %>%
    dplyr::rename(Tissue = BrainRegion.ref)),
  use.names = T, fill = T) %>%
  dplyr::mutate(Model = 'ApoE4', Gender = 'ALL')
apoe.de$Direction[apoe.de$logFC > 0] = 'UP'
apoe.de$Direction[apoe.de$logFC < 0] = 'DOWN'
############################################

############################################
# Get github commit
thisFileName <- 'code/R/collateDifferentialExpressionResults.R'
thisRepo <- getRepo(repository = "th1vairam/Brain_Reg_Net", ref="branch", refName='AMPAD_reprocessing')
thisFile <- getPermlink(repository = thisRepo, repositoryPath=thisFileName)
############################################

############################################
# Store in synapse
tmp = rbindlist(list(all.de, gender.de, aod.de, apoe.de), use.names = T, fill = T) %>%
  write.table('differentialExpressionSummary.tsv', sep = '\t', row.names = F, quote = F)
obj = File('differentialExpressionSummary.tsv', name = 'All Differential Expression (Merged)', parentId = 'syn8672415')
obj = synStore(obj, used = all.used.ids, executed = thisFile, activityName = 'Merge all differential expression')
############################################