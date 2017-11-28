## Plot differential expression results

# Load libraries
library(data.table)
library(tidyr)
library(plyr)
library(dplyr)
library(CovariateAnalysis) # devtools::install_github(th1vairam/CovariateAnalysis@dev)

library(ggplot2)
library(ggpubr)

library(synapseClient)
library(githubr)

synapseLogin()

# Set ggplot theme
my_theme <- function (base_size = 11, base_family = ""){
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA, colour = "grey20"), 
          panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "grey92", size = 0.25), 
          strip.background = element_rect(fill = "grey85", colour = "grey20"), 
          legend.key = element_rect(fill = "white", colour = NA), complete = TRUE,
          legend.position = 'top',
          plot.title=element_text(hjust=0.5))
}

# Get differential expression results from synapse
de = downloadFile('syn11180450') %>%
  dplyr::filter(Gender == 'ALL', Comparison == 'AD-CONTROL', Model == 'Diagnosis')
de$Direction = 'NONE'
de$Direction[de$logFC >= log2(1.2) & de$adj.P.Val <= 0.05] = 'UP'
de$Direction[de$logFC <= -log2(1.2) & de$adj.P.Val <= 0.05] = 'DOWN'

p = ggplot(de, aes(x = logFC, y = -log10(adj.P.Val), color = Direction)) + geom_point() + facet_grid(. ~ Tissue + Study)
p = p + scale_color_manual(values = c('green','grey','red')) + my_theme() + ggtitle('(a) Differential Expression')

t1 = de %>%
  dplyr::group_by(Tissue, Direction) %>%
  dplyr::summarise(count = length(unique(ensembl_gene_id))) %>%
  tidyr::spread(Direction, count) %>%
  ggtexttable(rows = NULL) + ggtitle('(e) Number of differentially expressed genes')

# Forest plots of most upregulated and most down regulated gene
de$Direction[de$logFC >= log2(1) & de$adj.P.Val <= 0.05] = 'UP'
de$Direction[de$logFC <= -log2(1) & de$adj.P.Val <= 0.05] = 'DOWN'
## Most up
de.up = dplyr::filter(de, Direction == 'UP') %>%
  dplyr::group_by(ensembl_gene_id) %>%
  summarise(count = length(unique(Tissue))) %>%
  dplyr::filter(count >= 7)
de.up = dplyr::filter(de, Direction == 'UP', ensembl_gene_id %in% de.up$ensembl_gene_id) %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::summarise(all.fc = sum(logFC)) %>%
  dplyr::arrange(desc(all.fc))
de.up = dplyr::filter(de, ensembl_gene_id %in% de.up$ensembl_gene_id[1])
p1 = de.up %>%
  ggplot(aes(x=logFC, y=Tissue, color = 'red')) + geom_point(aes(size=-log10(adj.P.Val)), color = 'red')
p1 = p1 + geom_errorbarh(aes(xmin=CI.L, xmax=CI.R), height = 0.2) + geom_vline(xintercept=0, linetype="dashed")
p1 = p1 + my_theme() + scale_color_manual(values = 'red', guide = FALSE) + ggtitle(paste('(b)', unique(de.up$hgnc_symbol)))
p1 = p1 + scale_size_continuous(range = c(0,5))
p1

## Most down
de.down = dplyr::filter(de, Direction == 'DOWN', hgnc_symbol != '') %>%
  dplyr::group_by(ensembl_gene_id) %>%
  summarise(count = length(unique(Tissue))) %>%
  dplyr::filter(count >= 7)
de.down = dplyr::filter(de, Direction == 'DOWN', ensembl_gene_id %in% de.down$ensembl_gene_id) %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::summarise(all.fc = sum(logFC)) %>%
  dplyr::arrange(desc(all.fc))
de.down = dplyr::filter(de, ensembl_gene_id %in% de.down$ensembl_gene_id[1])
p2 =  de.down %>%
  ggplot(aes(x=logFC, y=Tissue, color = 'darkgreen')) + geom_point(aes(size=-log10(adj.P.Val)), color = 'darkgreen')
p2 = p2 + geom_errorbarh(aes(xmin=CI.L, xmax=CI.R), height = 0.2) + geom_vline(xintercept=0, linetype="dashed")
p2 = p2 + my_theme() + scale_color_manual(values = 'darkgreen', guide = FALSE) + ggtitle(paste('(c)', unique(de.down$hgnc_symbol)))
p2 = p2  + scale_size_continuous(range = c(0,5))
p2

## Most varied
de.var = dplyr::filter(de) %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::filter(count >= 7)
de.var = dplyr::filter(de, ensembl_gene_id %in% de.var$ensembl_gene_id) %>%
  dplyr::select(-hgnc_symbol) %>% unique() %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::summarise(variability = var(logFC)) %>%
  dplyr::arrange((abs(variability)))
de.var = dplyr::filter(de, ensembl_gene_id %in% de.var$ensembl_gene_id[1])
p3 =  de.var %>%
  ggplot(aes(x=logFC, y=Tissue, color = Direction)) + geom_point(aes(size=-log10(adj.P.Val)))
p3 = p3 + geom_errorbarh(aes(xmin=CI.L, xmax=CI.R), height = 0.2) + geom_vline(xintercept=0, linetype="dashed")
p3 = p3 + my_theme() + scale_color_manual(values = c('grey'), guide = FALSE) + ggtitle(paste('(d)',unique(de.var$hgnc_symbol)))
p3 = p3 + scale_size_continuous(range = c(0,5))
p3

p4 = ggpubr::ggarrange(plotlist = list(t1, p1,p2,p3), ncol = 4, nrow = 1)
p4
p5 = ggpubr::ggarrange(plotlist = list(p,p4), ncol = 1, nrow = 2)
svg(file = 'DEFigure.svg', height = 10, width = 14)
p5
dev.off()
###################################################################
## Get all differential expression tables
# Get differential expression results from synapse
de = downloadFile('syn11180450') 
de$Direction = 'NONE'
de$Direction[de$logFC >= log2(1.2) & de$adj.P.Val <= 0.05] = 'UP'
de$Direction[de$logFC <= -log2(1.2) & de$adj.P.Val <= 0.05] = 'DOWN'

t.all = de %>%
  dplyr::filter(Gender == 'ALL', Comparison == 'AD-CONTROL', Model == 'Diagnosis') %>%
  dplyr::group_by(Tissue, Direction) %>%
  dplyr::summarise(count = length(unique(ensembl_gene_id))) %>%
  tidyr::spread(Direction, count) %>%
  ggtexttable(rows = NULL)

t.aod = de %>%
  dplyr::filter(Gender == 'ALL', Comparison == 'AD-CONTROL', Model == 'Diagnosis.AOD') %>%
  dplyr::group_by(Tissue, Direction) %>%
  dplyr::summarise(count = length(unique(ensembl_gene_id))) %>%
  tidyr::spread(Direction, count) %>%
  ggtexttable(rows = NULL)

t.male = de %>%
  dplyr::filter(Gender == 'MALE', Comparison == 'AD-CONTROL', Model == 'Diagnosis.Gender') %>%
  dplyr::group_by(Tissue, Direction) %>%
  dplyr::summarise(count = length(unique(ensembl_gene_id))) %>%
  tidyr::spread(Direction, count) %>%
  ggtexttable(rows = NULL)

t.female = de %>%
  dplyr::filter(Gender == 'FEMALE', Comparison == 'AD-CONTROL', Model == 'Diagnosis.Gender') %>%
  dplyr::group_by(Tissue, Direction) %>%
  dplyr::summarise(count = length(unique(ensembl_gene_id))) %>%
  tidyr::spread(Direction, count) %>%
  ggtexttable(rows = NULL)

t.20 = de %>%
  dplyr::filter(Gender == 'ALL', Comparison == '2-0', Model == 'ApoE4') %>%
  dplyr::group_by(Tissue, Direction) %>%
  dplyr::summarise(count = length(unique(ensembl_gene_id))) %>%
  tidyr::spread(Direction, count) %>%
  ggtexttable(rows = NULL)

t.21 = de %>%
  dplyr::filter(Gender == 'ALL', Comparison == '2-1', Model == 'ApoE4') %>%
  dplyr::group_by(Tissue, Direction) %>%
  dplyr::summarise(count = length(unique(ensembl_gene_id))) %>%
  tidyr::spread(Direction, count) %>%
  ggtexttable(rows = NULL)

t.10 = de %>%
  dplyr::filter(Gender == 'ALL', Comparison == '1-0', Model == 'ApoE4') %>%
  dplyr::group_by(Tissue, Direction) %>%
  dplyr::summarise(count = length(unique(ensembl_gene_id))) %>%
  tidyr::spread(Direction, count) %>%
  ggtexttable(rows = NULL)