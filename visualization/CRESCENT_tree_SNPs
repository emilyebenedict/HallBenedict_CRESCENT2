# File name:      CRESCENT_figures_eeb.R
# Author:         Emily Benedict, ebenedict@wustl.edu
# Last Modified:  2026-01-30
# Description:    This script generates figures for CRESCENT isolates

# loading packages
library(readxl)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(NatParksPalettes)
library(RColorBrewer)
library(phytools)
library(ggtree)
library(pheatmap)
library(tidytree)
library(ggnewscale)
library(ggprism)
library(igraph)

# setting wd
setwd('~/Box Sync/EEB/Dantas/Crescent/')

#### Reading in trees ####
# all CRESCENT isolates
tree = ape::read.tree('./manuscript/Github_scripts/data/CRESCENT_all_core_genome_tree')
tree$tip.label = gsub('K','K-',tree$tip.label)
tree$tip.label = gsub('--','-',tree$tip.label)
# CRESCENT Escherichia isolates
esch_tree = ape::read.tree('./manuscript/Github_scripts/data/CRESCENT_escerichia_core_genome_tree_eeb')
esch_tree$tip.label = gsub('K','K-',esch_tree$tip.label)
esch_tree$tip.label = gsub('--','-',esch_tree$tip.label)
# CRESCENT Enterobacter isolates
ent_tree = ape::read.tree('./manuscript/Github_scripts/data/CRESCENT_enterobacter_core_genome_tree_eeb')
ent_tree$tip.label = gsub('K','K-',ent_tree$tip.label)
ent_tree$tip.label = gsub('--','-',ent_tree$tip.label)
# CRESCENT Klebsiella isolates
kleb_tree = ape::read.tree('./manuscript/Github_scripts/data/CRESCENT_klebsiella_core_genome_tree_eeb')
kleb_tree$tip.label = gsub('K','K-',kleb_tree$tip.label)
kleb_tree$tip.label = gsub('--','-',kleb_tree$tip.label)

# CRESCENT Escherichia isolates + NCBI
ext_esch_tree = ape::read.tree('./manuscript/Github_scripts/data/CRESCENT_NCBI_escherichia_core_genome_tree')
ext_esch_tree$tip.label = gsub('K','K-',ext_esch_tree$tip.label)
ext_esch_tree$tip.label = gsub('--','-',ext_esch_tree$tip.label)
# CRESCENT Klebsiella isolates + NCBI
ext_kleb_tree = ape::read.tree('./manuscript/Github_scripts/data/CRESCENT_NCBI_klebsiella_core_genome_tree')
ext_kleb_tree$tip.label = gsub('K','K-',ext_kleb_tree$tip.label)
ext_kleb_tree$tip.label = gsub('--','-',ext_kleb_tree$tip.label)
# CRESCENT Enterobacter isolates + NCBI
ext_ent_tree = ape::read.tree('./manuscript/Github_scripts/data/CRESCENT_NCBI_enterobacter_core_genome_tree')
ext_ent_tree$tip.label = gsub('K','K-',ext_ent_tree$tip.label)
ext_ent_tree$tip.label = gsub('--','-',ext_ent_tree$tip.label)

#### Reading in data ####
# reading in CRESCENT study participants and isolates to include
final_list = read.csv('./manuscript/Github_scripts/data/CRESCENT_pts_and_samples.csv', stringsAsFactors = F)
final_pts = final_list
final_list = final_list[!is.na(final_list$WGS_ID),,drop = F]
final_list$WGS_ID = gsub('K','K-',final_list$WGS_ID)
final_list$WGS_ID = gsub('--','-',final_list$WGS_ID)
tree = drop.tip(tree, setdiff(tree$tip.label, final_list$WGS_ID))
# dropping any samples not present in the final CRESCENT core genome tree
final_list = final_list[final_list$WGS_ID %in% tree$tip.label,,drop = F]

# reading in colors
colors = data.frame(read_excel('./manuscript/Github_scripts/data/colors.xlsx'), stringsAsFactors = F)

all_mlst_colors = colors$color[grepl('MLST',colors$topic)]
names(all_mlst_colors) = colors$name[grepl('MLST',colors$topic)]

carb_colors = c(colors$color[colors$topic == 'Carbapenemase'])
names(carb_colors) = c(colors$name[colors$topic == 'Carbapenemase'])

# reading in pt/genomic data
pt_data = data.frame(read_excel('./manuscript/Supplemental_Tables/250828/250828_SupplementalTable3-Genomics_LRH.xlsx'), stringsAsFactors = F)
pt_data$Sample = gsub('EG','EG-',pt_data$Sample)
colnames(final_list)[ncol(final_list)] = 'Sample'
pt_data = left_join(pt_data, final_list)
pt_data$Category[pt_data$Category == 'Hospital-acquired'] = 'Hospital acquired'
pt_data$genus = gsub(' .*','',pt_data$WGS.identification)

# rooting trees
tree = drop.tip(tree, setdiff(tree$tip.label, pt_data$Sample))
tree = midpoint.root(tree)

kleb_tree = drop.tip(kleb_tree, setdiff(kleb_tree$tip.label, pt_data$Sample))
kleb_tree = midpoint.root(kleb_tree)

ent_tree = drop.tip(ent_tree, setdiff(ent_tree$tip.label, pt_data$Sample))
ent_tree = midpoint.root(ent_tree)

esch_tree = drop.tip(esch_tree, setdiff(esch_tree$tip.label, pt_data$Sample))
esch_tree = midpoint.root(esch_tree)
# plotting tree
treeplot = ggtree(tree)
circle_treeplot = ggtree(tree, layout = 'circular')

# making pt_dataframe in order of tree$tiplabels
pt_df = data.frame(matrix(NA, nrow = length(tree$tip.label), ncol = 1), stringsAsFactors = F)
colnames(pt_df) = 'Sample'
pt_df$Sample = tree$tip.label
for(i in 1:length(tree$tip.label)){
  pt_df$sink_pt[i] = pt_data$Sink_or_patient[pt_data$Sample == pt_df$Sample[i]]
  pt_df$mash[i] = pt_data$WGS.identification[pt_data$Sample == pt_df$Sample[i]]
  pt_df$category[i] = pt_data$Category[pt_data$Sample == pt_df$Sample[i]]
}
# breaking these into specific plotting dfs
sinkpt_df = pt_df[,'sink_pt', drop = F]
rownames(sinkpt_df) = pt_df$Sample
colnames(sinkpt_df) = 'Sink_patient'

mash_df = pt_df[,c('Sample','mash')]
rownames(mash_df) = mash_df$Sample
colnames(mash_df) = c('Sample','Mash')

epi_df = pt_df[,'category', drop = F]
rownames(epi_df) = pt_df$Sample
colnames(epi_df) = 'Category'

# reading in ASTs
ast_data = read.csv('./manuscript/Github_scripts/data/CRESECENT_asts.csv', stringsAsFactors = F)
# ordering ASTs for tree annotation
ast_df = data.frame(matrix(NA, nrow = length(tree$tip.label), ncol = (ncol(ast_data) - 2)))
colnames(ast_df) = c('Sample',colnames(ast_data)[c(1:(ncol(ast_data) - 3))])
colnames(ast_df) = gsub('_Interpretation','', colnames(ast_df))
ast_df$Sample = tree$tip.label
for(i in 1:nrow(ast_df)){
  if(ast_df$Sample[i] %in% ast_data$Sample){
    ast_df[i,2:ncol(ast_df)] = ast_data[ast_data$Sample == ast_df$Sample[i], 1:(ncol(ast_data)-3)]
  } 
}

ast_df[is.na(ast_df)] = 'no data'

# adding astdf to ptdf
pt_df = full_join(pt_df, ast_df)
# formatting astdf for plotting
rownames(ast_df) = ast_df$Sample
ast_df = ast_df[,-1]

# separating imipenem and meropenem for tree annotation
carbast_df = ast_df[,c('Imipenem','Meropenem'), drop = F]
rownames(carbast_df) = rownames(ast_df)
colnames(carbast_df) = c('Imipenem','Meropenem')

# reading in amrfinder results
amrfinder_data = data.frame(read_excel('./manuscript/Github_scripts/data/SupplementalTable3-Genomics.xlsx', sheet = 3), stringsAsFactors = F)
amrfinder_data$Name.x = gsub('K_','K-',amrfinder_data$Name.x)
amrfinder_data$Name.x = gsub('K','K-',amrfinder_data$Name.x)
amrfinder_data$Name.x = gsub('--','-',amrfinder_data$Name.x)
amrfinder_data$Name.x = gsub('EG_','EG-',amrfinder_data$Name.x)
amrfinder_data = amrfinder_data[amrfinder_data$Name.x %in% tree$tip.label,]
# correcting amrfinderdata's Category column
for(i in 1:nrow(amrfinder_data)){
  amrfinder_data$Category[i] = pt_data$Category[pt_data$Sample == amrfinder_data$Name.x[i]]
}

# adding carbapenemase data to all_data
carb_data = data.frame(matrix('not present', nrow = nrow(pt_data), ncol = (3 + length(table(amrfinder_data$Gene.symbol[amrfinder_data$Subclass == 'CARBAPENEM'])))), 
                       stringsAsFactors = F)
colnames(carb_data) = c('Sample','blaNDM','blaNDM1','blaNDM5','blaNDM7','blaOXA181','blaOXA232','blaOXA484','blaNDM_all','blaOXA_all')
carb_data$Sample = pt_data$Sample  
for(i in 1:nrow(carb_data)){
  if(carb_data$Sample[i] %in% amrfinder_data$Name.x[amrfinder_data$Gene.symbol == 'blaNDM']){
    carb_data$blaNDM[i] = amrfinder_data$Gene.symbol[amrfinder_data$Gene.symbol == 'blaNDM' & amrfinder_data$Name.x == carb_data$Sample[i]][1]
    carb_data$blaNDM_all[i] = amrfinder_data$Gene.symbol[amrfinder_data$Gene.symbol == 'blaNDM' & amrfinder_data$Name.x == carb_data$Sample[i]][1]
  }
  
  if(carb_data$Sample[i] %in% amrfinder_data$Name.x[amrfinder_data$Gene.symbol == 'blaNDM-1']){
    carb_data$blaNDM1[i] = amrfinder_data$Gene.symbol[amrfinder_data$Gene.symbol == 'blaNDM-1' & amrfinder_data$Name.x == carb_data$Sample[i]][1]
    carb_data$blaNDM_all[i] = amrfinder_data$Gene.symbol[amrfinder_data$Gene.symbol == 'blaNDM-1' & amrfinder_data$Name.x == carb_data$Sample[i]][1]
  }
  
  if(carb_data$Sample[i] %in% amrfinder_data$Name.x[amrfinder_data$Gene.symbol == 'blaNDM-5']){
    carb_data$blaNDM5[i] = amrfinder_data$Gene.symbol[amrfinder_data$Gene.symbol == 'blaNDM-5' & amrfinder_data$Name.x == carb_data$Sample[i]][1]
    carb_data$blaNDM_all[i] = amrfinder_data$Gene.symbol[amrfinder_data$Gene.symbol == 'blaNDM-5' & amrfinder_data$Name.x == carb_data$Sample[i]][1]
  }
  
  if(carb_data$Sample[i] %in% amrfinder_data$Name.x[amrfinder_data$Gene.symbol == 'blaNDM-7']){
    carb_data$blaNDM7[i] = amrfinder_data$Gene.symbol[amrfinder_data$Gene.symbol == 'blaNDM-7' & amrfinder_data$Name.x == carb_data$Sample[i]][1]
    carb_data$blaNDM_all[i] = amrfinder_data$Gene.symbol[amrfinder_data$Gene.symbol == 'blaNDM-7' & amrfinder_data$Name.x == carb_data$Sample[i]][1]
  }
  
  if(carb_data$Sample[i] %in% amrfinder_data$Name.x[amrfinder_data$Gene.symbol == 'blaOXA-181']){
    carb_data$blaOXA181[i] = amrfinder_data$Gene.symbol[amrfinder_data$Gene.symbol == 'blaOXA-181' & amrfinder_data$Name.x == carb_data$Sample[i]][1]
    carb_data$blaOXA_all[i] = amrfinder_data$Gene.symbol[amrfinder_data$Gene.symbol == 'blaOXA-181' & amrfinder_data$Name.x == carb_data$Sample[i]][1]
  }
  
  if(carb_data$Sample[i] %in% amrfinder_data$Name.x[amrfinder_data$Gene.symbol == 'blaOXA-232']){
    carb_data$blaOXA232[i] = amrfinder_data$Gene.symbol[amrfinder_data$Gene.symbol == 'blaOXA-232' & amrfinder_data$Name.x == carb_data$Sample[i]][1]
    carb_data$blaOXA_all[i] = amrfinder_data$Gene.symbol[amrfinder_data$Gene.symbol == 'blaOXA-232' & amrfinder_data$Name.x == carb_data$Sample[i]][1]
  }
  
  if(carb_data$Sample[i] %in% amrfinder_data$Name.x[amrfinder_data$Gene.symbol == 'blaOXA-484']){
    carb_data$blaOXA484[i] = amrfinder_data$Gene.symbol[amrfinder_data$Gene.symbol == 'blaOXA-484' & amrfinder_data$Name.x == carb_data$Sample[i]][1]
    carb_data$blaOXA_all[i] = amrfinder_data$Gene.symbol[amrfinder_data$Gene.symbol == 'blaOXA-484' & amrfinder_data$Name.x == carb_data$Sample[i]][1]
  }
}

carb_df = carb_data[,c('blaNDM_all','blaOXA_all'), drop = F]
rownames(carb_df) = carb_data$Sample
carb_df[carb_df == ''] = 'not present'
colnames(carb_df) = c('blaNDM','blaOXA')
carb_colors = colors$color[colors$topic == 'Carbapenemase']
names(carb_colors) = colors$name[colors$topic == 'Carbapenemase']

pt_data = full_join(pt_data, carb_data)

colnames(amrfinder_data)[1] = 'Sample'
amrfinder_data = full_join(amrfinder_data, pt_data)

## want to plot x-axis as epi category, y-axis as carbapenemase genes, points sized by count of isolates
carbapenemase_count_df = data.frame(matrix(NA, nrow = 4*length(table(amrfinder_data$Gene.symbol[amrfinder_data$Subclass == 'CARBAPENEM'])),
                                           ncol = 6), stringsAsFactors = F)
colnames(carbapenemase_count_df) = c('Carbapenemase','Status','Freq','EcFreq','EntFreq','KlebFreq')
carbapenemase_count_df$Carbapenemase = rep(names(table(amrfinder_data$Gene.symbol[amrfinder_data$Subclass == 'CARBAPENEM'])),4)
carbapenemase_count_df$Status = c(rep('Hospital-acquired',7),rep('Healthcare-associated',7),rep('Community-associated',7),
                                  rep('Environmental swab',7))
for(i in 1:nrow(carbapenemase_count_df)){
  carbapenemase_count_df$EcFreq[i] = sum(amrfinder_data$Gene.symbol == carbapenemase_count_df$Carbapenemase[i] & amrfinder_data$Category == carbapenemase_count_df$Status[i] &
                                           amrfinder_data$genus == 'Escherichia')
  carbapenemase_count_df$KlebFreq[i] = sum(amrfinder_data$Gene.symbol == carbapenemase_count_df$Carbapenemase[i] & amrfinder_data$Category == carbapenemase_count_df$Status[i] &
                                             amrfinder_data$genus == 'Klebsiella')
  carbapenemase_count_df$EntFreq[i] = sum(amrfinder_data$Gene.symbol == carbapenemase_count_df$Carbapenemase[i] & amrfinder_data$Category == carbapenemase_count_df$Status[i] & 
                                            amrfinder_data$genus == 'Enterobacter')
  carbapenemase_count_df$Freq[i] = sum(amrfinder_data$Gene.symbol == carbapenemase_count_df$Carbapenemase[i] & amrfinder_data$Category == carbapenemase_count_df$Status[i])
}

ndmoxa_df = data.frame(matrix(NA, nrow = 12, ncol = 3), stringsAsFactors = F)
colnames(ndmoxa_df) = c('blaNDM','blaOXA','Freq')
ndmoxa_df$blaOXA = c(rep('blaOXA-181',4), rep('blaOXA-232',4), rep('blaOXA-484',4))
ndmoxa_df$blaNDM = rep(c('blaNDM','blaNDM-1','blaNDM-5','blaNDM-7'),3)
colnames(carb_data) = c('Sample','blaNDMa','blaNDM-1','blaNDM-5','blaNDM-7','blaOXA-181','blaOXA-232','blaOXA-484','blaNDM','blaOXA')
for(i in 1:nrow(ndmoxa_df)){
  ndmoxa_df$Freq[i] = sum(carb_data[colnames(carb_data) == ndmoxa_df$blaNDM[i]] == ndmoxa_df$blaNDM[i] & carb_data[colnames(carb_data) == ndmoxa_df$blaOXA[i]] == ndmoxa_df$blaOXA[i])
}


# reading in plasmidfinder results
# reading in plasmid data
plasmid_data = data.frame(read_excel('./manuscript/Github_scripts/data/SupplementalTable3-Genomics.xlsx', sheet = 4), stringsAsFactors = F)
plasmid_data$Sample = gsub('K','K-',plasmid_data$Sample)
plasmid_data$Sample = gsub('--','-',plasmid_data$Sample)
# subsetting to those in our tree data
plasmid_data = plasmid_data[plasmid_data$Sample %in% tree$tip.label,]
# adding simplified_plasmid column
plasmid_data$simplified_plasmid = gsub('Col.*','Col',plasmid_data$Plasmid)
plasmid_data$simplified_plasmid = gsub('Inc.*','Inc',plasmid_data$simplified_plasmid)
plasmid_data$simplified_plasmid = gsub('FIA.*','FIA',plasmid_data$simplified_plasmid)

# generating plasmid colors
simple_plasmid_colors = c(natparks.pals('Volcanoes',length(table(plasmid_data$simplified_plasmid))), 'grey95')
names(simple_plasmid_colors) = c(names(table(plasmid_data$simplified_plasmid)), 'none')

# adding Inc column
plasmid_data$Inc = plasmid_data$Plasmid
plasmid_data$Inc[plasmid_data$simplified_plasmid != 'Inc'] = 'none'

# generating Inc plasmid colors
inc_plasmid_colors = c(natparks.pals('Olympic',(length(table(plasmid_data$Inc) - 1))), 'grey97')
names(inc_plasmid_colors) = names(table(plasmid_data$Inc))

# adding Col column
plasmid_data$Col = plasmid_data$Plasmid
plasmid_data$Col[plasmid_data$simplified_plasmid != 'Col'] = 'none'

# adding noncol noninc column
plasmid_data$otherplasmids = plasmid_data$simplified_plasmid
plasmid_data$otherplasmids = gsub('Col','none',plasmid_data$otherplasmids)
plasmid_data$otherplasmids = gsub('Inc','none',plasmid_data$otherplasmids)

# generating Col plasmid colors
col_plasmid_colors = c(natparks.pals('GrandCanyon', (length(table(plasmid_data$Col)) - 1)),'grey97')
names(col_plasmid_colors) = c(names(table(plasmid_data$Col)))


# making ordered dfs for plotting
col_plasmid_df = data.frame(matrix(NA, nrow = length(tree$tip.label), ncol = 2))
colnames(col_plasmid_df) = c('Isolate','Col')
col_plasmid_df$Isolate = tree$tip.label
col_plasmid_df$Col = 'none'
for(i in 1:nrow(col_plasmid_df)){
  if('Col' %in% plasmid_data$simplified_plasmid[plasmid_data$Sample == col_plasmid_df$Isolate[i]]){
    col_plasmid_df$Col[i] = 'Col'
  }
}
col_plasmid_df = data.frame(col_plasmid_df$Col)
rownames(col_plasmid_df) = tree$tip.label
colnames(col_plasmid_df) = 'Col'

inc_plasmid_df = data.frame(matrix(NA, nrow = length(tree$tip.label), ncol = 2))
colnames(inc_plasmid_df) = c('Isolate','Inc')
inc_plasmid_df$Isolate = tree$tip.label
inc_plasmid_df$Inc = 'none'
for(i in 1:nrow(inc_plasmid_df)){
  if('Inc' %in% plasmid_data$simplified_plasmid[plasmid_data$Sample == inc_plasmid_df$Isolate[i]]){
    inc_plasmid_df$Inc[i] = 'Inc'
  }
}
inc_plasmid_df = data.frame(inc_plasmid_df$Inc)
rownames(inc_plasmid_df) = tree$tip.label
colnames(inc_plasmid_df) = 'Inc'

other_plasmid_df = data.frame(matrix(NA, nrow = length(tree$tip.label), ncol = 2))
colnames(other_plasmid_df) = c('Isolate','other_plasmid')
other_plasmid_df$Isolate = tree$tip.label
other_plasmid_df$other_plasmid = 'none'
for(i in 1:nrow(other_plasmid_df)){
  if('FIA' %in% plasmid_data$simplified_plasmid[plasmid_data$Sample == other_plasmid_df$Isolate[i]]){
    other_plasmid_df$other_plasmid[i] = 'FIA'
  }else if('p0111' %in% plasmid_data$simplified_plasmid[plasmid_data$Sample == other_plasmid_df$Isolate[i]]){
    other_plasmid_df$other_plasmid[i] = 'p0111'
  }else if('rep36' %in% plasmid_data$simplified_plasmid[plasmid_data$Sample == other_plasmid_df$Isolate[i]]){
    other_plasmid_df$other_plasmid[i] = 'rep36'
  }
}
other_plasmid_df = data.frame(other_plasmid_df$other_plasmid)
rownames(other_plasmid_df) = tree$tip.label
colnames(other_plasmid_df) = 'other_plasmid'

plasmids_per_isolate = data.frame(table(plasmid_data$Sample), stringsAsFactors = F)
colnames(plasmids_per_isolate) = c('Sample','All')
plasmids_per_isolate$genus = 'Klebsiella'
plasmids_per_isolate$genus[plasmids_per_isolate$Sample %in% plasmid_data$Sample[plasmid_data$genus == 'Escherichia']] = 'Escherichia'
plasmids_per_isolate$genus[plasmids_per_isolate$Sample %in% plasmid_data$Sample[plasmid_data$genus == 'Enterobacter']] = 'Enterobacter'
plasmids_per_isolate$status = 'Community'
plasmids_per_isolate$status[plasmids_per_isolate$Sample %in% pt_data$Sample[pt_data$Category == 'Environmental swab']] = 'Environmental swab'
plasmids_per_isolate$status[plasmids_per_isolate$Sample %in% pt_data$Sample[pt_data$Category == 'Healthcare-associated']] = 'Healthcare associated'
plasmids_per_isolate$status[plasmids_per_isolate$Sample %in% pt_data$Sample[pt_data$Category == 'Hospital-acquired']] = 'Hospital acquired'

## Preparing genus-specific plasmid dataframes
## Klebsiella
kleb_tree = drop.tip(kleb_tree, setdiff(kleb_tree$tip.label, pt_data$Sample))
kleb_col_df = data.frame(matrix(NA, nrow = length(kleb_tree$tip.label), ncol = 2), stringsAsFactors = F)
colnames(kleb_col_df) = c('Sample','Col')
kleb_col_df$Sample = kleb_tree$tip.label
for(i in 1:nrow(kleb_col_df)){
  kleb_col_df$Col[i] = col_plasmid_df$Col[rownames(col_plasmid_df) == kleb_col_df$Sample[i]]
}
rownames(kleb_col_df) = kleb_col_df$Sample
kleb_col_df = kleb_col_df[,'Col', drop = F]

kleb_inc_df = data.frame(matrix(NA, nrow = length(kleb_tree$tip.label), ncol = 2), stringsAsFactors = F)
colnames(kleb_inc_df) = c('Sample','Inc')
kleb_inc_df$Sample = kleb_tree$tip.label
for(i in 1:nrow(kleb_inc_df)){
  kleb_inc_df$Inc[i] = inc_plasmid_df$Inc[rownames(inc_plasmid_df) == kleb_inc_df$Sample[i]]
}
rownames(kleb_inc_df) = kleb_inc_df$Sample
kleb_inc_df = kleb_inc_df[,'Inc', drop = F]

kleb_otherplasmids_df = data.frame(matrix(NA, nrow = length(kleb_tree$tip.label), ncol = 2))
colnames(kleb_otherplasmids_df) = c('Sample','Other')
kleb_otherplasmids_df$Sample = kleb_tree$tip.label
for(i in 1:nrow(kleb_otherplasmids_df)){
  kleb_otherplasmids_df$Other[i] = other_plasmid_df$other_plasmid[rownames(other_plasmid_df) == kleb_otherplasmids_df$Sample[i]]
}
rownames(kleb_otherplasmids_df) = kleb_otherplasmids_df$Sample
kleb_otherplasmids_df = kleb_otherplasmids_df[,'Other', drop = F]

## Enterobacter
ent_col_df = data.frame(matrix(NA, nrow = length(ent_tree$tip.label), ncol = 2), stringsAsFactors = F)
colnames(ent_col_df) = c('Sample','Col')
ent_col_df$Sample = ent_tree$tip.label
for(i in 1:nrow(ent_col_df)){
  ent_col_df$Col[i] = col_plasmid_df$Col[rownames(col_plasmid_df) == ent_col_df$Sample[i]]
}
rownames(ent_col_df) = ent_col_df$Sample
ent_col_df = ent_col_df[,'Col', drop = F]

ent_inc_df = data.frame(matrix(NA, nrow = length(ent_tree$tip.label), ncol = 2), stringsAsFactors = F)
colnames(ent_inc_df) = c('Sample','Inc')
ent_inc_df$Sample = ent_tree$tip.label
for(i in 1:nrow(ent_inc_df)){
  ent_inc_df$Inc[i] = inc_plasmid_df$Inc[rownames(inc_plasmid_df) == ent_inc_df$Sample[i]]
}
rownames(ent_inc_df) = ent_inc_df$Sample
ent_inc_df = ent_inc_df[,'Inc', drop = F]

ent_otherplasmids_df = data.frame(matrix(NA, nrow = length(ent_tree$tip.label), ncol = 2))
colnames(ent_otherplasmids_df) = c('Sample','Other')
ent_otherplasmids_df$Sample = ent_tree$tip.label
for(i in 1:nrow(ent_otherplasmids_df)){
  ent_otherplasmids_df$Other[i] = other_plasmid_df$other_plasmid[rownames(other_plasmid_df) == ent_otherplasmids_df$Sample[i]]
}
rownames(ent_otherplasmids_df) = ent_otherplasmids_df$Sample
ent_otherplasmids_df = ent_otherplasmids_df[,'Other', drop = F]

## Escherichia
esch_col_df = data.frame(matrix(NA, nrow = length(esch_tree$tip.label), ncol = 2), stringsAsFactors = F)
colnames(esch_col_df) = c('Sample','Col')
esch_col_df$Sample = esch_tree$tip.label
for(i in 1:nrow(esch_col_df)){
  esch_col_df$Col[i] = col_plasmid_df$Col[rownames(col_plasmid_df) == esch_col_df$Sample[i]]
}
rownames(esch_col_df) = esch_col_df$Sample
esch_col_df = esch_col_df[,'Col', drop = F]

esch_inc_df = data.frame(matrix(NA, nrow = length(esch_tree$tip.label), ncol = 2), stringsAsFactors = F)
colnames(esch_inc_df) = c('Sample','Inc')
esch_inc_df$Sample = esch_tree$tip.label
for(i in 1:nrow(esch_inc_df)){
  esch_inc_df$Inc[i] = inc_plasmid_df$Inc[rownames(inc_plasmid_df) == esch_inc_df$Sample[i]]
}
rownames(esch_inc_df) = esch_inc_df$Sample
esch_inc_df = esch_inc_df[,'Inc', drop = F]

esch_otherplasmids_df = data.frame(matrix(NA, nrow = length(esch_tree$tip.label), ncol = 2))
colnames(esch_otherplasmids_df) = c('Sample','Other')
esch_otherplasmids_df$Sample = esch_tree$tip.label
for(i in 1:nrow(esch_otherplasmids_df)){
  esch_otherplasmids_df$Other[i] = other_plasmid_df$other_plasmid[rownames(other_plasmid_df) == esch_otherplasmids_df$Sample[i]]
}
rownames(esch_otherplasmids_df) = esch_otherplasmids_df$Sample
esch_otherplasmids_df = esch_otherplasmids_df[,'Other', drop = F]

# Reading in SNP information
# starting with Escherichia
esch_snps = read.table('./data/esch_pairwise_core_snps.tsv', sep = '\t')
colnames(esch_snps) = c('sample1','sample2','core_snps')
# removing 'pseudo reference' rows
esch_snps = esch_snps[esch_snps$sample2 != 'pseudo_reference_sequence' & esch_snps$sample1 != 'pseudo_reference_sequence',]

# Klebsiella
kleb_snps = read.table('./data/kleb_pairwise_core_snps.tsv', sep = '\t')
colnames(kleb_snps) = c('sample1','sample2','core_snps')
# removing 'pseudo reference' rows
kleb_snps = kleb_snps[kleb_snps$sample2 != 'pseudo_reference_sequence' & kleb_snps$sample1 != 'pseudo_reference_sequence',]

external_mash = read.csv('./data/241213_CRESCENT_hq_external_genomes_mash_checkm_targettaxa.csv', stringsAsFactors = F)

#### Plotting - all CRESCENT isolates ####
## plotting with Mash species identification call on tips
# setting mash colors
mash_colors = c(colors$color[colors$topic == 'Mash'])
names(mash_colors) = c(colors$name[colors$topic == 'Mash'])

mash_df = mash_df[mash_df$Sample %in% tree$tip.label,]
treeplot = ggtree(tree) %<+% mash_df
# generating mlst-tip treeplot
tip_treeplot = treeplot + geom_tiplab(align = TRUE, size = 0)+ geom_tippoint(aes(color = Mash), size = 2.3) +
  scale_color_manual(values = mash_colors, breaks = names(mash_colors))+
  geom_treescale(x = 0.1, y = 20)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'))

tip_treeplot

# plotting epi status + Mash species identification
haca_colors = c(colors$color[colors$topic == 'HA/CA status'])
names(haca_colors) = c(colors$name[colors$topic == 'HA/CA status'])
epi_status_tree = gheatmap(tip_treeplot, epi_df, colnames_angle = 85, colnames_offset_y = -15,
                           offset = 0, width = 0.025, font.size = 5)+
  scale_fill_manual(values = haca_colors, breaks = names(haca_colors))+ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Epidemiological status')

#pdf(paste0('./figures/',Sys.Date(),'_CRESCENT_figure2A.pdf'), width = 15, height = 17)
epi_status_tree
#dev.off()

#### plots - initial escherichia isolate trees ####
# subsetting mash data to escherichia-only
esch_data = pt_data[grep('Escherichia', pt_data$WGS.identification),]
esch_tree = drop.tip(esch_tree, setdiff(esch_tree$tip.label, pt_data$Sample))
# making mash_df ordered by tree tip labels
esch_df = data.frame(matrix(NA, nrow = length(esch_tree$tip.label), ncol = 2), stringsAsFactors = F)
colnames(esch_df) = c('Sample','Mash')
esch_df$Sample = esch_tree$tip.label
for(i in 1:length(esch_tree$tip.label)){
  esch_df$Mash[i] = esch_data$WGS.identification[esch_data$Sample == esch_df$Sample[i]]
}

rownames(esch_df) = esch_data$Sample
# subsetting tree
esch_tree = drop.tip(esch_tree,setdiff(esch_tree$tip.label, pt_data$Sample))
# rooting tree
esch_tree = midpoint.root(esch_tree)
# plotting tree
esch_treeplot = ggtree(esch_tree)
circle_esch_treeplot = ggtree(esch_tree, layout = 'circular')

# making pt_dataframe in order of tree$tiplabels
esch_pt_df = data.frame(matrix(NA, nrow = length(esch_tree$tip.label), ncol = 1), stringsAsFactors = F)
colnames(esch_pt_df) = 'Sample'
esch_pt_df$Sample = esch_tree$tip.label
for(i in 1:length(esch_tree$tip.label)){
  esch_pt_df$location_swab[i] = esch_data$location_swab[esch_data$WGS_ID == esch_pt_df$Sample[i]]
  esch_pt_df$location_bed[i] = esch_data$location_bed[esch_data$WGS_ID == esch_pt_df$Sample[i]]
  esch_pt_df$hospital_days_before_swab[i] = esch_data$hospital_days_before_swab[esch_data$WGS_ID == esch_pt_df$Sample[i]]
  esch_pt_df$swab_day[i] = esch_data$swab_day[esch_data$WGS_ID == esch_pt_df$Sample[i]]
  esch_pt_df$icu_los[i] = esch_data$icu_los[esch_data$WGS_ID == esch_pt_df$Sample[i]]
  esch_pt_df$total_los[i] = esch_data$total_los[esch_data$WGS_ID == esch_pt_df$Sample[i]]
  esch_pt_df$isolation_source[i] = esch_data$isolation_source[esch_data$WGS_ID == esch_pt_df$Sample[i]]
  esch_pt_df$hosp_last_1yr[i] = esch_data$hosp_last_1yr[esch_data$WGS_ID == esch_pt_df$Sample[i]]
  esch_pt_df$curr_abx_exp[i] = esch_data$curr_abx_exp[esch_data$WGS_ID == esch_pt_df$Sample[i]]
}
# ordering epi status
esch_epi_df = data.frame(matrix(NA, nrow = length(esch_tree$tip.label), ncol = 2), stringsAsFactors = F)
colnames(esch_epi_df) = c('Sample','Epidemiolgical status')
esch_epi_df$Sample = esch_tree$tip.label
for(i in 1:length(esch_tree$tip.label)){
  esch_epi_df$`Epidemiolgical status`[i] = esch_data$Category[esch_data$Sample == esch_epi_df$Sample[i]]
}
esch_epi_df = data.frame(esch_epi_df$`Epidemiolgical status`)
colnames(esch_epi_df) = 'Epidemiological status'
rownames(esch_epi_df) = esch_tree$tip.label


## plotting with mash call on tips
# setting mash colors
esch_mash_colors = mash_colors
names(esch_mash_colors) = names(mash_colors)

## plotting with MLST tip annotation
# adding mlst to circle treeplot
circle_esch_treeplot = ggtree(esch_tree, layout = 'circular', open.angle = T) %<+% esch_df
# generating mlst-tip treeplot
tip_circle_esch_treeplot = circle_esch_treeplot + geom_tiplab(align = TRUE, size = 0)+geom_tippoint(aes(color = Mash), size = 2.2) +
  scale_color_manual(values = esch_mash_colors, breaks = names(esch_mash_colors))+
  geom_treescale(x = 0.0011, y = 10, offset = 0)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'))
tip_circle_esch_treeplot

# plotting epi status
# adding HAHOCA ring
esch_epi_status_tree_fan = gheatmap(tip_circle_esch_treeplot, esch_epi_df, colnames_angle = 85, colnames_offset_y = -9,
                                    offset = 0, width = 0.04, font.size = 5)+
  scale_fill_manual(values = haca_colors, breaks = names(haca_colors))+ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Epidemiological status')

#pdf(paste0('./figures/',Sys.Date(),'_crescent_epistatus_mash_escherichia_internal_fan.pdf'), width = 15)
esch_epi_status_tree_fan
#dev.off()

esch_epi_status_tree_fan2 = esch_epi_status_tree_fan + new_scale_fill()

#### plots - NCBI + CRESCENT isolates, Escherichia - Supplemental Figure 5A ####
ext_esch_samples = c(pt_data$Sample[pt_data$genus == 'Escherichia'], external_mash$Sample)
ext_esch_tree = drop.tip(ext_esch_tree, setdiff(ext_esch_tree$tip.label, ext_esch_samples))
# making mash_df ordered by tree tip labels
ext_esch_df = data.frame(matrix(NA, nrow = length(ext_esch_tree$tip.label), ncol = 3), stringsAsFactors = F)
colnames(ext_esch_df) = c('Sample','Mash','Category')
ext_esch_df$Sample = ext_esch_tree$tip.label
for(i in 1:length(ext_esch_tree$tip.label)){
  if(ext_esch_df$Sample[i] %in% pt_data$Sample){
  ext_esch_df$Mash[i] = pt_data$WGS.identification[pt_data$Sample == ext_esch_df$Sample[i]]
  ext_esch_df$Category[i] = pt_data$Category[pt_data$Sample == ext_esch_df$Sample[i]]
  }else{
    ext_esch_df$Mash[i] = external_mash$...7[external_mash$Sample == ext_esch_df$Sample[i]]
    ext_esch_df$Category[i] = 'External isolate'
  }
  }

rownames(ext_esch_df) = ext_esch_df$Sample
# rooting tree
ext_esch_tree = midpoint.root(ext_esch_tree)
# adding Mash to circle treeplot
circle_ext_esch_treeplot = ggtree(ext_esch_tree, layout = 'circular', open.angle = T) %<+% ext_esch_df
# generating mlst-tip treeplot
tip_circle_ext_esch_treeplot = circle_ext_esch_treeplot + geom_tiplab(align = TRUE, size = 0)+geom_tippoint(aes(color = Mash), size = 2.2) +
  scale_color_manual(values = esch_mash_colors, breaks = names(esch_mash_colors))+
  geom_treescale(x = 0.0011, y = 10, offset = 0)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'))
tip_circle_ext_esch_treeplot

# plotting epi status
# adding HAHOCA ring
ext_esch_epi_df = ext_esch_df[,'Category',drop = F]
rownames(ext_esch_epi_df) = ext_esch_df$Sample
ext_esch_epi_status_tree_fan = gheatmap(tip_circle_ext_esch_treeplot, ext_esch_epi_df, colnames_angle = 85, colnames_offset_y = -9,
                                    offset = 0, width = 0.04, font.size = 5)+
  scale_fill_manual(values = haca_colors, breaks = names(haca_colors))+ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Epidemiological status')

#pdf(paste0('./figures/',Sys.Date(),'_crescent_NCBI_epistatus_mash_escherichia_fan.pdf'), width = 15)
ext_esch_epi_status_tree_fan
#dev.off()
#### plots - initial enterobacter isolate ####
# subsetting mash data to ent-only
ent_mash = pt_data[grep('Enterobacter', pt_data$WGS.identification),]
ent_tree = drop.tip(ent_tree, setdiff(ent_tree$tip.label, pt_data$Sample))

# making mash_df ordered by tree tip labels
ent_mash_df = data.frame(matrix(NA, nrow = length(ent_tree$tip.label), ncol = 2), stringsAsFactors = F)
colnames(ent_mash_df) = c('Sample','Mash')
ent_mash_df$Sample = ent_tree$tip.label
for(i in 1:length(ent_tree$tip.label)){
  ent_mash_df$Mash[i] = ent_mash$WGS.identification[ent_mash$Sample == ent_mash_df$Sample[i]]
}

rownames(ent_mash_df) = ent_mash_df$Sample
# rooting tree
ent_tree = midpoint.root(ent_tree)
# plotting tree
circle_ent_treeplot = ggtree(ent_tree, layout = 'circular')
# subsetting pt_data
ent_pt_data = pt_data[pt_data$Sample %in% ent_tree$tip.label,]

# making pt_dataframe in order of tree$tiplabels
ent_pt_df = data.frame(matrix(NA, nrow = length(ent_tree$tip.label), ncol = 1), stringsAsFactors = F)
colnames(ent_pt_df) = 'Sample'
ent_pt_df$Sample = ent_tree$tip.label
for(i in 1:length(ent_tree$tip.label)){
  ent_pt_df$location_swab[i] = ent_pt_data$location_swab[ent_pt_data$WGS_ID == ent_pt_df$Sample[i]]
  ent_pt_df$location_bed[i] = ent_pt_data$location_bed[ent_pt_data$WGS_ID == ent_pt_df$Sample[i]]
  ent_pt_df$swab_day[i] = ent_pt_data$swab_coll_date[ent_pt_data$WGS_ID == ent_pt_df$Sample[i]]
  ent_pt_df$hosp_last_1yr[i] = ent_pt_data$hosp_last_1yr[ent_pt_data$WGS_ID == ent_pt_df$Sample[i]]
  ent_pt_df$curr_abx_exp[i] = ent_pt_data$curr_abx_exp[ent_pt_data$WGS_ID == ent_pt_df$Sample[i]]
}
## plotting with mash call on tips
# setting mash colors
ent_mash_colors = mash_colors
names(ent_mash_colors) = names(mash_colors)

## plotting with Mash identification tip annotation
circle_ent_treeplot = ggtree(ent_tree, layout = 'circular', open.angle = T) %<+% ent_mash_df
# generating mlst-tip treeplot
tip_circle_ent_treeplot = circle_ent_treeplot + geom_tiplab(align = TRUE, size = 0)+geom_tippoint(aes(color = Mash), size = 2.2) +
  geom_treescale()+
  scale_color_manual(values = ent_mash_colors, breaks = names(ent_mash_colors))

# ordering epi status
ent_epi_df = data.frame(matrix(NA, nrow = length(ent_tree$tip.label), ncol = 2), stringsAsFactors = F)
colnames(ent_epi_df) = c('Sample','Epidemiolgical status')
ent_epi_df$Sample = ent_tree$tip.label
for(i in 1:length(ent_tree$tip.label)){
  ent_epi_df$`Epidemiolgical status`[i] = ent_pt_data$Category[ent_pt_data$Sample == ent_epi_df$Sample[i]]
}
ent_epi_df = data.frame(ent_epi_df$`Epidemiolgical status`)
colnames(ent_epi_df) = 'Epidemiological status'
rownames(ent_epi_df) = ent_tree$tip.label

ent_epi_tree_fan = gheatmap(tip_circle_ent_treeplot, ent_epi_df, colnames_angle = 85, colnames_offset_y = -3,
                            offset = 0, width = 0.04, font.size = 5)+
  scale_fill_manual(values = haca_colors, breaks = names(haca_colors))+ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Epidemiological status')

#pdf(paste0('./figures/',Sys.Date(),'_crescent_epistatus_mash_ent_internal_fan.pdf'), width = 15, height = 15)
ent_epi_tree_fan
#dev.off()

ent_epi_tree_fan2 = ent_epi_tree_fan + new_scale_fill()

#### plots - NCBI + CRESCENT isolates, Enterobacter - Supplemental Figure 5C ####
ext_ent_samples = c(pt_data$Sample[pt_data$genus == 'Enterobacter'], external_mash$Sample)
ext_ent_tree = drop.tip(ext_ent_tree, setdiff(ext_ent_tree$tip.label, ext_ent_samples))
# making mash_df ordered by tree tip labels
ext_ent_df = data.frame(matrix(NA, nrow = length(ext_ent_tree$tip.label), ncol = 3), stringsAsFactors = F)
colnames(ext_ent_df) = c('Sample','Mash','Category')
ext_ent_df$Sample = ext_ent_tree$tip.label
for(i in 1:length(ext_ent_tree$tip.label)){
  if(ext_ent_df$Sample[i] %in% pt_data$Sample){
    ext_ent_df$Mash[i] = pt_data$WGS.identification[pt_data$Sample == ext_ent_df$Sample[i]]
    ext_ent_df$Category[i] = pt_data$Category[pt_data$Sample == ext_ent_df$Sample[i]]
  }else{
    ext_ent_df$Mash[i] = external_mash$...7[external_mash$Sample == ext_ent_df$Sample[i]]
    ext_ent_df$Category[i] = 'External isolate'
  }
}

rownames(ext_ent_df) = ext_ent_df$Sample
# rooting tree
ext_ent_tree = midpoint.root(ext_ent_tree)
# adding Mash to circle treeplot
circle_ext_ent_treeplot = ggtree(ext_ent_tree, layout = 'circular', open.angle = T) %<+% ext_ent_df
# generating mlst-tip treeplot
tip_circle_ext_ent_treeplot = circle_ext_ent_treeplot + geom_tiplab(align = TRUE, size = 0)+geom_tippoint(aes(color = Mash), size = 2.2) +
  scale_color_manual(values = esch_mash_colors, breaks = names(esch_mash_colors))+
  geom_treescale(x = 0.0011, y = 10, offset = 0)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'))
tip_circle_ext_ent_treeplot

# plotting epi status
# adding HAHOCA ring
ext_ent_epi_df = ext_ent_df[,'Category',drop = F]
rownames(ext_ent_epi_df) = ext_ent_df$Sample
ext_ent_epi_status_tree_fan = gheatmap(tip_circle_ext_ent_treeplot, ext_ent_epi_df, colnames_angle = 85, colnames_offset_y = -9,
                                        offset = 0, width = 0.04, font.size = 5)+
  scale_fill_manual(values = haca_colors, breaks = names(haca_colors))+ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Epidemiological status')

#pdf(paste0('./figures/',Sys.Date(),'_crescent_NCBI_epistatus_mash_enterobacter_fan.pdf'), width = 15)
ext_ent_epi_status_tree_fan
#dev.off()
#### plots - initial klebsiella isolates ####
# subsetting mash data to kleb-only
kleb_mash = pt_data[grep('Klebsiella', pt_data$WGS.identification),]
kleb_tree = drop.tip(kleb_tree, setdiff(kleb_tree$tip.label, pt_data$Sample))

# making mash_df ordered by tree tip labels
kleb_mash_df = data.frame(matrix(NA, nrow = length(kleb_tree$tip.label), ncol = 2), stringsAsFactors = F)
colnames(kleb_mash_df) = c('Sample','Mash')
kleb_mash_df$Sample = kleb_tree$tip.label
for(i in 1:length(kleb_tree$tip.label)){
  kleb_mash_df$Mash[i] = kleb_mash$WGS.identification[kleb_mash$Sample == kleb_mash_df$Sample[i]]
}

rownames(kleb_mash_df) = kleb_mash_df$Sample
# subsetting tree
kleb_tree = drop.tip(kleb_tree,setdiff(kleb_tree$tip.label, pt_data$Sample))
# rooting tree
kleb_tree = midpoint.root(kleb_tree)
# plotting tree
kleb_treeplot = ggtree(kleb_tree)
circle_kleb_treeplot = ggtree(kleb_tree, layout = 'circular')
# subsetting pt_data
kleb_pt_data = pt_data[pt_data$Sample %in% kleb_tree$tip.label,]

# ordering epi status
kleb_epi_df = data.frame(matrix(NA, nrow = length(kleb_tree$tip.label), ncol = 2), stringsAsFactors = F)
colnames(kleb_epi_df) = c('Sample','Epidemiolgical status')
kleb_epi_df$`Epidemiolgical status` = 'External isolate'
kleb_epi_df$Sample = kleb_tree$tip.label
for(i in 1:length(kleb_tree$tip.label)){
  if(kleb_epi_df$Sample[i] %in% kleb_pt_data$Sample){
    kleb_epi_df$`Epidemiolgical status`[i] = kleb_pt_data$Category[kleb_pt_data$Sample == kleb_epi_df$Sample[i]]
  }
}

kleb_epi_df = data.frame(kleb_epi_df$`Epidemiolgical status`)
colnames(kleb_epi_df) = 'Epidemiological status'
rownames(kleb_epi_df) = kleb_tree$tip.label

## plotting with mash call on tips
# setting mash colors
kleb_mash_colors = mash_colors
names(kleb_mash_colors) = names(mash_colors)
circle_kleb_treeplot = ggtree(kleb_tree, layout = 'circular', open.angle = T) %<+% kleb_mash_df
tip_circle_kleb_treeplot = circle_kleb_treeplot + geom_tiplab(align = TRUE, size = 0)+geom_tippoint(aes(color = Mash), size = 2.2) +
  geom_treescale()+
  scale_color_manual(values = kleb_mash_colors, breaks = names(kleb_mash_colors))
tip_circle_kleb_treeplot

# plotting epi status
kleb_epi_status_tree_fan = gheatmap(tip_circle_kleb_treeplot, kleb_epi_df, colnames_angle = 85, colnames_offset_y = -9,
                                    offset = 0.003, width = 0.06, font.size = 5)+
  scale_fill_manual(values = haca_colors, breaks = names(haca_colors))+ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Epidemiological status')

#pdf(paste0('./figures/',Sys.Date(),'_crescent_epistatus_mash_klebsiella_internal_fan.pdf'), width = 15, height = 15)
kleb_epi_status_tree_fan
#dev.off()

kleb_epi_status_tree_fan2 = kleb_epi_status_tree_fan + new_scale_fill()

#### plots - NCBI + CRESCENT isolates, Klebsiella - Supplemental Figure 5B ####
ext_kleb_samples = c(pt_data$Sample[pt_data$genus == 'Klebsiella'], external_mash$Sample)
ext_kleb_tree = drop.tip(ext_kleb_tree, setdiff(ext_kleb_tree$tip.label, ext_kleb_samples))
# making mash_df ordered by tree tip labels
ext_kleb_df = data.frame(matrix(NA, nrow = length(ext_kleb_tree$tip.label), ncol = 3), stringsAsFactors = F)
colnames(ext_kleb_df) = c('Sample','Mash','Category')
ext_kleb_df$Sample = ext_kleb_tree$tip.label
for(i in 1:length(ext_kleb_tree$tip.label)){
  if(ext_kleb_df$Sample[i] %in% pt_data$Sample){
    ext_kleb_df$Mash[i] = pt_data$WGS.identification[pt_data$Sample == ext_kleb_df$Sample[i]]
    ext_kleb_df$Category[i] = pt_data$Category[pt_data$Sample == ext_kleb_df$Sample[i]]
  }else{
    ext_kleb_df$Mash[i] = external_mash$...7[external_mash$Sample == ext_kleb_df$Sample[i]]
    ext_kleb_df$Category[i] = 'External isolate'
  }
}

rownames(ext_kleb_df) = ext_kleb_df$Sample
# rooting tree
ext_kleb_tree = midpoint.root(ext_kleb_tree)
# adding Mash to circle treeplot
circle_ext_kleb_treeplot = ggtree(ext_kleb_tree, layout = 'circular', open.angle = T) %<+% ext_kleb_df
# generating mlst-tip treeplot
tip_circle_ext_kleb_treeplot = circle_ext_kleb_treeplot + geom_tiplab(align = TRUE, size = 0)+geom_tippoint(aes(color = Mash), size = 2.2) +
  scale_color_manual(values = esch_mash_colors, breaks = names(esch_mash_colors))+
  geom_treescale(x = 0.0011, y = 10, offset = 0)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'))
tip_circle_ext_kleb_treeplot

# plotting epi status
# adding HAHOCA ring
ext_kleb_epi_df = ext_kleb_df[,'Category',drop = F]
rownames(ext_kleb_epi_df) = ext_kleb_df$Sample
ext_kleb_epi_status_tree_fan = gheatmap(tip_circle_ext_kleb_treeplot, ext_kleb_epi_df, colnames_angle = 85, colnames_offset_y = -9,
                                        offset = 0, width = 0.04, font.size = 5)+
  scale_fill_manual(values = haca_colors, breaks = names(haca_colors))+ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Epidemiological status')

#pdf(paste0('./figures/',Sys.Date(),'_crescent_NCBI_epistatus_mash_escherichia_fan.pdf'), width = 15)
ext_kleb_epi_status_tree_fan
#dev.off()
#### adding carbapenemases to Enterobacter, Klebsiella, and Escherichia trees ####
names(carb_colors)[length(carb_colors)] = ''
carb_colors = colors$color[colors$topic == 'Carbapenemase']
names(carb_colors) = colors$name[colors$topic == 'Carbapenemase']

## want to plot x-axis as epi category, y-axis as carbapenemase genes, points sized by count of isolates
carbapenemase_count_df = data.frame(matrix(NA, nrow = 4*length(table(amrfinder_data$Gene.symbol[amrfinder_data$Subclass == 'CARBAPENEM'])),
                                           ncol = 6), stringsAsFactors = F)
colnames(carbapenemase_count_df) = c('Carbapenemase','Status','Freq','EcFreq','EntFreq','KlebFreq')
carbapenemase_count_df$Carbapenemase = rep(names(table(amrfinder_data$Gene.symbol[amrfinder_data$Subclass == 'CARBAPENEM'])),4)
carbapenemase_count_df$Status = c(rep('Hospital acquired',7),rep('Healthcare associated',7),rep('Community',7),
                                  rep('Environmental swab',7))
for(i in 1:nrow(carbapenemase_count_df)){
  carbapenemase_count_df$EcFreq[i] = sum(amrfinder_data$Gene.symbol == carbapenemase_count_df$Carbapenemase[i] & amrfinder_data$Category == carbapenemase_count_df$Status[i] &
                                           amrfinder_data$genus == 'Escherichia')
  carbapenemase_count_df$KlebFreq[i] = sum(amrfinder_data$Gene.symbol == carbapenemase_count_df$Carbapenemase[i] & amrfinder_data$Category == carbapenemase_count_df$Status[i] &
                                             amrfinder_data$genus == 'Klebsiella')
  carbapenemase_count_df$EntFreq[i] = sum(amrfinder_data$Gene.symbol == carbapenemase_count_df$Carbapenemase[i] & amrfinder_data$Category == carbapenemase_count_df$Status[i] & 
                                            amrfinder_data$genus == 'Enterobacter')
  carbapenemase_count_df$Freq[i] = sum(amrfinder_data$Gene.symbol == carbapenemase_count_df$Carbapenemase[i] & amrfinder_data$Category == carbapenemase_count_df$Status[i])
}

ndmoxa_df = data.frame(matrix(NA, nrow = 12, ncol = 3), stringsAsFactors = F)
colnames(ndmoxa_df) = c('blaNDM','blaOXA','Freq')
ndmoxa_df$blaOXA = c(rep('blaOXA-181',4), rep('blaOXA-232',4), rep('blaOXA-484',4))
ndmoxa_df$blaNDM = rep(c('blaNDM','blaNDM-1','blaNDM-5','blaNDM-7'),3)
for(i in 1:nrow(ndmoxa_df)){
  ndmoxa_df$Freq[i] = sum(carb_data[colnames(carb_data) == ndmoxa_df$blaNDM[i]] == ndmoxa_df$blaNDM[i] & carb_data[colnames(carb_data) == ndmoxa_df$blaOXA[i]] == ndmoxa_df$blaOXA[i])
}

## Enterobacter trees
ent_carbdf = carb_df[rownames(carb_df) %in% ent_tree$tip.label,]
ordered_entcarbdf = data.frame(matrix(NA, nrow = length(ent_tree$tip.label), ncol = 2))
colnames(ordered_entcarbdf) = c('blaNDM','blaOXA')
rownames(ordered_entcarbdf) = ent_tree$tip.label
for(i in 1:nrow(ordered_entcarbdf)){
  ordered_entcarbdf$blaNDM[i] = ent_carbdf$blaNDM[rownames(ent_carbdf) == rownames(ordered_entcarbdf)[i]]
  ordered_entcarbdf$blaOXA[i] = ent_carbdf$blaOXA[rownames(ent_carbdf) == rownames(ordered_entcarbdf)[i]]
}

ent_epi_carb_tree_fan = gheatmap(ent_epi_tree_fan2, ordered_entcarbdf, colnames_angle = 85, colnames_offset_y = -1,
                             offset = 0.009, width = 0.06, font.size = 5)+
  scale_fill_manual(values = carb_colors, breaks = names(carb_colors))+ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Carbapenemase')


#pdf(paste0('./figures/', Sys.Date(),'_crescent_enttree_epistatus_carbapenemase_fan.pdf'),width = 15, height = 15)
ent_epi_carb_tree_fan
#dev.off()

ent_epi_carb_tree_fan1 = ent_epi_carb_tree_fan + new_scale_fill()

## Escherichia
esch_carbdf = carb_df[rownames(carb_df) %in% esch_tree$tip.label,]
ordered_eschcarbdf = data.frame(matrix(NA, nrow = length(esch_tree$tip.label), ncol = 2))
colnames(ordered_eschcarbdf) = c('blaNDM','blaOXA')
rownames(ordered_eschcarbdf) = esch_tree$tip.label
for(i in 1:nrow(ordered_eschcarbdf)){
  ordered_eschcarbdf$blaNDM[i] = esch_carbdf$blaNDM[rownames(esch_carbdf) == rownames(ordered_eschcarbdf)[i]]
  ordered_eschcarbdf$blaOXA[i] = esch_carbdf$blaOXA[rownames(esch_carbdf) == rownames(ordered_eschcarbdf)[i]]
}

esch_epi_carb_tree_fan = gheatmap(esch_epi_status_tree_fan2, ordered_eschcarbdf, colnames_angle = 85, colnames_offset_y = -1,
                                  offset = 0.0015, width = 0.09, font.size = 5)+
  scale_fill_manual(values = carb_colors, breaks = names(carb_colors))+ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Carbapenemase')


#pdf(paste0('./figures/', Sys.Date(),'_crescent_ecolitree_epistatus_carbapenemase_fan.pdf'),width = 15, height = 15)
esch_epi_carb_tree_fan
#dev.off()

esch_epi_carb_tree_fan1 = esch_epi_carb_tree_fan + new_scale_fill()

## Klebsiella
kleb_carbdf = carb_df[rownames(carb_df) %in% kleb_tree$tip.label,]
ordered_klebcarbdf = data.frame(matrix(NA, nrow = length(kleb_tree$tip.label), ncol = 2))
colnames(ordered_klebcarbdf) = c('blaNDM','blaOXA')
rownames(ordered_klebcarbdf) = kleb_tree$tip.label
for(i in 1:nrow(ordered_klebcarbdf)){
  ordered_klebcarbdf$blaNDM[i] = kleb_carbdf$blaNDM[rownames(kleb_carbdf) == rownames(ordered_klebcarbdf)[i]]
  ordered_klebcarbdf$blaOXA[i] = kleb_carbdf$blaOXA[rownames(kleb_carbdf) == rownames(ordered_klebcarbdf)[i]]
}

kleb_epi_carb_tree_fan = gheatmap(kleb_epi_status_tree_fan2, ordered_klebcarbdf, colnames_angle = 85, colnames_offset_y = -1,
                                  offset = 0.007, width = 0.09, font.size = 5)+
  scale_fill_manual(values = carb_colors, breaks = names(carb_colors))+ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Carbapenemase')


#pdf(paste0('./figures/', Sys.Date(),'_crescent_klebtree_epistatus_carbapenemase_fan.pdf'),width = 15, height = 15)
kleb_epi_carb_tree_fan
#dev.off()

kleb_epi_carb_tree_fan1 = kleb_epi_carb_tree_fan + new_scale_fill()

#### carbapenemase figures - Supplemental Figure 4A, Figure 3A ####
## want to plot x-axis as epi category, y-axis as carbapenemase genes, points sized by count of isolates
carbapenemase_count_df = data.frame(matrix(NA, nrow = 4*length(table(amrfinder_data$Gene.symbol[amrfinder_data$Subclass == 'CARBAPENEM'])),
                                           ncol = 6), stringsAsFactors = F)
colnames(carbapenemase_count_df) = c('Carbapenemase','Status','Freq','EcFreq','EntFreq','KlebFreq')
carbapenemase_count_df$Carbapenemase = rep(names(table(amrfinder_data$Gene.symbol[amrfinder_data$Subclass == 'CARBAPENEM'])),4)
carbapenemase_count_df$Status = c(rep('Hospital acquired',7),rep('Healthcare associated',7),rep('Community',7),
                                  rep('Environmental swab',7))
for(i in 1:nrow(carbapenemase_count_df)){
  carbapenemase_count_df$EcFreq[i] = sum(amrfinder_data$Gene.symbol == carbapenemase_count_df$Carbapenemase[i] & amrfinder_data$Category == carbapenemase_count_df$Status[i] &
                                           amrfinder_data$genus == 'Escherichia')
  carbapenemase_count_df$KlebFreq[i] = sum(amrfinder_data$Gene.symbol == carbapenemase_count_df$Carbapenemase[i] & amrfinder_data$Category == carbapenemase_count_df$Status[i] &
                                             amrfinder_data$genus == 'Klebsiella')
  carbapenemase_count_df$EntFreq[i] = sum(amrfinder_data$Gene.symbol == carbapenemase_count_df$Carbapenemase[i] & amrfinder_data$Category == carbapenemase_count_df$Status[i] & 
                                            amrfinder_data$genus == 'Enterobacter')
  carbapenemase_count_df$Freq[i] = sum(amrfinder_data$Gene.symbol == carbapenemase_count_df$Carbapenemase[i] & amrfinder_data$Category == carbapenemase_count_df$Status[i])
}

ndmoxa_df = data.frame(matrix(NA, nrow = 12, ncol = 3), stringsAsFactors = F)
colnames(ndmoxa_df) = c('blaNDM','blaOXA','Freq')
ndmoxa_df$blaOXA = c(rep('blaOXA-181',4), rep('blaOXA-232',4), rep('blaOXA-484',4))
ndmoxa_df$blaNDM = rep(c('blaNDM','blaNDM-1','blaNDM-5','blaNDM-7'),3)
for(i in 1:nrow(ndmoxa_df)){
  ndmoxa_df$Freq[i] = sum(carb_data[colnames(carb_data) == ndmoxa_df$blaNDM[i]] == ndmoxa_df$blaNDM[i] & carb_data[colnames(carb_data) == ndmoxa_df$blaOXA[i]] == ndmoxa_df$blaOXA[i])
}

## Supplemental Figure 4A
#pdf(paste0('./figures/', Sys.Date(),'_crescent_allisolates_blaNDM_blaOXA_bubbleplot.pdf'), height = 8, width = 10)
ggplot(ndmoxa_df, aes(x = factor(blaOXA, levels = c('blaOXA-484','blaOXA-232','blaOXA-181')),
                      y = factor(blaNDM, levels = c('blaNDM-7','blaNDM','blaNDM-1','blaNDM-5')), 
                      size = Freq, color = blaOXA))+
  geom_point(alpha = 0.95)+
  scale_size_area(max_size = 30)+
  geom_text(aes(label = Freq),  # Label by gene  # Place text in the middle of each segment
            color = "white", size = 4, fontface = 'bold')+  
  labs(y = 'blaNDM', x = 'blaOXA', size = 'Number of isolates')+
  theme_prism()+
  scale_color_manual(values = natparks.pals('Arches',3))+
  theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 0.9, hjust = 0.9), axis.text.y = element_text(size = 14),
        axis.ticks = element_blank(),
        title = element_text(size = 12, face = 'bold'), legend.text = element_text(size = 12))+
  ggtitle('Carbapenemases of all Crescent isolates')
#dev.off()
carbapenemase_count_df$Status[carbapenemase_count_df$Status == 'Community'] = 'CA'
carbapenemase_count_df$Status[carbapenemase_count_df$Status == 'Environmental swab'] = 'E'
carbapenemase_count_df$Status[carbapenemase_count_df$Status == 'Healthcare associated'] = 'HCA'
carbapenemase_count_df$Status[carbapenemase_count_df$Status == 'Hospital acquired'] = 'HA'

## Figure 3A
#pdf(paste0('./figures/', Sys.Date(),'_crescent_allisolates_carbapenemase_ndmvsoxa_bubbleplot.pdf'), height = 8, width = 10)
ggplot(carbapenemase_count_df, aes(x = factor(Status, levels = c('CA','E','HCA','HA')), 
                                   y = factor(Carbapenemase, levels = c('blaNDM-7','blaNDM','blaOXA-484','blaOXA-181','blaOXA-232',
                                                                        'blaNDM-1','blaNDM-5')), 
                                   size = Freq, color = Carbapenemase))+
  geom_point(alpha = 0.85)+
  scale_size_area(max_size = 30)+
  geom_text(aes(label = Freq),  # Label by gene  # Place text in the middle of each segment
            color = "white", size = 4, fontface = 'bold')+  
  labs(y = 'Carbapenemase', x = 'Epidemiological status', size = 'Number of isolates')+
  theme_bw()+
  scale_color_manual(values = carb_colors, breaks = names(carb_colors))+
  theme(axis.text.x = element_text(size = 14, vjust = 0.9, hjust = 0.9), axis.text.y = element_text(size = 14),
        axis.ticks = element_blank(),
        title = element_text(size = 12, face = 'bold'), legend.text = element_text(size = 12))+
  ggtitle('Carbapenemases of all Crescent isolates')
#dev.off()


carbapenemase_count_df = data.frame(matrix(NA, nrow = 4*length(table(amrfinder_data$Gene.symbol[amrfinder_data$Subclass == 'CARBAPENEM'])),
                                           ncol = 6), stringsAsFactors = F)
colnames(carbapenemase_count_df) = c('Carbapenemase','Status','Freq','EcFreq','EntFreq','KlebFreq')
carbapenemase_count_df$Carbapenemase = rep(names(table(amrfinder_data$Gene.symbol[amrfinder_data$Subclass == 'CARBAPENEM'])),4)
carbapenemase_count_df$Status = c(rep('Hospital acquired',7),rep('Healthcare associated',7),rep('Community',7),
                                  rep('Environmental swab',7))
for(i in 1:nrow(carbapenemase_count_df)){
  carbapenemase_count_df$Freq[i] = sum(amrfinder_data$Gene.symbol == carbapenemase_count_df$Carbapenemase[i] & amrfinder_data$Category == carbapenemase_count_df$Status[i])
  carbapenemase_count_df$EcFreq[i] = sum(amrfinder_data$Gene.symbol == carbapenemase_count_df$Carbapenemase[i] & amrfinder_data$Category == carbapenemase_count_df$Status[i] &
                                           amrfinder_data$genus == 'Escherichia')
  carbapenemase_count_df$KlebFreq[i] = sum(amrfinder_data$Gene.symbol == carbapenemase_count_df$Carbapenemase[i] & amrfinder_data$Category == carbapenemase_count_df$Status[i] &
                                             amrfinder_data$genus == 'Klebsiella')
  carbapenemase_count_df$EntFreq[i] = sum(amrfinder_data$Gene.symbol == carbapenemase_count_df$Carbapenemase[i] & amrfinder_data$Category == carbapenemase_count_df$Status[i] & 
                                            amrfinder_data$genus == 'Enterobacter')
}


#### adding month to full tree - Figure 2A ####
daytomonth = read.csv('./data/crescent_day_to_month.csv', stringsAsFactors = F)
study_day_data = data.frame(read_excel('./data/240822_CRESCENT_metadata_eeb.xlsx'), stringsAsFactors = F)
study_day_data$WGS_ID = gsub('K','K-',study_day_data$WGS_ID)
study_day_data$WGS_ID = gsub('--','-',study_day_data$WGS_ID)
for(i in 1:nrow(pt_data)){
  if(pt_data$Sample[i] %in% study_day_data$WGS_ID){
    pt_data$swab_day[i] = study_day_data$swab_day[study_day_data$WGS_ID == pt_data$Sample[i]]
  }else{
    pt_data$swab_day[i] = 'NA'
  }
}

for(i in 1:nrow(pt_data)){
  if(pt_data$swab_day[i] %in% daytomonth$Study.day.range){
    pt_data$month[i] = daytomonth$Month[daytomonth$Study.day.range == pt_data$swab_day[i]]
  }else{
    pt_data$month[i] = 'NA'
  }
}

month_df = data.frame(matrix(NA, nrow = length(tree$tip.label), ncol = 2), stringsAsFactors = F)
colnames(month_df) = c('Sample','Month')
month_df$Sample = tree$tip.label
for(i in 1:nrow(month_df)){
  month_df$Month[i] = pt_data$month[pt_data$Sample == month_df$Sample[i]]
} 
rownames(month_df) = month_df$Sample
month_df = month_df[,'Month', drop = F]

month_colors = colors$color[colors$topic == 'Month']
names(month_colors) = colors$name[colors$topic == 'Month']

epi_status_tree1 = epi_status_tree + new_scale_fill()

epistatus_month_tree = gheatmap(epi_status_tree1, month_df, colnames_angle = 85, colnames_offset_y = -7,
                                offset = 0.008, width = 0.029, font.size = 5)+
  scale_fill_manual(values = month_colors, breaks = names(month_colors))+ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Month of swab collection')

epistatus_month_tree2 = epistatus_month_tree + new_scale_fill()

## Figure 2A
#pdf(paste0('./figures/',Sys.Date(),'_crescent_epistatus_month_allisolates.pdf'), width = 15, height = 15)
epistatus_month_tree
#dev.off()


## Preparing genus-specific plasmid dataframes
## Klebsiella
kleb_col_df = data.frame(matrix(NA, nrow = length(kleb_tree$tip.label), ncol = 2), stringsAsFactors = F)
colnames(kleb_col_df) = c('Sample','Col')
kleb_col_df$Sample = kleb_tree$tip.label
for(i in 1:nrow(kleb_col_df)){
  kleb_col_df$Col[i] = col_plasmid_df$Col[rownames(col_plasmid_df) == kleb_col_df$Sample[i]]
}
rownames(kleb_col_df) = kleb_col_df$Sample
kleb_col_df = kleb_col_df[,'Col', drop = F]

kleb_inc_df = data.frame(matrix(NA, nrow = length(kleb_tree$tip.label), ncol = 2), stringsAsFactors = F)
colnames(kleb_inc_df) = c('Sample','Inc')
kleb_inc_df$Sample = kleb_tree$tip.label
for(i in 1:nrow(kleb_inc_df)){
  kleb_inc_df$Inc[i] = inc_plasmid_df$Inc[rownames(inc_plasmid_df) == kleb_inc_df$Sample[i]]
}
rownames(kleb_inc_df) = kleb_inc_df$Sample
kleb_inc_df = kleb_inc_df[,'Inc', drop = F]

kleb_otherplasmids_df = data.frame(matrix(NA, nrow = length(kleb_tree$tip.label), ncol = 2))
colnames(kleb_otherplasmids_df) = c('Sample','Other')
kleb_otherplasmids_df$Sample = kleb_tree$tip.label
for(i in 1:nrow(kleb_otherplasmids_df)){
  kleb_otherplasmids_df$Other[i] = other_plasmid_df$other_plasmid[rownames(other_plasmid_df) == kleb_otherplasmids_df$Sample[i]]
}
rownames(kleb_otherplasmids_df) = kleb_otherplasmids_df$Sample
kleb_otherplasmids_df = kleb_otherplasmids_df[,'Other', drop = F]

kleb_epi_carb_tree_fan2 = kleb_epi_carb_tree_fan + new_scale_fill()


#### AST data addition to Enterobacter, Klebsiella, and Escherichia trees - Figure 4A-C ####
ast_data = ast_data[ast_data$Sample %in% tree$tip.label,]
carb_ast_data = ast_data[,c('Sample','Imipenem_Interpretation','Meropenem_Interpretation')]

ordered_kleb_carb_asts = data.frame(matrix(NA, nrow = length(kleb_tree$tip.label), ncol = 1), stringsAsFactors = F)
colnames(ordered_kleb_carb_asts) = 'Sample'
ordered_kleb_carb_asts$Sample = kleb_tree$tip.label
for(i in 1:nrow(ordered_kleb_carb_asts)){
  if(ordered_kleb_carb_asts$Sample[i] %in% ast_data$Sample){
    ordered_kleb_carb_asts$Imipenem[i] = ast_data$Imipenem_Interpretation[ast_data$Sample == ordered_kleb_carb_asts$Sample[i]]
    ordered_kleb_carb_asts$Meropenem[i] = ast_data$Meropenem_Interpretation[ast_data$Sample == ordered_kleb_carb_asts$Sample[i]]
  }
}
rownames(ordered_kleb_carb_asts) = ordered_kleb_carb_asts$Sample
ordered_kleb_carb_asts = ordered_kleb_carb_asts[,c(2,3), drop = F]

ordered_esch_carb_asts = data.frame(matrix(NA, nrow = length(esch_tree$tip.label), ncol = 1), stringsAsFactors = F)
colnames(ordered_esch_carb_asts) = 'Sample'
ordered_esch_carb_asts$Sample = esch_tree$tip.label
for(i in 1:nrow(ordered_esch_carb_asts)){
  ordered_esch_carb_asts$Imipenem[i] = ast_data$Imipenem_Interpretation[ast_data$Sample == ordered_esch_carb_asts$Sample[i]]
  ordered_esch_carb_asts$Meropenem[i] = ast_data$Meropenem_Interpretation[ast_data$Sample == ordered_esch_carb_asts$Sample[i]]
}
rownames(ordered_esch_carb_asts) = ordered_esch_carb_asts$Sample
ordered_esch_carb_asts = ordered_esch_carb_asts[,c(2,3), drop = F]

ordered_ent_carb_asts = data.frame(matrix(NA, nrow = length(ent_tree$tip.label), ncol = 1), stringsAsFactors = F)
colnames(ordered_ent_carb_asts) = 'Sample'
ordered_ent_carb_asts$Sample = ent_tree$tip.label
for(i in 1:nrow(ordered_ent_carb_asts)){
  ordered_ent_carb_asts$Imipenem[i] = ast_data$Imipenem_Interpretation[ast_data$Sample == ordered_ent_carb_asts$Sample[i]]
  ordered_ent_carb_asts$Meropenem[i] = ast_data$Meropenem_Interpretation[ast_data$Sample == ordered_ent_carb_asts$Sample[i]]
}
rownames(ordered_ent_carb_asts) = ordered_ent_carb_asts$Sample
ordered_ent_carb_asts = ordered_ent_carb_asts[,c(2,3), drop = F]

ast_colors = colors$color[colors$topic == 'AST']
names(ast_colors) = colors$name[colors$topic == 'AST']

kleb_epi_carb_AST_tree_fan = gheatmap(kleb_epi_carb_tree_fan2, ordered_kleb_carb_asts, colnames_angle = 85, colnames_offset_y = -1,
                                      offset = 0.013, width = 0.1, font.size = 5)+
  scale_fill_manual(values = ast_colors, breaks = names(ast_colors))+
  ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Carbapenem susceptibilities')

#pdf(paste0('./figures/',Sys.Date(),'_crescent_epistatus_carb_AST_tree_kleb.pdf'), width = 15, height = 15)
kleb_epi_carb_AST_tree_fan
#dev.off()

kleb_epi_carb_AST_tree_fan2 = kleb_epi_carb_AST_tree_fan + new_scale_fill()

kleb_epi_carb_AST_col_tree_fan = gheatmap(kleb_epi_carb_AST_tree_fan2, kleb_col_df, colnames_angle = 85, colnames_offset_y = -1,
                                          offset = 0.019, width = 0.06, font.size = 5)+
  scale_fill_manual(values = simple_plasmid_colors, breaks = names(simple_plasmid_colors))+
  ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Col plasmid')

kleb_epi_carb_AST_col_tree_fan2 = kleb_epi_carb_AST_col_tree_fan + new_scale_fill()

kleb_epi_carb_AST_col_inc_tree_fan = gheatmap(kleb_epi_carb_AST_col_tree_fan2, kleb_inc_df, colnames_angle = 85, colnames_offset_y = -1,
                                              offset = 0.022, width = 0.06, font.size = 5)+
  scale_fill_manual(values = simple_plasmid_colors, breaks = names(simple_plasmid_colors))+
  ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Inc plasmid')

kleb_epi_carb_AST_col_inc_tree_fan2 = kleb_epi_carb_AST_col_inc_tree_fan + new_scale_fill()

kleb_epi_carb_AST_plasmid_tree_fan = gheatmap(kleb_epi_carb_AST_col_inc_tree_fan2, kleb_otherplasmids_df, colnames_angle = 85, colnames_offset_y = -1,
                                              offset = 0.0255, width = 0.06, font.size = 5)+
  scale_fill_manual(values = simple_plasmid_colors, breaks = names(simple_plasmid_colors))+
  ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Other plasmids')

## Figure 4B
#pdf(paste0('./figures/',Sys.Date(),'_crescent_epistatus_carb_AST_plasmid_tree_kleb.pdf'), width = 15, height = 15)
kleb_epi_carb_AST_plasmid_tree_fan
#dev.off()

## E. coli
esch_epi_carb_tree_fan2 = esch_epi_carb_tree_fan + new_scale_fill()

esch_epi_carb_AST_tree_fan = gheatmap(esch_epi_carb_tree_fan2, ordered_esch_carb_asts, colnames_angle = 85, colnames_offset_y = -1,
                                      offset = 0.0045, width = 0.08, font.size = 5)+
  scale_fill_manual(values = ast_colors, breaks = names(ast_colors))+
  ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Carbapenem susceptibilities')

#pdf(paste0('./figures/',Sys.Date(),'_crescent_epistatus_carb_AST_tree_esch.pdf'), width = 15, height = 15)
esch_epi_carb_AST_tree_fan
#dev.off()

esch_epi_carb_AST_tree_fan2 = esch_epi_carb_AST_tree_fan + new_scale_fill()

esch_epi_carb_AST_col_tree_fan = gheatmap(esch_epi_carb_AST_tree_fan2, esch_col_df, colnames_angle = 85, colnames_offset_y = -1,
                                          offset = 0.007, width = 0.06, font.size = 5)+
  scale_fill_manual(values = simple_plasmid_colors, breaks = names(simple_plasmid_colors))+
  ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Col plasmid')

esch_epi_carb_AST_col_tree_fan2 = esch_epi_carb_AST_col_tree_fan + new_scale_fill()

esch_epi_carb_AST_col_inc_tree_fan = gheatmap(esch_epi_carb_AST_col_tree_fan2, esch_inc_df, colnames_angle = 85, colnames_offset_y = -1,
                                              offset = 0.0087, width = 0.06, font.size = 5)+
  scale_fill_manual(values = simple_plasmid_colors, breaks = names(simple_plasmid_colors))+
  ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Inc plasmid')

esch_epi_carb_AST_col_inc_tree_fan2 = esch_epi_carb_AST_col_inc_tree_fan + new_scale_fill()

esch_epi_carb_AST_plasmid_tree_fan = gheatmap(esch_epi_carb_AST_col_inc_tree_fan2, esch_otherplasmids_df, colnames_angle = 85, colnames_offset_y = -1,
                                              offset = 0.0105, width = 0.06, font.size = 5)+
  scale_fill_manual(values = simple_plasmid_colors, breaks = names(simple_plasmid_colors))+
  ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Other plasmids')

## Figure 4A
#pdf(paste0('./figures/',Sys.Date(),'_crescent_epistatus_carb_AST_plasmid_tree_esch.pdf'), width = 15, height = 15)
esch_epi_carb_AST_plasmid_tree_fan
#dev.off()


## Enterobacter
ent_epi_carb_tree_fan2 = ent_epi_carb_tree_fan + new_scale_fill()
ent_epi_carb_AST_tree_fan = gheatmap(ent_epi_carb_tree_fan2, ordered_ent_carb_asts, colnames_angle = 85, colnames_offset_y = -1,
                                     offset = 0.018, width = 0.1, font.size = 5)+
  scale_fill_manual(values = ast_colors, breaks = names(ast_colors))+
  ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Carbapenem susceptibilities')

#pdf(paste0('./figures/',Sys.Date(),'_crescent_epistatus_carb_AST_tree_ent.pdf'), width = 15, height = 15)
ent_epi_carb_AST_tree_fan
#dev.off()

ent_epi_carb_AST_tree_fan2 = ent_epi_carb_AST_tree_fan + new_scale_fill()

ent_epi_carb_AST_col_tree_fan = gheatmap(ent_epi_carb_AST_tree_fan2, ent_col_df, colnames_angle = 85, colnames_offset_y = -1,
                                         offset = 0.033, width = 0.06, font.size = 5)+
  scale_fill_manual(values = simple_plasmid_colors, breaks = names(simple_plasmid_colors))+
  ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Col plasmid')

ent_epi_carb_AST_col_tree_fan2 = ent_epi_carb_AST_col_tree_fan + new_scale_fill()

ent_epi_carb_AST_col_inc_tree_fan = gheatmap(ent_epi_carb_AST_col_tree_fan2, ent_inc_df, colnames_angle = 85, colnames_offset_y = -1,
                                             offset = 0.042, width = 0.06, font.size = 5)+
  scale_fill_manual(values = simple_plasmid_colors, breaks = names(simple_plasmid_colors))+
  ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Inc plasmid')

ent_epi_carb_AST_col_inc_tree_fan2 = ent_epi_carb_AST_col_inc_tree_fan + new_scale_fill()

ent_epi_carb_AST_plasmid_tree_fan = gheatmap(ent_epi_carb_AST_col_inc_tree_fan2, ent_otherplasmids_df, colnames_angle = 85, colnames_offset_y = -1,
                                             offset = 0.05, width = 0.06, font.size = 5)+
  scale_fill_manual(values = simple_plasmid_colors, breaks = names(simple_plasmid_colors))+
  ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Other plasmids')

## Figure 4C
#pdf(paste0('./figures/',Sys.Date(),'_crescent_epistatus_carb_AST_plasmid_tree_ent.pdf'), width = 15, height = 15)
ent_epi_carb_AST_plasmid_tree_fan
#dev.off()


#### plotting all ARGs on full tree - Supplemental Figures 3A and 4B ####
amrfinder_subset = amrfinder_data
amrfinder_subset$classsample = paste0(amrfinder_subset$Class,'-',amrfinder_subset$Sample)
amrfinder_subset = amrfinder_subset[!duplicated(amrfinder_subset$classsample),]

ordered_amrfinder_class = data.frame(matrix('not present', nrow = length(tree$tip.label), ncol = 1+length(table(amrfinder_subset$Class))), stringsAsFactors = F)
colnames(ordered_amrfinder_class) = c('Sample',names(table(amrfinder_subset$Class)))
ordered_amrfinder_class$Sample = tree$tip.label
for(i in 1:nrow(ordered_amrfinder_class)){
  for(j in 2:ncol(ordered_amrfinder_class)){
    if(paste0(colnames(ordered_amrfinder_class)[j],'-',ordered_amrfinder_class$Sample[i]) %in% amrfinder_subset$classsample){
      ordered_amrfinder_class[i,j] = 'present'
    }
  }
}
rownames(ordered_amrfinder_class) = ordered_amrfinder_class$Sample
ordered_amrfinder_class = ordered_amrfinder_class[,c(2:ncol(ordered_amrfinder_class)), drop = F]

## Figure 3A
#pdf(paste0('./figures/',Sys.Date(),'_crescent_alltree_epistatus_month_argclass.pdf'), width = 20, height = 20)
gheatmap(epistatus_month_tree2, ordered_amrfinder_class, colnames_angle = 85, colnames_offset_y = -25,
         offset = 0.017, width = 0.6, font.size = 5)+
  scale_fill_manual(values = c('violetred4','grey90'), breaks = c('present','not present'))+
  ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Antibiotic resistance')
#dev.off()

amrfinder_subsetgene = amrfinder_data
amrfinder_subsetgene$genesample = paste0(amrfinder_subsetgene$Gene.symbol,'-',amrfinder_subsetgene$Sample)
amrfinder_subsetgene = amrfinder_subsetgene[!duplicated(amrfinder_subsetgene$genesample),]

ordered_amrfinder_gene = data.frame(matrix('not present', nrow = length(tree$tip.label), ncol = 1+length(table(amrfinder_subsetgene$Gene.symbol))), stringsAsFactors = F)
colnames(ordered_amrfinder_gene) = c('Sample',names(table(amrfinder_subsetgene$Gene.symbol)))
ordered_amrfinder_gene$Sample = tree$tip.label
for(i in 1:nrow(ordered_amrfinder_gene)){
  for(j in 2:ncol(ordered_amrfinder_gene)){
    if(paste0(colnames(ordered_amrfinder_gene)[j],'-',ordered_amrfinder_gene$Sample[i]) %in% amrfinder_subsetgene$genesample){
      ordered_amrfinder_gene[i,j] = 'present'
    }
  }
}
rownames(ordered_amrfinder_gene) = ordered_amrfinder_gene$Sample
ordered_amrfinder_gene = ordered_amrfinder_gene[,c(2:ncol(ordered_amrfinder_gene)), drop = F]


# plotting

thick_epistatus_tree = gheatmap(tip_treeplot, epi_df, colnames_angle = 85, colnames_offset_y = -5,
                                offset = 0, width = 0.7, font.size = 5)+
  scale_fill_manual(values = haca_colors, breaks = names(haca_colors))+
  ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Epidemiological status')

thick_epistatus_tree2 = thick_epistatus_tree + new_scale_fill()

## Supplemental Figure 4B 
#pdf(paste0('./figures/',Sys.Date(),'_crescent_alltree_epistatus_month_args.pdf'), width = 30, height = 20)
gheatmap(thick_epistatus_tree2, ordered_amrfinder_gene, colnames_angle = 85, colnames_offset_y = -5,
         offset = 0.3, width = 100, font.size = 5)+
  scale_fill_manual(values = c('darkslateblue','grey90'), breaks = c('present','not present'))+
  ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Antibiotic resistance gene/mutation')
#dev.off()


#### geom tile for environmental isolates by month - Figure 1B ####
env_data = data.frame(read_excel('./manuscript/Github_scripts/data/SupplementalTable1-SwabResults.xlsx', sheet = 2), stringsAsFactors = F)
env_data$Date = paste0(env_data$Month.of.collection, ' 2023')
env_data$Date = gsub('December 2023','December 2022', env_data$Date)
vitek_colors = colors$color[colors$topic == 'Vitek']
names(vitek_colors) = colors$name[colors$topic == 'Vitek']

## Figure 1B
#pdf(paste0('./figures/',Sys.Date(),'_CRESCENT_environmental_samples_tile.pdf'), width = 15, height = 10)
ggplot(env_data, aes(x = factor(Date, levels = c('December 2022','January 2023','March 2023','April 2023','August 2023',
                                                 'September 2023','October 2023','November 2023')), 
                     y = factor(Sink.source, levels = c('Sink 2 drain','Sink 3 drain','Sink 4 drain','Sink 6 drain','Sink 6 surface')), 
                     fill = Vitek.identification, color = Vitek.identification))+
  geom_point(alpha = 0.8, size = 15)+
  scale_fill_manual(values = vitek_colors, breaks = vitek_colors, guide = 'legend')+
  scale_color_manual(values = vitek_colors, breaks = vitek_colors, guide = 'legend')+
  theme_bw()+
  labs(fill = 'Vitek identification', y = 'Swab source', x = 'Swab month', color = 'Vitek identification')+
  theme(axis.text.y = element_text(size = 18, face = 'bold'), 
        axis.text.x = element_text(size = 18, angle = 45, vjust = 0.9, hjust = 0.9, face = 'bold'), 
        axis.title = element_text(size = 18, face = 'bold'),
        legend.position = 'bottom')
#dev.off()


#### number of plasmids per isolate by genus - Figure 4D and 4E ####
plasmids_per_isolate = data.frame(table(plasmid_data$Sample), stringsAsFactors = F)
colnames(plasmids_per_isolate) = c('Sample','All')
plasmids_per_isolate$genus = 'Klebsiella'

colnames(plasmid_data)[1] = 'Sample'
for(i in 1:nrow(plasmid_data)){
  plasmid_data$genus[i] = pt_data$genus[pt_data$Sample == plasmid_data$Sample[i]]
}

plasmids_per_isolate$genus[plasmids_per_isolate$Sample %in% plasmid_data$Sample[plasmid_data$genus == 'Escherichia']] = 'Escherichia'
plasmids_per_isolate$genus[plasmids_per_isolate$Sample %in% plasmid_data$Sample[plasmid_data$genus == 'Enterobacter']] = 'Enterobacter'
plasmids_per_isolate$status = 'Community'
plasmids_per_isolate$status[plasmids_per_isolate$Sample %in% pt_data$Sample[pt_data$Category == 'Environmental swab']] = 'Environmental swab'
plasmids_per_isolate$status[plasmids_per_isolate$Sample %in% pt_data$Sample[pt_data$Category == 'Healthcare associated']] = 'Healthcare-associated'
plasmids_per_isolate$status[plasmids_per_isolate$Sample %in% pt_data$Sample[pt_data$Category == 'Hospital acquired']] = 'Hospital-acquired'


genus_colors = c(colors$color[colors$name == 'Klebsiella pneumoniae'][1], colors$color[colors$name == 'Escherichia coli'][1],
                 colors$color[colors$name == 'Enterobacter cloacae'][1])
names(genus_colors) = c('Klebsiella','Escherichia','Enterobacter')


col_isolates = plasmid_data[plasmid_data$Col != 'none',]
col_isolates = col_isolates[!duplicated(col_isolates$Sample),]
col_isolates = col_isolates[,c('Sample','Plasmid','Col','genus')]
col_isolates$type = 'Col'

## Figure 4D
#pdf(paste0('./figures/',Sys.Date(),'_CRESCENT_Colplasmids.pdf'), width = 15, height = 10)
ggplot(col_isolates, aes(x = Col, fill = type, color = type))+
  geom_histogram(stat = 'count')+
  facet_grid(~genus, scales = 'free_y')+
  scale_fill_manual(values = simple_plasmid_colors, breaks = names(simple_plasmid_colors), guide = 'legend')+
  scale_color_manual(values = simple_plasmid_colors, breaks = names(simple_plasmid_colors), guide = 'legend')+
  theme_bw()+
  labs(fill = element_blank(), y = 'Isolate count', x = 'Col type', color = element_blank())+
  theme(axis.text.y = element_text(size = 18, face = 'bold'), 
        axis.text.x = element_text(size = 18, angle = 45, vjust = 0.9, hjust = 0.9, face = 'bold'), 
        axis.title = element_text(size = 18, face = 'bold'))
#dev.off()

# count of all isolates with an inc plasmid
inc_isolates = plasmid_data[plasmid_data$Inc != 'none',]
inc_isolates$IncSamp = paste0(inc_isolates$Inc, '-', inc_isolates$Sample)
inc_isolates = inc_isolates[!duplicated(inc_isolates$IncSamp),]
inc_isolates = inc_isolates[,c('Sample','Plasmid_Type','Inc','genus')]

## Figure 4E
#pdf(paste0('./figures/',Sys.Date(),'_CRESCENT_Incplasmids.pdf'), width = 15, height = 10)
ggplot(inc_isolates, aes(x = Inc, fill = Plasmid_Type, color = Plasmid_Type))+
  geom_histogram(stat = 'count')+
  facet_grid(~genus, scales = 'free_y')+
  scale_fill_manual(values = simple_plasmid_colors, breaks = names(simple_plasmid_colors), guide = 'legend')+
  scale_color_manual(values = simple_plasmid_colors, breaks = names(simple_plasmid_colors), guide = 'legend')+
  theme_bw()+
  labs(fill = element_blank(), y = 'Isolate count', x = 'Col type', color = element_blank())+
  theme(axis.text.y = element_text(size = 18, face = 'bold'), 
        axis.text.x = element_text(size = 18, angle = 45, vjust = 0.9, hjust = 0.9, face = 'bold'), 
        axis.title = element_text(size = 18, face = 'bold'))
#dev.off()


#### SNP visualization work - Figures 2B, 2C, and supplemental Figure 2 ####
## Figure 2A
# SNP clusters
esch_nodes = esch_data
# imputation to avoid errors when getting inverse
esch_snps$core_snps = esch_snps$core_snps + 1
esch_snps$sample1 = gsub('K','K-', esch_snps$sample1)
esch_snps$sample1 = gsub('--','',esch_snps$sample1)
esch_snps$sample2 = gsub('K','K-', esch_snps$sample2)
esch_snps$sample2 = gsub('--','',esch_snps$sample2)

esch_snps = esch_snps[esch_snps$sample1 %in% pt_data$Sample & esch_snps$sample2 %in% pt_data$Sample,]

colnames(esch_nodes)[1] = 'id'
esch_graph = graph_from_data_frame(d = esch_snps, vertices = esch_nodes, directed = F)

V(esch_graph)$color = haca_colors[V(esch_graph)$Category]

E(esch_graph)$width = (1 / E(esch_graph)$core_snps)*10

#pdf(paste0('./figures/',Sys.Date(),'_crescent_esch_clusters.pdf'), width = 14, height = 14)
plot(esch_graph, vertex.size = 4, vertex.label = NA)
#dev.off()

## Figure 2B
# SNP clusters
kleb_nodes = kleb_pt_data
kleb_snps$core_snps = kleb_snps$core_snps + 1
kleb_snps$sample1 = gsub('K','K-', kleb_snps$sample1)
kleb_snps$sample1 = gsub('--','',kleb_snps$sample1)
kleb_snps$sample2 = gsub('K','K-', kleb_snps$sample2)
kleb_snps$sample2 = gsub('--','',kleb_snps$sample2)

kleb_snps = kleb_snps[kleb_snps$sample1 %in% pt_data$Sample & kleb_snps$sample2 %in% pt_data$Sample,]

colnames(kleb_nodes)[1] = 'id'
kleb_graph = graph_from_data_frame(d = kleb_snps, vertices = kleb_nodes, directed = F)

V(kleb_graph)$color = haca_colors[V(kleb_graph)$Category]

E(kleb_graph)$width = (1 / E(kleb_graph)$core_snps)*10

#pdf(paste0('./figures/',Sys.Date(),'_crescent_kleb_clusters.pdf'), width = 14, height = 14)
plot(kleb_graph, vertex.size = 8, vertex.label = NA)
#dev.off()

## Supplemental Figure 2A
# adding column for genus
esch_snps$genus = 'Escerichia'

#pdf(paste0('./figures/',Sys.Date(),'_escherichia_coresnps_hist.pdf'))
ggplot(esch_snps, aes(x = core_snps, fill = genus))+
  geom_histogram()+
  scale_fill_manual(values = colors$color[colors$name == 'Escherichia coli'])+
  theme_bw()+
  labs(fill = 'Genus', x = 'Core genome SNPs', y = 'Count')+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        legend.title = element_text(size = 14, face = 'bold'), legend.text = element_text(size = 12))
#dev.off()

#pdf(paste0('./figures/',Sys.Date(),'_escherichia_coresnps_hist_inset.pdf'))
ggplot(esch_snps[esch_snps$core_snps < 100,], aes(x = core_snps, fill = genus))+
  geom_histogram()+
  scale_fill_manual(values = colors$color[colors$name == 'Escherichia coli'])+
  theme_bw()+
  labs(fill = 'Genus', x = 'Core genome SNPs', y = 'Count')+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        legend.title = element_text(size = 14, face = 'bold'), legend.text = element_text(size = 12))
#dev.off()

## Supplemental Figure 2B
# adding column for genus
kleb_snps$genus = 'Klebsiella'
#pdf(paste0('./figures/',Sys.Date(),'_klebsiella_coresnps_hist.pdf'))
ggplot(kleb_snps, aes(x = core_snps, fill = genus))+
  geom_histogram()+
  scale_fill_manual(values = colors$color[colors$name == 'Klebsiella pneumoniae'])+
  theme_bw()+
  labs(fill = 'Genus', x = 'Core genome SNPs', y = 'Count')+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        legend.title = element_text(size = 14, face = 'bold'), legend.text = element_text(size = 12))
#dev.off()

#pdf(paste0('./figures/',Sys.Date(),'_klebsiella_coresnps_hist_inset.pdf'))
ggplot(kleb_snps[kleb_snps$core_snps < 100,], aes(x = core_snps, fill = genus))+
  geom_histogram()+
  scale_fill_manual(values = colors$color[colors$name == 'Klebsiella pneumoniae'])+
  theme_bw()+
  labs(fill = 'Genus', x = 'Core genome SNPs', y = 'Count')+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        legend.title = element_text(size = 14, face = 'bold'), legend.text = element_text(size = 12))
#dev.off()

