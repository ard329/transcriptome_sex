#!/usr/bin/env Rscript

#######################
# load data and results
#######################

e.keep = readRDS('filtered_expression_matrix.rds')
keep.genes = readRDS('keep_genes.rds')
meta = readRDS('cayo_bulkbrain_combined_metadata.rds')
meta = meta[order(match(meta$LID, colnames(e.keep))),]

emma.results = readRDS('emma_results.rds')
emma.ct.results = readRDS('emma_results_cell_type.rds')
emma.t.results = readRDS('emma_results_transcript.rds')
mash.results = readRDS('mashr_results.rds')
mash.ct.results = readRDS('mashr_results_cell_type.rds')
mash.t.results = readRDS('mashr_results_transcript.rds')

#################################
# extract coefficients and p vals
#################################

emma.beta = emma.results[,paste('beta',predictor,sep='.'),]
emma.pval = emma.results[,paste('pval',predictor,sep='.'),]
emma.qval = apply(emma.pval,2,function(x) p.adjust(x,'fdr'))
emma.beta = emma.beta[,region.levels]
emma.pval = emma.pval[,region.levels]
emma.qval = emma.qval[,region.levels]

emma.beta.ct = emma.ct.results[,paste('beta',predictor,sep='.'),]
emma.pval.ct = emma.ct.results[,paste('pval',predictor,sep='.'),]
emma.qval.ct = apply(emma.pval.ct,2,function(x) p.adjust(x,'fdr'))
emma.beta.ct = emma.beta.ct[,region.levels]
emma.pval.ct = emma.pval.ct[,region.levels]
emma.qval.ct = emma.qval.ct[,region.levels]

get_pm=function(x){x$result$PosteriorMean}
get_lfsr=function(x){x$result$lfsr}
get_psd=function(x){x$result$PosteriorSD}

mash.beta = get_pm(mash.results)
mash.lfsr = get_lfsr(mash.results)
mash.sbet = mash.beta / get_psd(mash.results)
mash.beta = mash.beta[,region.levels]
mash.lfsr = mash.lfsr[,region.levels]
mash.sbet = mash.sbet[,region.levels]

mash.beta.ct = get_pm(mash.ct.results)
mash.lfsr.ct = get_lfsr(mash.ct.results)
mash.sbet.ct = mash.beta.ct / get_psd(mash.ct.results)
mash.beta.ct = mash.beta.ct[,region.levels]
mash.lfsr.ct = mash.lfsr.ct[,region.levels]
mash.beta.ct = mash.beta.ct[,region.levels]

ensembl.gene.names = unique(unlist(keep.genes))

#################
## get gene lists
#################

mashr.genes = rownames(mash.beta)
names(mashr.genes) = rownames(mash.beta)
all.region = numeric(length=length(mashr.genes))
names(all.region) = mashr.genes
betanow = mash.beta
lsfrnow = mash.lfsr

all.region[names(which(unlist(lapply(mashr.genes,function(x) {
  (sum(lsfrnow[x,] < fsr.cutoff) >= 1/15 * length(keep.genes)) && sum(betanow[x,][lsfrnow[x,] < fsr.cutoff] > 0) >= 1/15 * length(keep.genes)
}))))] = 1
all.region[names(which(unlist(lapply(mashr.genes,function(x) {
  (sum(lsfrnow[x,] < fsr.cutoff) >= 1/15 * length(keep.genes)) && sum(betanow[x,][lsfrnow[x,] < fsr.cutoff] < 0) >= 1/15 * length(keep.genes)
}))))] = -1
table(all.region)

mb = names(all.region[which(all.region == 1)])
fb = names(all.region[which(all.region == -1)])
sb = c(mb, fb)
all = names(all.region)

######################################
## cell type marker data from BRETIGEA
######################################

library(BRETIGEA)
library(biomaRt)
library(expss)

markers = markers_df_brain

hsap = useEnsembl(biomart = 'ensembl',dataset='hsapiens_gene_ensembl',mirror='www') 
gene_conv = getBM(attributes=c('ensembl_gene_id','external_gene_name','mmulatta_homolog_ensembl_gene','mmulatta_homolog_orthology_type'), mart = hsap)
gene_conv = subset(gene_conv, mmulatta_homolog_orthology_type == 'ortholog_one2one')

for (i in 1:length(markers$markers)){
  markers$id[i] = vlookup(markers$markers[i], dict = gene_conv, result_column = 3,lookup_column = 2)}
markers = markers[which(!is.na(markers$id)),]
colnames(markers) = c("name","cell","markers")
markers = markers[markers$markers %in% row.names(e.keep),]
markers = markers[markers$markers %!in% markers$markers[duplicated(markers$markers)],]
table(markers$cell)
markers = markers[,c(3,2)]
colnames(markers) = c('macID', 'cellName')

ct_markers_all = markers
head(ct_markers_all)
table(ct_markers_all$cellName)
cells = levels(as.factor(ct_markers_all$cellName))
cells

###################
## fisher tests
###################

celltypes = data.frame()

for (i in 1:length(cells)){
  cellnow = cells[i]
  
  sexnow = fb
  markers = unique(subset(ct_markers_all, cellName == cellnow)$macID)
  combo.cell = length(which(markers %in% sexnow))
  cell.all = markers[which(markers %in% all)]
  cell.only = length(which(cell.all %!in% sexnow))
  sex.all = length(sexnow)
  sex.only = sex.all- combo.cell
  neither = length(all) - combo.cell - cell.only - sex.only
  mat = matrix(c(combo.cell, sex.only, cell.only, neither), nrow = 2, ncol = 2)
  t1 = fisher.test(mat, alternative = 'greater')
  
  sexnow = mb
  markers = unique(subset(ct_markers_all, cellName == cellnow)$macID)
  combo.cell = length(which(markers %in% sexnow))
  cell.all = markers[which(markers %in% all)]
  cell.only = length(which(cell.all %!in% sexnow))
  sex.all = length(sexnow)
  sex.only = sex.all- combo.cell
  neither = length(all) - combo.cell - cell.only - sex.only
  mat = matrix(c(combo.cell, sex.only, cell.only, neither), nrow = 2, ncol = 2)
  t2 = fisher.test(mat, alternative = 'greater')
  
  sexnow = sb
  markers = unique(subset(ct_markers_all, cellName == cellnow)$macID)
  combo.cell = length(which(markers %in% sexnow))
  cell.all = markers[which(markers %in% all)]
  cell.only = length(which(cell.all %!in% sexnow))
  sex.all = length(sexnow)
  sex.only = sex.all- combo.cell
  neither = length(all) - combo.cell - cell.only - sex.only
  mat = matrix(c(combo.cell, sex.only, cell.only, neither), nrow = 2, ncol = 2)
  t3 = fisher.test(mat, alternative = 'greater')
  
  celltypes[i,1] = t1$p.value
  celltypes[i,2] = t1$estimate
  celltypes[i,3] = t2$p.value
  celltypes[i,4] = t2$estimate
  celltypes[i,5] = t3$p.value
  celltypes[i,6] = t3$estimate
  
}

rownames(celltypes) = cells
celltypes

celltypes$fbadj = p.adjust(celltypes$V1, method = 'BH')
celltypes$mbadj = p.adjust(celltypes$V3, method = 'BH')
celltypes$sbadj = p.adjust(celltypes$V5, method = 'BH')
celltypes = celltypes[,c(2,1,7,4,3,8,6,5,9)]
colnames(celltypes) = c('female OR', 'female p val', 'female p adj', 'male OR', 'male p val', 'male p adj', 'overall OR', 'overall p val', 'overall p adj')
celltypes

#######
## plot
#######

library(ggplot2)

plotdata = data.frame(genes = c(fb, mb))
plotdata$bias = ifelse(plotdata$genes %in% fb, "Female", "Male")
table(plotdata$bias)
length(unique(plotdata$genes))

plotdata_try = plotdata
colnames(plotdata_try) = c('macID','bias')
plotdata_try = merge(plotdata_try, ct_markers_all, by = 'macID')
table(plotdata_try$bias)
length(unique(plotdata_try$macID))

for(i in 1:length(plotdata$genes)){plotdata$`Cell Type`[i] = vlookup(plotdata$genes[i], dict = ct_markers_all, lookup_column = 1, result_column = 2)}
plotdata2 <- aggregate(.~`Cell Type`+bias, plotdata, length)
plotdata2$Proportion = ifelse(plotdata2$bias == "Female", plotdata2$genes/sum(subset(plotdata2, bias == "Female")$genes), plotdata2$genes/sum(subset(plotdata2, bias == "Male")$genes))
plotdata2$`Cell Type` = as.factor(plotdata2$`Cell Type`)
levels(plotdata2$`Cell Type`) = c("Astrocytes","Endothelial Cells","Microglia","Neurons","Oligodendrocytes","OPCs")
ggplot(data=plotdata2, aes(x=bias, y=genes, fill=`Cell Type`)) + 
  geom_bar(stat="identity",position="stack") + theme_classic() + 
  xlab("Bias") + ylab("Count") + scale_fill_manual(values=region.colors)
ggplot(data=plotdata2, aes(x=bias, y=Proportion, fill=`Cell Type`)) + 
  geom_bar(stat="identity",position="stack") + theme_classic() + 
  xlab("Bias") + ylab("Proportion") + scale_fill_brewer(palette='PuOr') +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=18), legend.text = element_text(size=18), legend.title = element_text(size=18))

