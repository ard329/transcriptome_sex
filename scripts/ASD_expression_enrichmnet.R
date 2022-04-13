#!/usr/bin/env Rscript

##############
# load results
#############

e.keep = readRDS('filtered_expression_matrix.rds')

keep.genes = readRDS('keep_genes.rds')
meta = readRDS('cayo_bulkbrain_combined_metadata.rds')
meta = meta[order(match(meta$LID, colnames(e.keep))),]

mash.results = readRDS('mashr_results.rds')
mash.ct.results = readRDS('mashr_results_cell.rds')
mash.hsap = readRDS('gtex_mashr_results_sex.rds')

#################################
# extract coefficients and p vals
# select data set
#################################

library(mashr)

# macaque primary analyses
mash.beta = get_pm(mash.results)
mash.lfsr = get_lfsr(mash.results)
mash.sbet = mash.beta / get_psd(mash.results)
mash.beta = mash.beta[,region.levels]
mash.lfsr = mash.lfsr[,region.levels]
mash.sbet = mash.sbet[,region.levels]
ensembl.gene.names = unique(unlist(keep.genes))

# macaque cell type corrected
mash.beta.ct = get_pm(mash.ct.results)
mash.lfsr.ct = get_lfsr(mash.ct.results)
mash.sbet.ct = mash.beta.ct / get_psd(mash.ct.results)
mash.beta.ct = mash.beta.ct[,region.levels]
mash.lfsr.ct = mash.lfsr.ct[,region.levels]
mash.beta.ct = mash.beta.ct[,region.levels]
ensembl.gene.names = unique(unlist(keep.genes))

# human gtex
hsap.beta = get_pm(mash.hsap)
hsap.lfsr = get_lfsr(mash.hsap)
hsap.berr = get_psd(mash.hsap)
hsap.sbet = hsap.beta / hsap.berr
library(stringr)
ensembl.gene.names = str_sub(rownames(hsap.beta),1,15)

# select data set

mash.lfsr = mash.lfsr[which(rownames(mash.lfsr) %!in% Ychrom$ensembl_gene_id),]
mash.beta = mash.beta[which(rownames(mash.beta) %!in% Ychrom$ensembl_gene_id),]

mash.lfsr = mash.lfsr.ct[which(rownames(mash.lfsr.ct) %!in% Ychrom$ensembl_gene_id),]
mash.beta = mash.beta.ct[which(rownames(mash.beta.ct) %!in% Ychrom$ensembl_gene_id),]

rownames(hsap.lfsr) = str_sub(rownames(hsap.lfsr),1,15)
rownames(hsap.beta) = str_sub(rownames(hsap.beta),1,15)
mash.lfsr = hsap.lfsr[which(rownames(hsap.lfsr) %!in% Ychrom$ensembl_gene_id),]
mash.beta = hsap.beta[which(rownames(hsap.beta) %!in% Ychrom$ensembl_gene_id),]


##############################
# load disease expression data
##############################

diseases = read.csv('Haney2021.csv')

# get orthologs
ortho = readRDS('human_macaque_one2one.rds')
diseases2 = merge(diseases, ortho, by = 'ensembl_gene_id', all.x = T)

## remove Y chrom

library(biomaRt)

# macaques
mmul = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='mmulatta_gene_ensembl') 
gene2chrom = getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name'), filters = 'ensembl_gene_id', values = ensembl.gene.names, mart = mmul)
Ychrom = subset(gene2chrom, chromosome_name == "Y")
ensembl.gene.names = ensembl.gene.names[which(ensembl.gene.names %!in% Ychrom$ensembl_gene_id)]

# humans
hsap = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl') 
gene2chrom = getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name'), filters = 'ensembl_gene_id', values = ensembl.gene.names, mart = hsap)
Ychrom = subset(gene2chrom, chromosome_name == "Y")
ensembl.gene.names = ensembl.gene.names[which(ensembl.gene.names %!in% Ychrom$ensembl_gene_id)]

# limit to subset of regions

reg = c("dmPFC","dlPFC","vmPFC","vlPFC","ACCg","M1","STS","V1") # cortical only
reg = c("AMY","CA3","DG","CN","Pu","LGN","VMH") # subcortical only

mash.beta = mash.beta[,which(colnames(mash.beta) %in% reg)]
mash.lfsr = mash.lfsr[,which(colnames(mash.lfsr) %in% reg)]

# gene sig genes

mashr.genes = rownames(mash.beta)
names(mashr.genes) = rownames(mash.beta)

all.region.fet = numeric(length=length(mashr.genes))
names(all.region.fet) = mashr.genes

fsr.cutoff.now = 0.05
#fraction.cutoff.now = 1/10 
#fraction.cutoff.now = 1/15
fraction.cutoff.now = 1/length(reg)

all.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
  (sum(mash.lfsr[x,] < fsr.cutoff.now) >= fraction.cutoff.now * length(colnames(mash.beta))) && sum(mash.beta[x,][mash.lfsr[x,] < fsr.cutoff.now] > 0) >= fraction.cutoff.now * length(colnames(mash.beta))
}))))] = 1
all.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
  (sum(mash.lfsr[x,] < fsr.cutoff.now) >= fraction.cutoff.now * length(colnames(mash.beta))) && sum(mash.beta[x,][mash.lfsr[x,] < fsr.cutoff.now] < 0) >= fraction.cutoff.now * length(colnames(mash.beta))
}))))] = -1
table(all.region.fet)

all.region.join = data.frame(mmulatta_homolog_ensembl_gene = mashr.genes, direction = as.integer(all.region.fet))
all.region.do.pass = merge(all.region.join, diseases2, by='mmulatta_homolog_ensembl_gene')

all.region.join = data.frame(ensembl_gene_id = mashr.genes, direction = as.integer(all.region.fet))
all.region.do.pass = merge(all.region.join, diseases2, by='ensembl_gene_id')

# run on Werling

for(i in 1:4){
  
  #x = subset(all.region.do.pass, study == 'Voineagu')[,c(2,i+4)]
  x = subset(all.region.do.pass, study == 'Voineagu')[,c(2,i+3)]
  rownames(x) = subset(all.region.do.pass, study == 'Voineagu')$ensembl_gene_id
  
  table(x)

  f = table(x)[1,]
  n = table(x)[2,]
  m = table(x)[3,]

  contingency.matrix.m = matrix(c(m[2], f[2] + n[2], m[1], f[1] + n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))
  contingency.matrix.f = matrix(c(f[2], m[2] + n[2], f[1], m[1] + n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))
  contingency.matrix.m
  contingency.matrix.f
  
  inc.fet.test = fisher.test(contingency.matrix.m,alternative='greater')
  dec.fet.test = fisher.test(contingency.matrix.f,alternative='greater')
  
  print(colnames(x)[2])
  print(inc.fet.test)
  print(dec.fet.test)
  
}

for(i in 1:3){
  
  #x = subset(all.region.do.pass, study == 'Gupta')[,c(2,i+8)]
  x = subset(all.region.do.pass, study == 'Gupta')[,c(2,i+7)]
  rownames(x) = subset(all.region.do.pass, study == 'Gupta')$ensembl_gene_id
  
  f = table(x)[1,]
  n = table(x)[2,]
  m = table(x)[3,]
  
  contingency.matrix.m = matrix(c(m[2], f[2] + n[2], m[1], f[1] + n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))
  contingency.matrix.f = matrix(c(f[2], m[2] + n[2], f[1], m[1] + n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))
  contingency.matrix.m
  contingency.matrix.f
  
  inc.fet.test = fisher.test(contingency.matrix.m,alternative='greater')
  dec.fet.test = fisher.test(contingency.matrix.f,alternative='greater')
  
  print(colnames(x)[2])
  print(inc.fet.test)
  print(dec.fet.test)

}

# combine up and down

colnames(all.region.do.pass)

#x = all.region.do.pass[,c(1,5,8,10,11)] # up
#x = all.region.do.pass[,c(1,6,7,9)] #down
x = all.region.do.pass[,c(1,4,7,9,10)] # up
x = all.region.do.pass[,c(1,5,6,8)] #down

x$sum = rowSums(x[,c(-1)], na.rm = TRUE)
x2 = x %>% group_by(ensembl_gene_id) %>% summarise(association = sum(sum))
#x2 = x %>% group_by(mmulatta_homolog_ensembl_gene) %>% summarise(association = sum(sum))
x2$association = ifelse(x2$association > 0, 1, 0)
x2 = merge(x2, all.region.join, by = 'ensembl_gene_id')
#x2 = merge(x2, all.region.join, by = 'mmulatta_homolog_ensembl_gene')
rownames(x2) = x2$ensembl_gene_id
x2 = x2[,c(-1)]
x2 = x2[,c(2,1)]

table(x2)

f = table(x2)[1,]
n = table(x2)[2,]
m = table(x2)[3,]

contingency.matrix.m = matrix(c(m[2], f[2] + n[2], m[1], f[1] + n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))
contingency.matrix.m
contingency.matrix.f = matrix(c(f[2], m[2] + n[2], f[1], m[1] + n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))
contingency.matrix.f

fisher.test(contingency.matrix.m,alternative='greater')
fisher.test(contingency.matrix.f,alternative='greater')

##############
# run on Haney
##############

x = all.region.do.pass[,c(2,6)] #associated
x = all.region.do.pass[,c(2,7)] #up
x = all.region.do.pass[,c(2,8)] #down
x = all.region.do.pass[,c(2,5)] #associated
x = all.region.do.pass[,c(2,6)] #up
x = all.region.do.pass[,c(2,7)] #down

rownames(x) = all.region.do.pass$mmulatta_homolog_ensembl_gene
rownames(x) = all.region.do.pass$ensembl_gene_id

table(x)

f = table(x)[1,]
n = table(x)[2,]
m = table(x)[3,]

contingency.matrix.m = matrix(c(m[2], f[2] + n[2], m[1], f[1] + n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))
contingency.matrix.f = matrix(c(f[2], m[2] + n[2], f[1], m[1] + n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))
contingency.matrix.all = matrix(c(f[2] + m[2], n[2], f[1] + m[1], n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))

contingency.matrix.m
contingency.matrix.f
contingency.matrix.all

fisher.test(contingency.matrix.m,alternative='greater')
fisher.test(contingency.matrix.f,alternative='greater')
fisher.test(contingency.matrix.all,alternative='greater')

##################
## check cell type
##################

library(BRETIGEA)
library(biomaRt)
library(expss)

markers = markers_df_brain

hsap = useEnsembl(biomart = 'ensembl',dataset='hsapiens_gene_ensembl',mirror='www') 
gene_conv = getBM(attributes=c('ensembl_gene_id','external_gene_name'), mart = hsap)

for (i in 1:length(markers$markers)){
  markers$id[i] = vlookup(markers$markers[i], dict = gene_conv, result_column = 1,lookup_column = 2)}
markers = markers[which(!is.na(markers$id)),]
colnames(markers) = c("name","cell","markers")
table(markers$cell)
cells = levels(as.factor(markers$cell))
cells

# haney check

x = all.region.do.pass[,c(3,2,7)] # up
f = subset(x, direction == -1 & Up == 1)$ensembl_gene_id
m = subset(x, direction == 1 & Up == 1)$ensembl_gene_id

x = all.region.do.pass[,c(3,2,8)] # down
f = subset(x, direction == -1 & Down == 1)$ensembl_gene_id
m = subset(x, direction == 1 & Down == 1)$ensembl_gene_id

for(i in 1:length(cells)){
  print(cells[i])
  markersnow = subset(markers, cell == cells[i])$markers
  print(sum(f %in% markersnow))
  print(sum(m %in% markersnow))
}


