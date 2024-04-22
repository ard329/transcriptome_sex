#!/usr/bin/env Rscript

library(mashr)
library(biomaRt)

##############
# macaques
#############

mash.results = readRDS('mashr_results_d0.final.rds')
mash.ct.results = readRDS('mashr_results_cell_type_d0.final.rds')

# extract coefficients and p vals

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

# load disease expression data

diseases = read.csv('ASD-expression.csv')
ortho = readRDS('human_macaque_one2one.rds')
diseases2 = merge(diseases, ortho, by = 'ensembl_gene_id', all.x = T)

## remove Y chrom

mmul = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='mmulatta_gene_ensembl') 
gene2chrom = getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name'), filters = 'ensembl_gene_id', values = ensembl.gene.names, mart = mmul)
Ychrom = subset(gene2chrom, chromosome_name == "Y")
ensembl.gene.names = ensembl.gene.names[which(ensembl.gene.names %!in% Ychrom$ensembl_gene_id)]

# select data set (primary or cell type corrected)

mash.lfsr = mash.lfsr[which(rownames(mash.lfsr) %!in% Ychrom$ensembl_gene_id),]
mash.beta = mash.beta[which(rownames(mash.beta) %!in% Ychrom$ensembl_gene_id),]

mash.lfsr = mash.lfsr.ct[which(rownames(mash.lfsr.ct) %!in% Ychrom$ensembl_gene_id),]
mash.beta = mash.beta.ct[which(rownames(mash.beta.ct) %!in% Ychrom$ensembl_gene_id),]

# gene sex-biased genes

mashr.genes = rownames(mash.beta)
names(mashr.genes) = rownames(mash.beta)

all.region.fet = numeric(length=length(mashr.genes))
names(all.region.fet) = mashr.genes

fsr.cutoff.now = 0.05
fraction.cutoff.now = 1/15

all.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
  (sum(mash.lfsr[x,] < fsr.cutoff.now) >= fraction.cutoff.now * length(colnames(mash.beta))) && sum(mash.beta[x,][mash.lfsr[x,] < fsr.cutoff.now] > 0) >= fraction.cutoff.now * length(colnames(mash.beta))
}))))] = 1
all.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
  (sum(mash.lfsr[x,] < fsr.cutoff.now) >= fraction.cutoff.now * length(colnames(mash.beta))) && sum(mash.beta[x,][mash.lfsr[x,] < fsr.cutoff.now] < 0) >= fraction.cutoff.now * length(colnames(mash.beta))
}))))] = -1
table(all.region.fet)

all.region.join = data.frame(mmulatta_homolog_ensembl_gene = mashr.genes, direction = as.integer(all.region.fet))
all.region.do.pass = merge(all.region.join, diseases2, by='mmulatta_homolog_ensembl_gene')

# run

for(i in 1:4){
  
  x = subset(all.region.do.pass, study == 'Voineagu')[,c(2,i+4)]
  rownames(x) = subset(all.region.do.pass, study == 'Voineagu')$ensembl_gene_id
  
  table(x)
  
  f = table(x)[1,]
  n = table(x)[2,]
  m = table(x)[3,]
  
  contingency.matrix.m = matrix(c(m[2], f[2] + n[2], m[1], f[1] + n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))
  contingency.matrix.f = matrix(c(f[2], m[2] + n[2], f[1], m[1] + n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))
  contingency.matrix.m
  contingency.matrix.f
  
  inc.fet.test = fisher.test(contingency.matrix.m)
  dec.fet.test = fisher.test(contingency.matrix.f)
  
  print(colnames(x)[2])
  print(inc.fet.test)
  print(dec.fet.test)
  
}

for(i in 1:3){
  
  x = subset(all.region.do.pass, study == 'Gupta')[,c(2,i+8)]
  rownames(x) = subset(all.region.do.pass, study == 'Gupta')$ensembl_gene_id
  
  f = table(x)[1,]
  n = table(x)[2,]
  m = table(x)[3,]
  
  contingency.matrix.m = matrix(c(m[2], f[2] + n[2], m[1], f[1] + n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))
  contingency.matrix.f = matrix(c(f[2], m[2] + n[2], f[1], m[1] + n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))
  contingency.matrix.m
  contingency.matrix.f
  
  inc.fet.test = fisher.test(contingency.matrix.m)
  dec.fet.test = fisher.test(contingency.matrix.f)
  
  print(colnames(x)[2])
  print(inc.fet.test)
  print(dec.fet.test)
  
}

for(i in 1:2){
  
  x = subset(all.region.do.pass, study == 'Gandal')[,c(2,i+11)]
  rownames(x) = subset(all.region.do.pass, study == 'Gandal')$ensembl_gene_id
  
  f = table(x)[1,]
  n = table(x)[2,]
  m = table(x)[3,]
  
  contingency.matrix.m = matrix(c(m[2], f[2] + n[2], m[1], f[1] + n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))
  contingency.matrix.f = matrix(c(f[2], m[2] + n[2], f[1], m[1] + n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))
  contingency.matrix.m
  contingency.matrix.f
  
  inc.fet.test = fisher.test(contingency.matrix.m)
  dec.fet.test = fisher.test(contingency.matrix.f)
  
  print(colnames(x)[2])
  print(inc.fet.test)
  print(dec.fet.test)
  
}

##############
# humans
#############

mash.hsap = readRDS('gtex_mashr_results_sex.rds')

# extract coefficients and p vals

hsap.beta = get_pm(mash.hsap)
hsap.lfsr = get_lfsr(mash.hsap)
hsap.berr = get_psd(mash.hsap)
hsap.sbet = hsap.beta / hsap.berr

library(stringr)

ensembl.gene.names = str_sub(rownames(hsap.beta),1,15)

# load disease expression data

diseases = read.csv('ASD-expression.csv')
ortho = readRDS('human_macaque_one2one.rds')
diseases2 = merge(diseases, ortho, by = 'ensembl_gene_id', all.x = T)

## remove Y chrom

hsap = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl') 
gene2chrom = getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name'), filters = 'ensembl_gene_id', values = ensembl.gene.names, mart = hsap)
Ychrom = subset(gene2chrom, chromosome_name == "Y")
ensembl.gene.names = ensembl.gene.names[which(ensembl.gene.names %!in% Ychrom$ensembl_gene_id)]

rownames(hsap.lfsr) = str_sub(rownames(hsap.lfsr),1,15)
rownames(hsap.beta) = str_sub(rownames(hsap.beta),1,15)
mash.lfsr = hsap.lfsr[which(rownames(hsap.lfsr) %!in% Ychrom$ensembl_gene_id),]
mash.beta = hsap.beta[which(rownames(hsap.beta) %!in% Ychrom$ensembl_gene_id),]

# gene sig genes

mashr.genes = rownames(mash.beta)
names(mashr.genes) = rownames(mash.beta)

all.region.fet = numeric(length=length(mashr.genes))
names(all.region.fet) = mashr.genes

fsr.cutoff.now = 0.05
fraction.cutoff.now = 1/10 

all.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
  (sum(mash.lfsr[x,] < fsr.cutoff.now) >= fraction.cutoff.now * length(colnames(mash.beta))) && sum(mash.beta[x,][mash.lfsr[x,] < fsr.cutoff.now] > 0) >= fraction.cutoff.now * length(colnames(mash.beta))
}))))] = 1
all.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
  (sum(mash.lfsr[x,] < fsr.cutoff.now) >= fraction.cutoff.now * length(colnames(mash.beta))) && sum(mash.beta[x,][mash.lfsr[x,] < fsr.cutoff.now] < 0) >= fraction.cutoff.now * length(colnames(mash.beta))
}))))] = -1
table(all.region.fet)

all.region.join = data.frame(ensembl_gene_id = mashr.genes, direction = as.integer(all.region.fet))
all.region.do.pass = merge(all.region.join, diseases2, by='ensembl_gene_id')

for(i in 1:4){
  
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
  
  inc.fet.test = fisher.test(contingency.matrix.m)
  dec.fet.test = fisher.test(contingency.matrix.f)
  
  print(colnames(x)[2])
  print(inc.fet.test)
  print(dec.fet.test)
  
}

for(i in 1:3){
  
  x = subset(all.region.do.pass, study == 'Gupta')[,c(2,i+7)]
  rownames(x) = subset(all.region.do.pass, study == 'Gupta')$ensembl_gene_id
  
  f = table(x)[1,]
  n = table(x)[2,]
  m = table(x)[3,]
  
  contingency.matrix.m = matrix(c(m[2], f[2] + n[2], m[1], f[1] + n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))
  contingency.matrix.f = matrix(c(f[2], m[2] + n[2], f[1], m[1] + n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))
  contingency.matrix.m
  contingency.matrix.f
  
  inc.fet.test = fisher.test(contingency.matrix.m)
  dec.fet.test = fisher.test(contingency.matrix.f)
  
  print(colnames(x)[2])
  print(inc.fet.test)
  print(dec.fet.test)
  
}

for(i in 1:2){
  
  x = subset(all.region.do.pass, study == 'Gandal')[,c(2,i+10)]
  rownames(x) = subset(all.region.do.pass, study == 'Gandal')$ensembl_gene_id
  
  f = table(x)[1,]
  n = table(x)[2,]
  m = table(x)[3,]
  
  contingency.matrix.m = matrix(c(m[2], f[2] + n[2], m[1], f[1] + n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))
  contingency.matrix.f = matrix(c(f[2], m[2] + n[2], f[1], m[1] + n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))
  contingency.matrix.m
  contingency.matrix.f
  
  inc.fet.test = fisher.test(contingency.matrix.m)
  dec.fet.test = fisher.test(contingency.matrix.f)
  
  print(colnames(x)[2])
  print(inc.fet.test)
  print(dec.fet.test)
  
}
