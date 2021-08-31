#!/usr/bin/env Rscript

#############
## gene level
#############

# download metadata and kallisto summarized to the gene level

txi.gene = readRDS('kallisto_genes.RDS')
meta = readRDS('cayo_bulkbrain_combined_metadata.RDS')

txi.gene$counts = txi.gene$counts[,which(colnames(txi.gene$counts) %in% meta$LID)]
txi.gene$abundance = txi.gene$abundance[,which(colnames(txi.gene$abundance) %in% meta$LID)]
txi.gene$length = txi.gene$length[,which(colnames(txi.gene$length) %in% meta$LID)]

# remove low count genes (TPM)
# genes with mean >= 10 TPM (transcripts per million) 
# for males OR females within each region

samples.by.region.and.sex = split(meta$LID, list(meta$Region,meta$sex), drop = TRUE)

tpm.cutoff = 10

keep.genes = lapply(samples.by.region.and.sex,function(x) {
  names(which(rowMeans(txi.gene$abundance[,x]) >= tpm.cutoff))
})

regions = levels(as.factor(meta$Region))
keep.new = list()
keep.genes2 = unlist(keep.genes)
for (i in 1:length(regions)){
  now = keep.genes2[grepl(regions[i], names(keep.genes2))]
  keep.new[[i]] = unique(now)
}
names(keep.new) = regions

keep = Reduce(union,keep.genes)

saveRDS(keep.new, 'keep_genes.rds')

# normalize

library(limma)
library(edgeR)

counts = txi.gene$counts
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
head(d0$samples)
d0 = d0[keep,]

v=voom(d0)

e = v$E[keep,rownames(meta)]

saveRDS(e, 'filtered_expression_matrix.rds')

###################
## transcript level
###################

# download metadata and kallisto transcript level data

txi.transcripts = readRDS('kallisto_transcripts.RDS')
meta = readRDS('cayo_bulkbrain_combined_metadata.RDS')

txi.transcripts$counts = txi.transcripts$counts[,which(colnames(txi.transcripts$counts) %in% meta$LID)]
txi.transcripts$abundance = txi.transcripts$abundance[,which(colnames(txi.transcripts$abundance) %in% meta$LID)]
txi.transcripts$length = txi.transcripts$length[,which(colnames(txi.transcripts$length) %in% meta$LID)]

# only include genes in overall gene set

library(biomaRt)

mmul = useEnsembl(biomart = 'ensembl',dataset='mmulatta_gene_ensembl') ##run multiple times
tx2gene.chrom = getBM(attributes=c('ensembl_transcript_id_version', 'ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position'), filters = 'ensembl_transcript_id_version', values = rownames(txi.transcripts$abundance), mart = mmul)
keep.genes = readRDS('keep_genes.rds')
keep = Reduce(union,keep.genes)
keep.trans.list = subset(tx2gene.chrom, ensembl_gene_id %in% keep) # N=34824

# list of complete transcripts

library(Biostrings)

seqs = getSequence(id = keep.trans.list$ensembl_transcript_id_version, type = 'ensembl_transcript_id_version', seqType = 'coding', mart = mmul)
seqs = subset(seqs, coding != 'Sequence unavailable') # N=34636
for (i in 1:length(seqs$coding)){
  pepnow = try(translate(DNAString(seqs$coding[i])), TRUE)
  if(isTRUE(class(pepnow)=="try-error")){seqs$pep[i] = 'error'} else {seqs$pep[i] = as.character(pepnow)}
}

library(stringr)

keep.trans.list.complete = seqs
keep.trans.list.complete = subset(keep.trans.list.complete, pep != 'error') # N=34635
keep.trans.list.complete = keep.trans.list.complete[which(str_sub(keep.trans.list.complete$coding,1,3) == 'ATG'),] # N=32966
keep.trans.list.complete = keep.trans.list.complete[which(str_sub(keep.trans.list.complete$pep,-1) == '*'),] # N=32691

# transcripts with mean >= 1 TPM (transcripts per million) 
# for males OR females within each region
# and also part of gene that is in original filtering

samples.by.region.and.sex = split(meta$LID, list(meta$Region,meta$sex), drop = TRUE)

tpm.cutoff = 1

keep.trans = lapply(samples.by.region.and.sex,function(x) {
  names(which(rowMeans(txi.transcripts$abundance[,x]) >= tpm.cutoff 
              & rownames(txi.transcripts$abundance) %in% keep.trans.list.complete$ensembl_transcript_id_version))
})

regions = levels(as.factor(meta$Region))
keep.new = list()
keep.trans2 = unlist(keep.trans)
for (i in 1:length(regions)){
  now = keep.trans2[grepl(regions[i], names(keep.trans2))]
  keep.new[[i]] = unique(now)
}
names(keep.new) = regions

keep.t = Reduce(union,keep.trans)

saveRDS(keep.new, 'keep_transcripts_complete.rds') # N=27511

# normalize 

library(limma)
library(edgeR)

counts = txi.transcripts$counts
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
head(d0$samples)
d0 = d0[keep.t,]

v=voom(d0)

e = v$E[keep.t,rownames(meta)]

saveRDS(e, 'filtered_transcript_expression_matrix_complete.rds')





