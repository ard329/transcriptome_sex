library(limma)
library(biomaRt)
library(expss)
library(ggplot2)
library(ggrepel)
library(egg)

################
## run LMMs
################

## load human counts and metadata
# add AGE column to 'gtex_meta.rds' using 'SAMPID' column if not already completed

hum_exp = readRDS('human_gtex_d0.rds')

hum_meta = readRDS('gtex_meta.rds')
hum_meta = hum_meta[complete.cases(hum_meta[,c('SMRIN','SMNTRNRT','TRISCHD')]),]
hum_exp = hum_exp[,which(colnames(hum_exp) %in% hum_meta$SAMPID)]

## load macaque counts and metadata

mac_exp = readRDS('macaque_d0.rds')
mac_meta = readRDS('cayo_bulkbrain_combined_metadata.rds')

## convert macaque gene ID to human gene ID

mmul = useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='mmulatta_gene_ensembl')
mmul.hsap = getBM(attributes=c('ensembl_gene_id','chromosome_name','external_gene_name','hsapiens_homolog_ensembl_gene','hsapiens_homolog_associated_gene_name','hsapiens_homolog_orthology_type','hsapiens_homolog_chromosome'),mart = mmul)
mmul.hsap.o2o = subset(mmul.hsap,hsapiens_homolog_orthology_type == 'ortholog_one2one')

mac_genes = data.frame(ensembl_gene_id = rownames(mac_exp$counts))
for(i in 1:length(mac_genes$ensembl_gene_id)){
  mac_genes$external_gene_name[i] = vlookup(mac_genes$ensembl_gene_id[i],
                                            dict = mmul.hsap.o2o, 
                                            lookup_column = 1,
                                            result_column = 3)}
rownames(mac_exp$counts) = mac_genes$external_gene_name

## combine human and macaque counts
## calculate normalization factors
## voom normalization with quality weights

# run mergeDGEList.R

d0 = mergeDGEList(hum_exp, mac_exp)
d0 = calcNormFactors(d0, method = 'TMM')

cpm = cpm(d0$counts)
cpm_mean = rowMeans(cpm, na.rm = T)
keep = names(cpm_mean[which(cpm_mean > 10)])

d0 = d0[keep,]

v = voom(d0, plot = T)
saveRDS(v, 'combo_voom.rds')

v = voomWithQualityWeights(d0, plot = T)
saveRDS(v, 'combo_voom_qw.rds')

## remove covariate effects (within each species)

v = readRDS('combo_voom_qw.rds')

hum_exp2 = v[,which(colnames(v) %in% hum_meta$SAMPID)]
mac_exp2 = v[,which(colnames(v) %in% mac_meta$LID)]

hum_exp_adj = removeBatchEffect(hum_exp2, covariates = cbind(hum_meta$SMRIN,hum_meta$SMNTRNRT,hum_meta$TRISCHD,hum_meta$AGE))
mac_exp_adj = removeBatchEffect(mac_exp2, batch = mac_meta$Library.batch, batch2 = mac_meta$ordinal.rank, covariates = cbind(mac_meta$exact_age_years,mac_meta$RIN))

# combine adjusted expression and meta data

combined_exp_adj = cbind(hum_exp_adj, mac_exp_adj) 
idx = match(colnames(v$E), colnames(combined_exp_adj))
combined_exp_adj = combined_exp_adj[,idx]
v$E = combined_exp_adj

colnames(hum_meta)[1] = 'LID'
colnames(hum_meta)[3] = 'Region'
colnames(hum_meta)[8] = 'sex'
hum_meta$sex = as.character(hum_meta$sex)
hum_meta$sex[which(hum_meta$sex == 'female')] = 'f'
hum_meta$sex[which(hum_meta$sex == 'male')] = 'm'
combo_meta = rbind(cbind(hum_meta[,c('LID','sex','Region')], species = 'human'), 
                   cbind(mac_meta[,c('LID','sex','Region')], species = 'macaque'))
#check
table(table(combo_meta$LID))

## run linear mixed models for each region

comparisons = c(
  AMY = 'Brain - Amygdala',
  ACCg = 'Brain - Cortex',
  CN = 'Brain - Striatum',
  dmPFC = 'Brain - Cortex',
  DG = 'Brain - Hippocampus',
  CA3 = 'Brain - Hippocampus',
  VMH = 'Brain - Hypothalamus',
  Pu = 'Brain - Striatum'
)

res_all = data.frame()

for(i in 1:length(comparisons)){
  
  print(comparisons[i])
  
  meta_now = subset(combo_meta, Region == names(comparisons[i]) | Region == comparisons[i])
  
  v_now = v[,which(colnames(v) %in% meta_now$LID)]
  idx = match(colnames(v_now$E), meta_now$LID)
  meta_now = meta_now[idx,]
  
  dupcor = duplicateCorrelation(v_now, model.matrix(~ 0 + sex, meta_now), block = meta_now$species)
  fitdupcor = lmFit(v_now, model.matrix(~ 0 + sex, meta_now), block = meta_now$species, correlation = dupcor$consensus.correlation)
  contr = makeContrasts(sexm - sexf, levels = colnames(coef(fitdupcor)))
  tmp = contrasts.fit(fitdupcor, contr)
  tmp = eBayes(tmp)
  res = data.frame(beta = tmp$coefficients, pval = tmp$p.value, 
                   se = sqrt(tmp$s2.post) * tmp$stdev.unscaled,
                   tissue = names(comparisons[i]))
  colnames(res) = c('beta','pval','se','tissue')
  res$padj = p.adjust(res$pval, method = 'BH')
  res$gene = rownames(res)
  rownames(res) = NULL
  res_all = rbind(res_all, res)
  
}

head(res_all)
saveRDS(res_all, 'LMM_results.rds')


