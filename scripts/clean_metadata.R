#!/usr/bin/env Rscript

# Read in metadata
meta.technical = read.delim('cayo_brain_bulk_metadata_technical_final.txt',stringsAsFactors=FALSE)
meta.animal = read.delim('cayo_brain_bulk_metadata_animals_final.tsv',stringsAsFactors=FALSE)

# Split dataset by sex and scale ranks for each sex
meta.animal.split = split(meta.animal,meta.animal$sex)

# Scale so that highest ranked are 1 and lowest ranked are 0
#meta.animal = do.call(rbind,lapply(meta.animal.split,function(x) {
#  x = within(x,{
#    rank.scaled = 1 - ((rank - min(rank,na.rm=TRUE)) / (max(rank,na.rm=TRUE) - min(rank,na.rm=TRUE)))
#  })
#}))
#row.names(meta.animal) = NULL

# Factorize ordinal ranks
meta.animal$ordinal.rank = ordered(meta.animal$ordinal.rank,levels=c('L','M','H'))
meta.animal$integer.rank = as.integer(meta.animal$ordinal.rank)

fix.yn = function(YN,Y='Y',N='N') {
  # YN = yes/no vector
  as.logical(as.integer(factor(YN,levels=c(N,Y)))-1)
}

fix.time = function(TIME,DATE,CUTOFF,TZ='America/Puerto_Rico') {
  # TIME = time (03:00), DATE = date, CUTOFF = cutoff for AM/PM, TZ = time zone
  TIME[!is.na(TIME) & nchar(TIME) < 5] = NA
  TIME[!is.na(TIME) & as.integer(substr(TIME,1,2)) <6] = paste0(as.integer(substr(TIME[!is.na(TIME) & as.integer(substr(TIME,1,2)) <6],1,2))+12,substr(TIME[!is.na(TIME) & as.integer(substr(TIME,1,2)) <6],3,5))
  TIME[!is.na(TIME)] = paste0(TIME[!is.na(TIME)],':00')
  TIME[!is.na(TIME)] = paste(DATE[!is.na(TIME)],TIME[!is.na(TIME)])
  TIME = as.POSIXlt(TIME,tz=TZ)
  TIME
}

# Format certain columns
meta.animal = within(meta.animal,{
  sex = factor(sex)
  blood_in_brain = fix.yn(blood_in_brain)
  with_infant = fix.yn(with_infant)
  weight_kg = as.numeric(weight_kg)
  group.at.birth = factor(group.at.birth,levels=na.omit(unique(meta.animal$group.at.birth)[order(nchar(unique(meta.animal$group.at.birth)),unique(meta.animal$group.at.birth))]))
  euthanasia_time = fix.time(euthanasia_time,euthanasia_date)
  brain_extraction_time = fix.time(brain_extraction_time,euthanasia_date)
  brain_frozen_time = fix.time(brain_frozen_time,euthanasia_date)
  ketamine_time = fix.time(ketamine_time,euthanasia_date)
  birth_date = as.Date(birth_date)
  trapping_date = as.Date(trapping_date)
  euthanasia_date = as.Date(euthanasia_date)
  frozen_time = NA
  frozen_time[!is.na(brain_frozen_time) & !is.na(euthanasia_time)] = brain_frozen_time[!is.na(brain_frozen_time) & !is.na(euthanasia_time)] - euthanasia_time[!is.na(brain_frozen_time) & !is.na(euthanasia_time)]
})

cayo = merge(meta.technical,meta.animal,by.x='Individual',by.y='animal_id',all.x=TRUE)

rownames(cayo) = cayo$LID
cayo = cayo[order(rownames(cayo)),]

# Add in sequencing stats

library.stats = read.delim('library_sequencing_stats.txt',stringsAsFactors=FALSE)
rownames(library.stats) = library.stats$library
cayo = data.frame(cayo,library.stats[rownames(cayo),])

# Remove libraries
cayo = subset(cayo,!(qual_mean < 36 & perc_duplicate > 50) & perc_mapped > 50 & abs(perc_G-perc_C) < 4)

# Get lists of male and female samples
m = subset(cayo, sex == 'm')$LID
f = subset(cayo, sex == 'f')$LID
write.table(m, file = 'checkpoints/male_ids.txt', col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(f, file = 'checkpoints/female_ids.txt', col.names = FALSE, row.names = FALSE, quote = FALSE)

# Factorize regions
cayo$Region = factor(cayo$Region,levels=region.levels)

saveRDS(cayo,file='checkpoints/cayo_bulkbrain_combined_metadata.rds')
