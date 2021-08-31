# transcriptome_sex

#### Repository for the sex-biased brain transcriptome project from Cayo Santiago macaques (bulk-tissue RNA-seq)

This repository contains scripts used in the analysis of sex effects in bulk-tissue RNA-seq data for the Cayo Santiago rhesus macaque population. 

Note that we ran most steps on the University of Washington ([Mox](https://wiki.cac.washington.edu/display/hyakusers/Hyak+mox+Overview)) high-performance computing cluster. We have aimed to generalize the code here by removing system-specific references to installed software and modules. Instead, we document required software and version numbers below (excluding standard Unix programs and R). For HPC systems, the required scripts and binaries must be in the PATH. The easiest way to do this is to use an existing module or to install your own. In these cases, the modules should be loaded prior to running the appropriate code below.

As Mox uses the [slurm](https://slurm.schedmd.com/documentation.html) scheduler, most code below should run on slurm systems with little or no modification. For non-slurm HPC systems, slurm scripts and environmental variables will need to be adjusted, though hopefully without too much hassle.

We ran most analysis steps using [R](https://cran.r-project.org/) (v4.1). We recommend the following utility or visualization packages to extend base R's functionality.

# Inputs

The following files are expected:

* Demultiplexed pair-end fastq files should be compressed with gzip and placed in the ```fastq/``` folder with the naming convention ```<library ID>.R1.fastq.gz``` (read 1) and ```<library ID>.R2.fastq.gz``` (read 2).

* An animal metadata file should be placed in ```data/cayo_brain_bulk_metadata_animals.tsv```

* A library metadata file should be placed in ```data/cayo_brain_bulk_metadata_technical.tsv```

* An animal social metrics file should be placed in ```data/social_metrics.csv```
  
# Pipeline
  
### Map reads with splice-aware aligner

* **Required software**: STAR (v2.5), SAMtools (v1.9), GATK (v4.1.2.0)

```
# Download and index genome
sbatch scripts/star_index.sh

# Map each library using STAR
sbatch --array=1-$(tail -n+2 data/cayo_brain_bulk_metadata_technical.tsv | wc -l | xargs) scripts/star_map.sh
```

### Merge alignments per genotype

* **Required software**: SAMtools (v1.9)

```
# Merge bam files for each genotype
sbatch --array=1-$(tail -n+2 data/cayo_brain_bulk_metadata_animals.tsv | wc -l | xargs) scripts/samtools_merge.sh
```

### Call and filter genotypes

* **Required software**: SAMtools (v1.9), GATK (v4.1.2.0), VCFtools (v0.1.16)

```
# Split variants into chromosomes (per genotype)
sbatch --array=1-$(tail -n+2 data/cayo_brain_bulk_metadata_animals.tsv | wc -l | xargs) scripts/samtools_split.sh

# Clean variants
sbatch --array=1-$(($(tail -n+2 data/cayo_brain_bulk_metadata_animals.tsv | wc -l | xargs)*20)) scripts/gatk_clean_reads.sh

# Call variants
sbatch --array=1-20 scripts/gatk_call_variants.sh

# Filter variants
sbatch --array=1-20 scripts/gatk_filter_variants.sh

# Concatenate variants across chromosomes
scripts/vcftools_concat.sh
```

### Compute relatedness with lcMLkin

* **Required software**: VCFtools (v0.1.16), lcMLkin (v20190218)

```
# Thin variants
scripts/vcftools_thin_variants.sh

# Call kinship
scripts/lcmlkin_call_kinship.sh
```

### Calculate sequencing stats

* **Required software**: SAMtools (v1.9), ea-utils (v1.04.807), GNU parallel (v20171122)

```
# Calculate sequencing stats across libraries
sbatch scripts/sequencing_stats_parallel.sh

# Summarize sequencing stats
scripts/sequencing_stats_summarize.sh
```

### Import and clean sample and library metadata

```
# Read, format, and clean animal and library metadata
# Create male and female sample lists
scripts/clean_metadata.R
```

### Quantify transcripts with kallisto

* **Required software**: kallisto (v0.43.1)

```
# Create sex specific transcriptomes and index
scripts/kallisto_transcriptomes.sh

# Count transcripts
sbatch --array=1-$(wc -l checkpoints/male_ids.txt | cut -d ' ' -f 1) scripts/kallisto_count_males.sh
sbatch --array=1-$(wc -l checkpoints/female_ids.txt | cut -d ' ' -f 1) scripts/kallisto_count_females.sh
```

### Import and combine expression data

* **Key libraries:** tximport, rdf5, biomaRt

```
# Import kallisto results into R and combine
scripts/kallisto_import.R
```

### Filter expression matrix

* **Key libraries:** biomaRt, limma

```
# Apply filters to gene expression dataset
scripts/filter_expression.R
```

### Visualize expression data (pre-modeling)

* **Key libraries:** variancePartition, umap, dendextend, ape, phangorn

```
# Visualize expression data
scripts/visualize_expression.R
```

### Fit linear mixed model(s)

* **Key libraries:** EMMREML

```
# Fit linear mixed effect model
scripts/emma_model.R
```

### Apply adaptive shrinkage

* **Key libraries:** mashr

```
# Refine effects and significance values with adaptive shrinkage
scripts/mashr_model.R
```

### Cell type enrichment

* **Key libraries:** EMMREML

```
# Fit linear mixed effect model
scripts/cell_type_enrichment.R
```

### Adjust expression data for cell type proportions 

* **Key libraries:** BRETIGEA

```
# Fit linear mixed effect model
scripts/adjust_expression.R
```

### Fit linear mixed model(s) to cell type corrected data

* **Key libraries:** EMMREML

```
# Fit linear mixed effect model
scripts/emma_model_cell_type.R
```

### Apply adaptive shrinkage to cell type corrected data

* **Key libraries:** mashr

```
# Refine effects and significance values with adaptive shrinkage
scripts/mashr_model_cell_type.R
```

### Disease, motif, and functional enrichment analyses

