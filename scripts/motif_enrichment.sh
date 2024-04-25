#!/bin/bash

module load contrib/homer/4.10

loadGenome.pl -name mmul10 -org mmul -fasta mmul10.dna.assem.fa -gtf mmul10.anno.gtf

# macaque (primary)
findMotifs.pl list-sex-biased-genes.txt mmul10 output_directory -start -1000 -end 300 -bg list-all-expressed-genes.txt -mset vertebrates -nomotif -nogo
findMotifs.pl list-male-biased-genes.txt mmul10 output_directory -start -1000 -end 300 -bg list-all-expressed-genes.txt -mset vertebrates -nomotif -nogo
findMotifs.pl list-female-biased-genes.txt mmul10 output_directory -start -1000 -end 300 -bg list-all-expressed-genes.txt -mset vertebrates -nomotif -nogo

# macaque (cell-type corrected)
findMotifs.pl list-sex-biased-genes-cell-type-corrected.txt mmul10 output_directory -start -1000 -end 300 -bg list-all-expressed-genes.txt -mset vertebrates -nomotif -nogo
findMotifs.pl list-male-biased-genes-cell-type-corrected.txt mmul10 output_directory -start -1000 -end 300 -bg list-all-expressed-genes.txt -mset vertebrates -nomotif -nogo
findMotifs.pl list-female-biased-genes-cell-type-corrected.txt mmul10 output_directory -start -1000 -end 300 -bg list-all-expressed-genes.txt -mset vertebrates -nomotif -nogo
