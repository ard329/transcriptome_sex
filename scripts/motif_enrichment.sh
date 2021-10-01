#!/bin/bash

module load contrib/homer/4.10

findMotifs.pl sb_genes_13regions.txt mmul10 sb_genes_13regions -start -1000 -end 300 -bg all_exp_genes_no_Y.txt -mset vertebrates -nomotif -nogo
