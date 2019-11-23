#!/bin/bash
#SBATCH --job-name=salmonSamp1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=10000

module load salmon
salmon quant -p 6 -i /proj/milovelab/mccabe/proj/GTEx/data/gencodeV19/gencode.v19_salmon-0.13.1 -l A --gcBias -r $1.fastq.gz --validateMappings -o $1_quant
