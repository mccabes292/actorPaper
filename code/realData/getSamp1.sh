#!/bin/bash
#SBATCH --job-name=getSampMot
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=10000


module load sratoolkit/2.9.2

while read p; do fasterq-dump $p -p; done < SRRList.txt
