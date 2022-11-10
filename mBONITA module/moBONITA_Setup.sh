#!/bin/sh
#SBATCH --partition=debug
#SBATCH -J bonitaStep1
#SBATCH -o bonitaStep1.log
#SBATCH -t 1:00:00

module load anaconda3/2020.07
source activate scBonita

python3 moBonita_kegg_parser.py --sep , --org hsa --data concatenated_datasets.csv