#!/bin/sh
#SBATCH --partition=debug
#SBATCH -J moBONITA_Setup
#SBATCH -o moBONITA_Setup
#SBATCH -t 1:00:00

module load anaconda3/2020.07
source activate scBonita

python3 moBonita_kegg_parser.py -org hsa -sep , --data YOUR_DATAFILE