#!/bin/sh
#SBATCH --partition=debug
#SBATCH --mem 16G
#SBATCH -o Figure4.log
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH -c 4

module load anaconda3/2020.11
activate BONITA

python3 Figure4.py