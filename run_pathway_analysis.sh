#!/bin/sh
#SBATCH --partition=debug
#SBATCH -J PA
#SBATCH -o pathway_analysis_log.txt
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH -c 1

module load anaconda3/2020.11
activate BONITA
#Template:
#python3 pathway_analysis_score_pathways.py Your_omics_data Your_condition_matrix Your_desired_contrasts -sep Separator_used_in_gmt_and_omics_data

#Example
#python3 pathway_analysis_score_pathways.py concatenated_datasets.csv concatenated_conditions.csv contrasts.csv -sep ,

python3 pathway_analysis_score_pathways_mBonita.py