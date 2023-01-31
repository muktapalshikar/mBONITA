#!/bin/sh
#SBATCH --partition=standard
#SBATCH -a 1-5
#SBATCH -J GArun
#SBATCH -o GA_out_%A_%a.txt
#SBATCH -t 120:00:00
#SBATCH -n 1
#SBATCH -c 4

module load anaconda3/2020.07
source activate BONITA
python3 pathway_analysis_score_nodes.py 'metaNetwork.gpickle' $SLURM_ARRAY_TASK_ID
echo "ran GA"
