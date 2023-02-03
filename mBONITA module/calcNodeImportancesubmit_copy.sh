#!/bin/sh
#SBATCH --partition=preempt
#SBATCH -a 1
#SBATCH -J GArun
#SBATCH -o GA_out_%A_%a.txt
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -c 1

module load anaconda3/2020.07
source activate BONITA
python3 pathway_analysis_score_nodes.py $1 $SLURM_ARRAY_TASK_ID
echo "ran GA"
