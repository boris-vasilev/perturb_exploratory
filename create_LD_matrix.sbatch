#!/usr/bin/env bash
#SBATCH -A MRC-BSU2-SL2-CPU
#SBATCH -p icelake-himem
#SBATCH -N 1
#SBATCH --job-name=ld_matrix
#SBATCH --array=1-58
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

. /etc/profile.d/modules.sh

module load openjdk/17.0.11_9/gcc/szfinnph
source ~/.bashrc
conda activate magenpy_env

set -euo pipefail

CSV=/home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/perturb_exploratory/data/discordant_SNP_genes.csv

# Extract the line corresponding to this array task
# +1 because CSV has a header
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$CSV")

IFS=',' read -r SNP cis_gene <<< "$LINE"

echo "Task ${SLURM_ARRAY_TASK_ID}: SNP=${SNP}, gene=${cis_gene}"

python plot_cis_eQTL.py \
    --gene "${cis_gene}" \
    --window 1000
