#!/usr/bin/env bash
#SBATCH -A MRC-BSU2-SL2-CPU
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --job-name=finemap_INTERVAL
#SBATCH --array=1-58
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err


set -euo pipefail

CSV=/home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/perturb_exploratory/data/discordant_SNP_genes.csv

# Extract the line corresponding to this array task
# +1 because CSV has a header
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$CSV")

IFS=',' read -r SNP cis_gene <<< "$LINE"

echo "Task ${SLURM_ARRAY_TASK_ID}: SNP=${SNP}, gene=${cis_gene}"

Rscript fine_map_cis_eQTL_INTERVAL.R \
    --gene "${cis_gene}"
