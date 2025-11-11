#!/bin/bash
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

module load R
R --vanilla < TCGA_ECM_expr_rpkm_Thac_20251111.R 
