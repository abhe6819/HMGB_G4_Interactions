#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=pqsProm
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abhe6819@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=3:00:00
#SBATCH --output=pqsProm.out
#SBATCH --error=pqsProm.err
pwd; hostname; date
echo "pqsProm"
module load R/3.6
R CMD BATCH pqsProm.R
date
