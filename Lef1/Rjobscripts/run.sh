#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=pqsShuf
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=abhe6819@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=7:00:00
#SBATCH --output=pqsS.out
#SBATCH --error=pqsS.err
pwd; hostname; date
echo "pqsS"
module load R/3.6
R CMD BATCH pqsshuffled.R
date
