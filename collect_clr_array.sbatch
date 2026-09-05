#!/bin/bash
#SBATCH --job-name=clr-collect
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=TwinsUK_target_exempt_clr_scores/collect-%A_%a.log

set -euo pipefail
cd /projects/genomic-ml/da2343/necromass
/projects/genomic-ml/da2343/soak-r/bin/Rscript collect_clr_array.R collect
