#!/bin/bash
#
#SBATCH --cpus-per-task=32
#SBATCH --job-name=date_treesnotB117
#SBATCH --output=date_treesnotB117.out
#SBATCH --account=lomannj-covid-19-realtime-epidemiology
#SBATCH --qos=lomannj
#SBATCH --time=24:00:00
#SBATCH --nodes=1
eval "$(conda shell.bash hook)"
conda activate r-environment
Rscript /cephfs/covid/bham/climb-covid19-geidelbergl/R/compare_lineages_d1_NOTB.1.1.7.R
