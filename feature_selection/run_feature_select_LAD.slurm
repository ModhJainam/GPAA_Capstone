#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --partition=bii
#SBATCH --account=bii_gpaa
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000

module purge
module load anaconda

python feature_select_LAD.py
