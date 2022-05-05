#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --partition=bii
#SBATCH --account=bii_gpaa
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000

module purge
module load anaconda
module load gcc/7.1.0 openmpi/3.1.4 R/4.1.1

#python /project/gpaa/machine_learning/jainam_capstone/preprocessing/separate_AA_LAD_samples.py

#R --vanilla /project/gpaa/machine_learning/jainam_capstone/preprocessing/preprocess_GPAA_AA_and_LAD.r

python /project/gpaa/machine_learning/jainam_capstone/preprocessing/merge_metadata_and_genomicdata_AA_LAD_split.py
