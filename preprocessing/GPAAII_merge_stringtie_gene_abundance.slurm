#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --partition=bii
#SBATCH --account=bii_gpaa
#SBATCH --cpus-per-task=8

module purge
module load anaconda

python merge_stringtie_gene_abundance.py -d /project/gpaa/output/combined_FPFCFS_protein-coding-biotype-only_known_20200115/stringtie /project/gpaa/output/AA_batch8_protein-coding-biotype-only_known_stranded_20201012/work/stringtie /project/gpaa/output/AA-batch6_protein-coding-biotype-only_known_stranded_20201014/work/stringtie/ /project/gpaa/output/AA_batch12_protein-coding-biotype-only_known_stranded_20201008/work/stringtie/ /project/gpaa/output/GPAAII_batch7-16_stringtie/ -o /project/gpaa/machine_learning/jainam_capstone/preprocessing/GPAA_all_FPKM.csv

