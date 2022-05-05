#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --partition=bii
#SBATCH --account=bii_gpaa
#SBATCH --cpus-per-task=8


module load anaconda/2020.11-py3.8
/project/biocomplexity/gpaa/software/stringtie-2.0.6.Linux_x86_64/prepDE.py -l 125 -g /project/biocomplexity/gpaa/machine_learning/jainam_capstone/differential_expression/gene_count_matrix_LAD.csv -t /project/biocomplexity/gpaa/machine_learning/jainam_capstone/differential_expression/transcript_count_matrix_LAD.csv --legend /project/biocomplexity/gpaa/machine_learning/jainam_capstone/differential_expression/legend_LAD.csv -i LAD_prep_DE_input.txt