#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from sklearn.preprocessing import OneHotEncoder

genomic_data_AA = "GPAA_AA_samples_preprocessed.csv"
metadata_AA = "GPAA_AA_batches0-16_metadata.csv"

genomic_df = pd.read_csv(genomic_data_AA).rename(columns={"Unnamed: 0": "sample_id"})
metadata_df = pd.read_csv(metadata_AA, usecols = [3, 4, 5, 8, 23])

merged = metadata_df.merge(genomic_df, how='inner', left_on = "UVA_sample_id", right_on = "sample_id").drop(['UVA_sample_id'], axis=1)
merged["binary_pathology"] = np.where(merged['pathology (assigned category)'].isin(['nl', 'fs']), 0, 1)
merged = merged.drop(['pathology (assigned category)'], axis = 1)
final = pd.get_dummies(merged, columns = ['SEX', 'RACE'])
final.rename(columns = {'SEX_F': 'female', 
                         'SEX_M': 'male',
                        'RACE_A': 'asian',
                        'RACE_B': 'black',
                        'RACE_H': 'hispanic',
                        'RACE_W': 'white'}, inplace = True)
final.to_csv("AA_labeled_batches0-16.csv", index=False)

'''
genomic_data_AA = "GPAA_AA_genomic_b13_preprocessed_asr_DE_genes.csv"
genomic_df = pd.read_csv(genomic_data_AA).rename(columns={"Unnamed: 0": "sample_id"})

merged = metadata_df.merge(genomic_df, how='inner', on = "sample_id")
merged.to_csv("AA_labeled_b13_fixed_asr_DE_genes.csv", index=False)

genomic_data_AA = "GPAA_AA_genomic_b13_preprocessed_noconfound_DE_genes.csv"
genomic_df = pd.read_csv(genomic_data_AA).rename(columns={"Unnamed: 0": "sample_id"})

merged = metadata_df.merge(genomic_df, how='inner', on = "sample_id")
merged.to_csv("AA_labeled_b13_fixed_noconfound_DE_genes.csv", index=False)

'''
genomic_data_LAD = "GPAA_LAD_samples_preprocessed.csv"
metadata_LAD = "GPAA_LAD_batches0-16_metadata.csv"


genomic_df = pd.read_csv(genomic_data_LAD).rename(columns={"Unnamed: 0": "sample_id"})
metadata_df = pd.read_csv(metadata_LAD, usecols = [3, 4, 5, 8, 23])

merged = metadata_df.merge(genomic_df, how='inner', left_on = "UVA_sample_id", right_on = "sample_id").drop(['UVA_sample_id'], axis=1)
merged["binary_pathology"] = np.where(merged['pathology (assigned category)'].isin(['nl', 'fs']), 0, 1)
merged = merged.drop(['pathology (assigned category)'], axis = 1)
final = pd.get_dummies(merged, columns = ['SEX', 'RACE'])
final.rename(columns = {'SEX_F': 'female', 
                         'SEX_M': 'male',
                        'RACE_A': 'asian',
                        'RACE_B': 'black',
                        'RACE_H': 'hispanic',
                        'RACE_W': 'white'}, inplace = True)
final.to_csv("LAD_labeled_batches0-16.csv", index=False)

'''
genomic_data_LAD = "GPAA_LAD_genomic_b13_preprocessed_asr_DE_genes.csv"
genomic_df = pd.read_csv(genomic_data_LAD).rename(columns={"Unnamed: 0": "sample_id"})

merged = metadata_df.merge(genomic_df, how='inner', on = "sample_id")
merged.to_csv("LAD_labeled_b13_fixed_asr_DE_genes.csv", index=False)

genomic_data_LAD = "GPAA_LAD_genomic_b13_preprocessed_noconfound_DE_genes.csv"
genomic_df = pd.read_csv(genomic_data_LAD).rename(columns={"Unnamed: 0": "sample_id"})

merged = metadata_df.merge(genomic_df, how='inner', on = "sample_id")
merged.to_csv("LAD_labeled_b13_fixed_noconfound_DE_genes.csv", index=False)
'''
