#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

genomic_data_AA = "GPAA_AA_genomic_b13_preprocessed.csv"
metadata_AA = "GPAA_AA_b13_metadata_fixed.csv"

genomic_df = pd.read_csv(genomic_data_AA).rename(columns={"Unnamed: 0": "sample_id"})
metadata_df = pd.read_csv(metadata_AA)

merged = metadata_df.merge(genomic_df, how='inner', left_on = "UVA_sample_id", right_on = "sample_id")
merged.to_csv("AA_labeled_b13_fixed.csv", index=False)

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

genomic_data_LAD = "GPAA_LAD_genomic_b13_preprocessed.csv"
metadata_LAD = "GPAA_LAD_b13_metadata_fixed.csv"


genomic_df = pd.read_csv(genomic_data_LAD).rename(columns={"Unnamed: 0": "sample_id"})
metadata_df = pd.read_csv(metadata_LAD)


merged = metadata_df.merge(genomic_df, how='inner', on = "sample_id")
merged.to_csv("LAD_labeled_b13_fixed.csv", index=False)

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