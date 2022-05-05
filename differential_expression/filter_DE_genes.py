#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd

f = open('/gpfs/gpfs0/project/gpaa/machine_learning/may_2021_repo/differential_expression/noconfound_AA_DE_genes.txt', 'r')
noconfound_AA = f.readlines()
f = open('/gpfs/gpfs0/project/gpaa/machine_learning/may_2021_repo/differential_expression/noconfound_LAD_DE_genes.txt', 'r')
noconfound_LAD = f.readlines()
f = open('/gpfs/gpfs0/project/gpaa/machine_learning/may_2021_repo/differential_expression/asr_AA_DE_genes.txt', 'r')
asr_AA = f.readlines()
f = open('/gpfs/gpfs0/project/gpaa/machine_learning/may_2021_repo/differential_expression/asr_LAD_DE_genes.txt', 'r')
asr_LAD = f.readlines()

LAD_preprocessed = pd.read_csv('LAD_labeled_b13_fixed.csv')
AA_preprocessed = pd.read_csv('AA_labeled_b13_fixed.csv')

noconfound_AA_split = ['sample_id','age','white','black','hispanic','male','female','binary_pathology']
noconfound_LAD_split = ['sample_id','age','white','black','hispanic','male','female','binary_pathology']
asr_AA_split = ['sample_id','age','white','black','hispanic','male','female','binary_pathology']
asr_LAD_split = ['sample_id','age','white','black','hispanic','male','female','binary_pathology']

for line in noconfound_AA:
    noconfound_AA_split.append(line.rstrip().split('|')[0])
    
for line in noconfound_LAD:
    noconfound_LAD_split.append(line.rstrip().split('|')[0])

for line in asr_AA:
    asr_AA_split.append(line.rstrip().split('|')[0])
    
for line in asr_LAD:
    asr_LAD_split.append(line.rstrip().split('|')[0])
    
AA_asr_DE_preprocessed = AA_preprocessed.filter(items=asr_AA_split)
LAD_asr_DE_preprocessed = LAD_preprocessed.filter(items=asr_LAD_split)
AA_noconfound_DE_prepreocessed = AA_preprocessed.filter(items=noconfound_AA_split)
LAD_noconfound_DE_preprocessed = LAD_preprocessed.filter(items=noconfound_LAD_split)

AA_asr_DE_preprocessed.to_csv('AA_asr_DE_labeled.csv', index=False)
LAD_asr_DE_preprocessed.to_csv('LAD_asr_DE_labeled.csv', index=False)

AA_noconfound_DE_prepreocessed.to_csv('AA_noconfound_DE_labeled.csv', index=False)
LAD_noconfound_DE_preprocessed.to_csv('LAD_noconfound_DE_labeled.csv', index=False)
