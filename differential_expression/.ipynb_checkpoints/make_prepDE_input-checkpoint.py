#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os

# desired format: ERR188021 <PATH_TO_ERR188021.gtf>

with open('AA_stringtie_dirs.txt', 'r') as f:
    AA_dirs = [line.rstrip('\n') for line in f]
    
    
AA_lines = []

for directory in AA_dirs:
    for d in os.listdir(directory):
        if os.path.isdir(directory+'/'+d):
            AA_lines.append(d+' '+directory+'/'+d+'/transcripts.gtf\n')


AA_prep_DE_in = open('AA_prep_DE_input.txt','w')

AA_prep_DE_in.writelines(AA_lines)

AA_prep_DE_in.close()




with open('LAD_stringtie_dirs.txt', 'r') as f:
    LAD_dirs = [line.rstrip('\n') for line in f]
    
    
LAD_lines = []

for directory in LAD_dirs:
    for d in os.listdir(directory):
        if os.path.isdir(directory+'/'+d):
            LAD_lines.append(d+' '+directory+'/'+d+'/transcripts.gtf\n')


LAD_prep_DE_in = open('LAD_prep_DE_input.txt','w')

LAD_prep_DE_in.writelines(LAD_lines)

LAD_prep_DE_in.close()