#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 15:41:40 2021

@author: cwp5j
"""

import os
import pandas as pd
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description='Parse through input directories and merge stringtie gene abundance outputs')
parser.add_argument('-t', "--gene_type", default='FPKM', choices=['FPKM', 'TPM'], help='parse for FPKM or TPM, default FPKM')
parser.add_argument('-o', "--output_file", required=True, help='name of output csv, will be written to current directory unless otherwise specified')
parser.add_argument('-d', "--directories", nargs='+', required=True, help='absolute pathways to stringtie output directories')

args = parser.parse_args()

#print(args.gene_type)
#print(args.output_file)
#print(args.directories)

df_dict = defaultdict(list)
columns_list = []

last_key = ""
total_num_samples = 0

for out_dir in args.directories:
    print("processing: " + out_dir)
    for d in os.listdir(out_dir):
        if os.path.isdir(out_dir+'/'+d):
            tsv = pd.read_csv(out_dir+"/"+d+"/gene_abund.tab", sep='\t')
        
            fpkm = tsv[args.gene_type].tolist()
            fpkm_flt = [float(x) for x in fpkm]
            gene_ID = tsv['Gene ID'].tolist()
            gene_ID_str = [str(x) for x in gene_ID]
            columns_list.append(d)
            
            total_num_samples += 1
        
            for i in range(0, len(gene_ID_str)):
                df_dict[gene_ID_str[i]].append(fpkm_flt[i])
                last_key = gene_ID_str[i]
                

keys_to_remove = []    


for key in df_dict:
    if len(df_dict[key]) != total_num_samples: # if for any reason uneven amounts are read across samples, cannot be considered

        keys_to_remove.append(key)
        
for key in keys_to_remove:
    df_dict.pop(key)
    
final = pd.DataFrame.from_dict(df_dict, orient='index', columns=columns_list)
final_transposed = final.T
print("Processed " + str(total_num_samples) + " samples")
print(final_transposed.head())


final_transposed.to_csv(args.output_file, index=True, index_label = 'sample_id')
print("Wrote file " + args.output_file)