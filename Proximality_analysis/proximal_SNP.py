#!/usr/bin/env python
# coding: utf-8
#This scripts collects output files from Removed_positions directory output by proximal.py and zc_proximal. Output tsv file can be used with bcftools -filter to remove positions from all VCF files or can be included in pipeline-specific filters (e.g. vSNP Mbovis define filter file).
import pandas as pd
import numpy as np
import glob
import re

def proximal(z):
    vcf_file = pd.read_table('{s}_removed.vcf'.format(s=z),header=None,sep=",")
    vcf_file = vcf_file[:-1]
    vcf_file = vcf_file[1:]
    pos1 = vcf_file[0]
    pos2 = vcf_file[1]
    conc_pos = pd.concat([pos1,pos2])
    conc_pos.to_csv('{s}_proxpos.tsv'.format(s=z),header=False,index=False,sep="\t")

vcf_files = glob.glob('*.vcf')
vcf_files.sort()
prefixes = [re.sub(r'_removed.vcf','',vcf) for vcf in vcf_files]
prefixes.sort()
file_count=0
length = len(vcf_files)

while file_count < length:
    if file_count == length :
        pass
        file_count +=1
    else:
        sample = prefixes[file_count]
        print(sample)
        print('File count:{s}'.format(s=file_count))
        proximal_positions = proximal(sample)
        file_count = file_count + 1


prox_pos = glob.glob('*proxpos*')
df_prox_pos = [pd.read_table(vcf,header=None) for vcf in prox_pos]
df = pd.concat(df_prox_pos, axis = 0)
df = df.drop_duplicates(keep="first")
df = df.sort_values(by=[0])
df.to_csv('proximal_positions.tsv', header=None,index=None,sep='\t')

