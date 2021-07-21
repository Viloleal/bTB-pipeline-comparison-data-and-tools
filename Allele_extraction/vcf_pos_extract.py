#!/usr/bin/env python
# coding: utf-8



import numpy as np
import pandas as pd
import glob
import sys
import vcf
import csv
import os
from os import listdir
from os import path
import subprocess
import re
import itertools
import matplotlib
from matplotlib_venn import venn3, venn3_circles, venn3_unweighted
from matplotlib import pyplot as plt




def read(f):
    reader = vcf.Reader(open(f))
    df = pd.DataFrame([vars(r) for r in reader])
    out = df.merge(pd.DataFrame(df.INFO.tolist()),
                   left_index=True, right_index=True)
    return out

def vcf_merge(z):
    vsnp = read('{s}.vcf'.format(s=z))
    vsnp = vsnp[['CHROM','POS','REF','ALT']]
    vsnp["CHROM"] = vsnp["CHROM"].replace('NC_002945.4','{s}_vsnp'.format(s=sample))
#REMOVE LIST FROM ALT
    vsnp['ALT'] = vsnp['ALT'].apply(lambda x: ', '.join(map(str, x)))
    vsnp = vsnp.astype({'POS':str})
    vsnp["Allele_{s}_vsnp".format(s=sample)] = vsnp["POS"] + vsnp['ALT']
    vsnp_allele = vsnp["Allele_{s}_vsnp".format(s=sample)]
    vsnp_allele.to_csv('{s}_vsnp.tsv'.format(s=sample),header=False,index=False, sep="\t")

files = glob.glob("*.vcf")
prefixes = [re.sub(r'.vcf','',vcf) for vcf in files]
prefixes.sort()
file_count = 0
length = len(files)
Ex=[]

while file_count < length:
    if file_count == length :
        pass
        file_count +=1
    else:
        sample = prefixes[file_count]
        print(sample)
        print('File count:{s}'.format(s=file_count))
        combined_vcfs= vcf_merge(sample)
        #Ex.append(count_list)
        file_count = file_count + 1

vsnp = glob.glob('*_vsnp.tsv')
vsnp.sort()
print(vsnp)
dfs_vsnp = [pd.read_table(vcf, header=None) for vcf in vsnp]
#dfs = [df.set_index('POS') for df in df_list]
big_vsnp = pd.concat(dfs_vsnp, axis=0)
all_vsnp = pd.concat(dfs_vsnp, axis=1, join='outer')
all_vsnp.sort_index(axis=0,ascending=True,inplace=True)
big_vsnp.sort_index(axis=0,ascending=True,inplace=True)
big_vsnp2 = big_vsnp.drop_duplicates(keep="first")
big_vsnp2.to_csv('vsnp_positions.tsv',header=['vsnp'],index=None)
big_vsnp.to_csv('vsnp_allpositions.tsv',header=['vsnp'],index=None)
all_vsnp.to_csv('all_samples_vsnp.tsv')	   
