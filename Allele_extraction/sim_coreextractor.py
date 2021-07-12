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
from functools import reduce
import allel

def read(f):
    reader = vcf.Reader(open(f))
    df = pd.DataFrame([vars(r) for r in reader])
    out = df.merge(pd.DataFrame(df.INFO.tolist()),
                   left_index=True, right_index=True)
    return out


#Load vcf files
pd.options.display.max_rows = 999
vcfs = glob.glob('*.vcf')
vcfs.sort()
prefixes = [re.sub(r'_\w+.vcf','',vcf) for vcf in vcfs]
prefixes.sort()
file_count = 0
length = len(vcfs)
Ex=[]
titles = list(prefixes)
REF = {}

#Extract all reference positions in all vcf files
while file_count < length:
    if file_count == length :
        pass
        file_count +=1
    else:
        sample = vcfs[file_count]
        prefix = prefixes[file_count]
        print(sample)
        print('File count:{s}'.format(s=file_count))
        print('selected pipeline = bovtb')
        df = allel.vcf_to_dataframe('{s}'.format(s=sample),fields=['POS','REF','ALT'])
        ref_pos = df[['POS','REF']]
        dictio =pd.Series(ref_pos.REF.values,index=ref_pos.POS).to_dict()
        file_count = file_count + 1
    REF.update(dictio)
dfref = pd.DataFrame.from_dict([REF]).T
dfref.reset_index(inplace=True)
dfref.rename(columns={'index':'POS',0:'REF'}, inplace=True)
dfref.to_csv('reference.tsv',sep="\t",header=True,index=False)



file_count = 0
length = len(vcfs)
while file_count < length:
    if file_count == length :
        pass
        file_count +=1
    else:
        sample = vcfs[file_count]
        prefix = prefixes[file_count]
        print(sample)
        print('File count:{s}'.format(s=file_count))
        df = allel.vcf_to_dataframe('{s}'.format(s=sample),fields=['POS','REF','ALT'])
        df = df.merge(dfref,how='outer',on='POS')
        df['ALT_1'] = np.where((df['ALT_1'].isna()),
                               df['REF_y'], #We place REF_x values
                               df['ALT_1'])      #In ALT_1
        #Create filter file per sample
        df_final = df[['POS','ALT_1']]
        df_final.rename(columns={"ALT_1":'{p}_ALT'.format(p=prefix)},inplace=True)
        df_final.to_csv('{p}_alt.tsv'.format(p=prefix), sep='\t', index=False, header=True)
        file_count +=1



tsvs = glob.glob('*_alt.tsv')
tsvs.sort()
prefixes = [re.sub(r'_alt.tsv','',tsv) for tsv in tsvs]
prefixes.sort()
file_count = 0
length = len(vcfs)
Ex=[]
titles = list(prefixes)
df_list = [pd.read_table(tsv) for tsv in tsvs]
dfs = [df.set_index('POS') for df in df_list]
big_df = pd.concat(dfs, axis=1, join='outer')
big_df.to_csv('all_snps.tsv', sep='\t', header=True)
#big_df = big_df.T
#big_df = big_df.drop('REF', axis=1)
big_df2 = big_df[big_df.nunique(axis=1).ne(1)]
ref_df = pd.read_table('reference.tsv', sep="\t",index_col='POS')
big_df2 = big_df2.merge(ref_df, how="left", on='POS')
big_df2 = big_df2.T
big_df2.to_csv('informative_snps.tsv',sep='\t', header=True)
indexes = list(big_df2.index.values)
row_list = []
for i in range ((big_df2.shape[0])):
        row_list.append(list(big_df2.iloc[i,:]))
        #print(row_list)
for (indx, rw) in zip (indexes,row_list):
        #print ("index: ", indx, "; row: ", rw)
        listToStr1 = ''.join(map(str, indx))
        listToStr2 = ''.join(map(str, rw))
        file = open("core.fasta", "a")
        file.write(">%s\n%s\n" % (listToStr1,listToStr2))
        file.close()

