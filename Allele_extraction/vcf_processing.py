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
    simulated = read('{s}.sim.vcf'.format(s=z))
    snpgenie = read('{s}.snpgenie.vcf'.format(s=z))
    vsnp = read('{s}.vsnp.vcf'.format(s=z))
    mtbseq = read('{s}.mtbseq.vcf'.format(s=z))
    bovtb = read('{s}.bovtb.vcf'.format(s=z))
    simulated = simulated[['CHROM','POS','REF','ALT']]
    #EXTRACT COLUMNS AND SWITCH CHROM NAME TO SAMPLE NAME
    simulated["CHROM"] = simulated["CHROM"].replace('NC_002945.4','{s}_sim'.format(s=sample))
    snpgenie = snpgenie[['CHROM','POS','REF','ALT']]
    snpgenie["CHROM"] = snpgenie["CHROM"].replace('NC_002945.4','{s}_snpgenie'.format(s=sample))
    vsnp = vsnp[['CHROM','POS','REF','ALT']]
    vsnp["CHROM"] = vsnp["CHROM"].replace('NC_002945.4','{s}_vsnp'.format(s=sample))
    mtbseq = mtbseq[['CHROM','POS','REF','ALT']]
    mtbseq["CHROM"] = mtbseq["CHROM"].replace('NC_002945.4','{s}_mtbseq'.format(s=sample))
    bovtb = bovtb[['CHROM','POS','REF','ALT']]
    bovtb["CHROM"] = bovtb["CHROM"].replace('NC_002945.4','{s}_bovtb'.format(s=sample))
#REMOVE LIST FROM ALT
    simulated['ALT'] = simulated['ALT'].apply(lambda x: ', '.join(map(str, x)))
    snpgenie['ALT'] = snpgenie['ALT'].apply(lambda x: ', '.join(map(str, x)))
    vsnp['ALT'] = vsnp['ALT'].apply(lambda x: ', '.join(map(str, x)))
    mtbseq['ALT'] = mtbseq['ALT'].apply(lambda x: ', '.join(map(str, x)))
    bovtb['ALT'] = bovtb['ALT'].apply(lambda x: ', '.join(map(str, x)))

#ADD PREFIXES TO COLUMN NAMES
    #simulated = simulated.rename(columns = {'REF':'{s}_sim_filt_REF'.format(s=sample),'ALT':'{s}_sim_filt_ALT'.format(s=sample)})
    #snpgenie = snpgenie.rename(columns = {'REF':'{s}_snpgenie_filt_REF'.format(s=sample),'ALT':'{s}_snpgenie_filt_ALT'.format(s=sample)})
    #vsnp = vsnp.rename(columns = {'REF':'{s}_vsnp_filt_REF'.format(s=sample),'ALT':'{s}_vsnp_filt_ALT'.format(s=sample)})

    simulated = simulated.astype({'POS':str})
    simulated["Allele_{s}_sim".format(s=sample)] = simulated["POS"] + simulated['ALT']
    simulated_allele = simulated["Allele_{s}_sim".format(s=sample)]
    
    snpgenie = snpgenie.astype({'POS':str})
    snpgenie["Allele_{s}_snpgenie".format(s=sample)] = snpgenie["POS"] + snpgenie['ALT']
    snpgenie_allele = snpgenie["Allele_{s}_snpgenie".format(s=sample)]
        
    vsnp = vsnp.astype({'POS':str})
    vsnp["Allele_{s}_vsnp".format(s=sample)] = vsnp["POS"] + vsnp['ALT']
    vsnp_allele = vsnp["Allele_{s}_vsnp".format(s=sample)]
        
    mtbseq = mtbseq.astype({'POS':str})    
    mtbseq["Allele_{s}_mtbseq".format(s=sample)] = mtbseq["POS"] + mtbseq['ALT']
    mtbseq_allele = mtbseq["Allele_{s}_mtbseq".format(s=sample)]

    bovtb = bovtb.astype({'POS':str})    
    bovtb["Allele_{s}_bovtb".format(s=sample)] = bovtb["POS"] + bovtb['ALT']
    bovtb_allele = bovtb["Allele_{s}_bovtb".format(s=sample)]    
        
    simulated_allele.to_csv('{s}_simulated.tsv'.format(s=sample),header=False,index=False, sep="\t")
    snpgenie_allele.to_csv('{s}_snpgenie.tsv'.format(s=sample),header=False,index=False, sep="\t")
    vsnp_allele.to_csv('{s}_vsnp.tsv'.format(s=sample),header=False,index=False, sep="\t")	
    mtbseq_allele.to_csv('{s}_mtbseq.tsv'.format(s=sample),header=False,index=False, sep="\t")
    bovtb_allele.to_csv('{s}_bovtb.tsv'.format(s=sample),header=False,index=False, sep="\t")
#COMBINE DATAFRAMES BASED ON POSITION
    dataframe = pd.merge(simulated, snpgenie,  how='outer', on=['POS'])
    combined_vcfs = pd.merge(dataframe, vsnp, how='outer', on=['POS'])
    combined_vcfs = pd.merge(combined_vcfs, mtbseq, how='outer', on=['POS'])
    combined_vcfs = pd.merge(combined_vcfs, bovtb, how='outer', on=['POS'])
    combined_vcfs = combined_vcfs[['REF_x','Allele_{s}_sim'.format(s=sample), 'Allele_{s}_snpgenie'.format(s=sample),'Allele_{s}_vsnp'.format(s=sample), 'Allele_{s}_mtbseq'.format(s=sample), 'Allele_{s}_bovtb'.format(s=sample)]]
    combined_vcfs.to_csv('{s}_combined.tsv'.format(s=sample),header=False,index=False, sep="\t")
#EXTRACT & COMBINE POSITIONS
    sim_pos = simulated['POS']
    sim_pos = pd.DataFrame(sim_pos)
    sim_pos.columns = ["{s}_sim".format(s=sample)]
    snpgenie_pos = snpgenie['POS']
    snpgenie_pos = pd.DataFrame(snpgenie_pos)
    snpgenie_pos.columns = ["{s}_snpgenie".format(s=sample)]
    vsnp_pos = vsnp['POS']
    vsnp_pos = pd.DataFrame(vsnp_pos)
    vsnp_pos.columns = ["{s}_vsnp".format(s=sample)]
    mtbseq_pos = mtbseq['POS']
    mtbseq_pos = pd.DataFrame(mtbseq_pos)
    mtbseq_pos.columns = ["{s}_mtbseq".format(s=sample)]
    bovtb_pos = bovtb['POS']
    bovtb_pos = pd.DataFrame(bovtb_pos)
    bovtb_pos.columns = ["{s}_bovtb".format(s=sample)]
    
    dataframe = pd.merge(sim_pos, snpgenie_pos,  left_index=True, right_index=True, how='outer')
    combined_pos = pd.merge(dataframe, vsnp_pos, left_index=True, right_index=True, how='outer')
    combined_pos = pd.merge(combined_pos, mtbseq_pos, left_index=True, right_index=True, how='outer')
    combined_pos = pd.merge(combined_pos, bovtb_pos, left_index=True, right_index=True, how='outer')
    combined_pos.to_csv('{s}_combined.tsv'.format(s=sample),header=True,index=True, sep="\t")
    #numbers = list(range(0,combined_vcfs.shape[0]))
    #combined_vcfs = compare(numbers)
    #combined_vcfs.to_csv('{s}_combined.csv'.format(s=sample),index=True)
    
    return combined_vcfs, combined_pos


files = glob.glob("*sim.vcf")
prefixes = [re.sub(r'.sim.vcf','',vcf) for vcf in files]
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

simulated = glob.glob('*_simulated.tsv')
simulated.sort()
dfs_sim = [pd.read_table(sim, header=None) for sim in simulated]
#dfs = [df.set_index('POS') for df in df_list]
big_sim = pd.concat(dfs_sim, axis=0)
all_sim = pd.concat(dfs_sim, axis=1, join='outer')
all_sim.sort_index(axis=0,ascending=True,inplace=True)
big_sim.sort_index(axis=0,ascending=True,inplace=True)
big_sim2 = big_sim.drop_duplicates(keep="first")
big_sim2.to_csv('simulated_positions.tsv',header=['simulated'],index=None)
big_sim.to_csv('simulated_allpositions.tsv',header=['simulated'],index=None)
all_sim.to_csv('all_samples_simulated.tsv')


vsnp = glob.glob('*_vsnp.tsv')
vsnp.sort()
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


snpgenie = glob.glob('*_snpgenie.tsv')
snpgenie.sort()
dfs_snpgenie = [pd.read_table(vcf, header=None) for vcf in snpgenie]
#dfs = [df.set_index('POS') for df in df_list]
big_snpgenie = pd.concat(dfs_snpgenie, axis=0)
all_snpgenie = pd.concat(dfs_snpgenie, axis=1, join='outer')
all_snpgenie.sort_index(axis=0,ascending=True,inplace=True)
big_snpgenie.sort_index(axis=0,ascending=True,inplace=True)
big_snpgenie2 = big_snpgenie.drop_duplicates(keep="first")
big_snpgenie2.to_csv('genie_positions.tsv',header=['snpgenie'],index=None)
big_snpgenie.to_csv('genie_allpositions.tsv',header=['snpgenie'],index=None)
all_snpgenie.to_csv('all_samples_snpgenie.tsv')

mtbseq = glob.glob('*_mtbseq.tsv')
mtbseq.sort()
dfs_mtbseq = [pd.read_table(vcf, header=None) for vcf in mtbseq]
#dfs = [df.set_index('POS') for df in df_list]
big_mtbseq = pd.concat(dfs_mtbseq, axis=0)
all_mtbseq = pd.concat(dfs_mtbseq, axis=1, join='outer')
all_mtbseq.sort_index(axis=0,ascending=True,inplace=True)
big_mtbseq.sort_index(axis=0,ascending=True,inplace=True)
big_mtbseq2 = big_mtbseq.drop_duplicates(keep="first")
big_mtbseq2.to_csv('mtbseq_positions.tsv',header=['mtbseq'],index=None)
big_mtbseq.to_csv('mtbseq_allpositions.tsv',header=['mtbseq'],index=None)
all_mtbseq.to_csv('all_samples_mtbseq.tsv')

bovtb = glob.glob('*_bovtb.tsv')
bovtb.sort()
dfs_bovtb = [pd.read_table(vcf,header=None) for vcf in bovtb]
#dfs = [df.set_index('POS') for df in df_list]
big_bovtb = pd.concat(dfs_bovtb, axis=0)
all_bovtb = pd.concat(dfs_bovtb, axis=1, join='outer')
all_bovtb.sort_index(axis=0,ascending=True,inplace=True)
big_bovtb.sort_index(axis=0,ascending=True,inplace=True)
big_bovtb2 = big_bovtb.drop_duplicates(keep="first")
big_bovtb2.to_csv('bovtb_positions.tsv',header=['bovtb'],index=None)
big_bovtb.to_csv('bovtb_allpositions.tsv',header=['bovtb'],index=None)
all_bovtb.to_csv('all_samples_bovtb.tsv')




