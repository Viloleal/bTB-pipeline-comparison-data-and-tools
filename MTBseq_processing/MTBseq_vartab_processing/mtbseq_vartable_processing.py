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



def extract_snps(f):
    table = pd.read_table(sample,sep='\t')
    table = table.loc[table['Type'] == 'SNP']
    table = table.loc[table['Freq'] > 50]
    table = table[["#Pos", "Ref", "Allel", "Type", "CovFor", "CovRev", "Qual20", "Freq", "Cov"]]
    table.to_csv('{p}_snps.tsv'.format(p=prefix),header=True,index=False,sep='\t')
    return table
def extract_alleles(f):
    table2 = table.astype({'#Pos':str})
    table2["Allele"] = table2["#Pos"] + table2['Allel']
    table2["Allele"].to_csv('{p}_alleles.tsv'.format(p=prefix),header=False,index=False)
    

tabs = glob.glob("*.tab")
prefixes=[re.sub(r'_mtbseq_variants.tab','',tab) for tab in tabs]
tabs.sort()
prefixes.sort()
file_count = 0
length = len(tabs)
Ex=[]
while file_count < length:
    if file_count == length :
        pass
        file_count +=1
    else:
        sample = tabs[file_count]
        prefix = prefixes[file_count]
        print(sample)
        print('File count:{s}'.format(s=file_count))
        table = extract_snps(tabs)
        extract_alleles(table)
        file_count = file_count + 1


def proto_vcf(z):
    pvcf = pd.read_csv('{s}_snps.tsv'.format(s=z),sep="\t")
    pvcf = pvcf.drop(columns = ['Freq','CovFor','CovRev', 'Qual20', 'Cov', 'Type'])
    pvcf = pvcf.rename(columns = {'Pos':'POS', 'Ref':'REF', 'Allel':'ALT'})
    pvcf.insert(0, "#CHROM", "NC_002945.4")
    pvcf.insert(2, "ID", ".")
    pvcf.insert(5, "QUAL", 999)
    pvcf.insert(6, "FILTER", ".")
    pvcf.insert(7, "INFO", "DP=3552;ADF=0,2101;ADR=0,1364;AD=0,3465;VDB=1;SGB=-30.9728;MQSB=1;MQ0F=0;AC=45;AN=45;DP4=0,0,2101,1364;MQ=60")
    pvcf.insert(8, "FORMAT", "GT:PL:DP:SP:ADF:ADR:AD")
    pvcf.insert(9, "{s}".format(s=z), "1:255,0:87:0:0,71:0,16:0,87")
    return (pvcf)
def mtbseq_vcf(z):
    pvcf = pd.read_csv('{s}_snps.tsv'.format(s=z),sep="\t")
    pvcf['CovFor'] ='CovFor=' + pvcf['CovFor'].astype(str)
    pvcf['CovFor'] = pvcf['CovFor'].astype(str) + ';'
    pvcf['CovRev'] ='CovRev=' + pvcf['CovRev'].astype(str)
    pvcf['CovRev'] = pvcf['CovRev'].astype(str) + ';'
    pvcf['Qual20'] ='Qual20=' + pvcf['Qual20'].astype(str)
    pvcf['Qual20'] = pvcf['Qual20'].astype(str) + ';'
    pvcf['Cov'] ='Cov=' + pvcf['Cov'].astype(str)
    pvcf['Cov'] = pvcf['Cov'].astype(str) + ';'
    pvcf['Freq'] = 'Freq=' + pvcf['Freq'].astype(str)
    pvcf['Freq'] = pvcf['Freq'].astype(str) + ';'
    pvcf['INFO'] = pvcf['CovFor'] + pvcf['CovRev'] + pvcf['Qual20'] + pvcf['Cov'] + pvcf['Freq']
    pvcf = pvcf.drop(columns = ['Freq','CovFor','CovRev', 'Qual20', 'Cov', 'Type'])
    pvcf = pvcf.rename(columns = {'Pos':'POS', 'Ref':'REF', 'Allel':'ALT'})
    pvcf.insert(0, "#CHROM", "NC_002945.4")
    pvcf.insert(2, "ID", ".")
    pvcf.insert(5, "QUAL", 999)
    pvcf.insert(6, "FILTER", ".")
    pvcf.insert(8, "FORMAT", "GT:PL:DP:SP:ADF:ADR:AD")
    pvcf.insert(9, "{s}".format(s=z), "0:0,0:0:0:0,0:0,0:0,0")
    return (pvcf)
def create_headers(z):
    file_in = 'mtbseq_vcf_header.txt'
    file_out = '{s}_vcf_header.txt'.format(s=z)
    f_in = open(file_in, 'r')
    f_out = open(file_out, 'w') 
    for line in f_in:
        line = line.replace('3-10110', '{s}'.format(s=z))
        #print(line)
        f_out.write(line)
    f_out.close()
    f_in.close()
def append_headers(z):
    subprocess.call(['bash', '-c', 'cat {s}_vcf_header.txt {s}.vcf >> {s}.processed.vcf'.format(s=z)])

files = glob.glob("*snps.tsv")
prefixes = [re.sub(r'_snps.tsv','',tsvs) for tsvs in files]
file_count=0
length = len(files)
files.sort()
prefixes.sort()
while file_count < length:
    if file_count == length :
        pass
        file_count +=1
    else:
        sample = prefixes[file_count]
        print(sample)
        #pvcf = proto_vcf(sample)
        mtvcf = mtbseq_vcf(sample)
        #pvcf.to_csv('{s}.vcf'.format(s=sample),index=False,header=False, sep = "\t")
        mtvcf.to_csv('{s}.vcf'.format(s=sample),index=False,header=False, sep = "\t")
        create_headers(sample)
        append_headers(sample)
        file_count = file_count + 1


