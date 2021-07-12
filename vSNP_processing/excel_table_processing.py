#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import glob
import os
import re
import subprocess


def vcf_create(x):
    a = counter
    lis = ex.iloc[:,0],ex.iloc[:,1],ex.iloc[:,a]
    df = pd.DataFrame(lis)
    df = df.transpose()
    df.insert(0,'#CHROM','NC_002945.4')
    df.insert(2, "ID", ".")
    df = df.rename(columns = {'Position':'POS'})
    df = df.rename(columns = {'root':'REF'})
    df = df.rename(columns = {'{s}'.format(s=ex.columns[a]):'ALT'})
    df.insert(5, "QUAL", 999)
    df.insert(6, "FILTER", ".")
    df.insert(7, "INFO", "DP=3552;ADF=0,2101;ADR=0,1364;AD=0,3465;VDB=1;SGB=-30.9728;MQSB=1;MQ0F=0;AC=45;AN=45;DP4=0,0,2101,1364;MQ=60")
    df.insert(8, "FORMAT", "GT:PL:DP:SP:ADF:ADR:AD")
    df.insert(9, '{s}'.format(s=ex.columns[a]), "1:255,0:87:0:0,71:0,16:0,87")
    df = df[df['REF'] != df['ALT']]
    df.to_csv('{s}.vcf'.format(s=ex.columns[a]), index=False, sep = "\t")
def create_headers(z):
    file_in = 'vsnp_header.txt'
    file_out = '{s}_vcf_header.txt'.format(s=z)
    f_in = open(file_in, 'r')
    f_out = open(file_out, 'w') 
    for line in f_in:
        f_out.write(line)
    f_out.close()
    f_in.close()
def append_headers(z):
    subprocess.call(['bash', '-c', 'cat {s}_vcf_header.txt {s}_zc.vcf >> {s}.processed.vcf'.format(s=z)])

ex = pd.read_excel('all_vcf_sort_table-2021-06-28_09-56-11.xlsx',engine = 'openpyxl',index_col=None)
ex.columns = ex.columns.str.replace(r'NC_002945.4:', '')
ex = ex.drop([51])
ex = ex.transpose()
ex = ex.reset_index()
ex.drop(columns=[])
ex.columns=ex.iloc[0]
ex = ex.rename(columns={"Unnamed: 0": "Position"})
ex = ex.drop([0])
ex = ex.drop(columns=['af2122_zc','MQ','annotations'])

counter = 2
col_len = len(ex.columns)
print("Processing excel file:")
while counter < col_len:
    if counter == col_len :
        pass
        counter +=1
    else:
        sample = ex.columns[counter]
        print(sample)
        vcf_create(sample)
        counter = counter + 1

print("Generating VCF files")

files = glob.glob("*_zc.vcf")
prefixes = [re.sub(r'_zc.vcf','',vcfs) for vcfs in files]
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
        create_headers(sample)
        append_headers(sample)
        file_count = file_count + 1

