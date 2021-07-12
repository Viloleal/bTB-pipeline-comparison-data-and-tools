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
    df = df.rename(columns = {'#Position':'POS'})
    df = df.rename(columns = {'Ref':'REF'})
    df = df.rename(columns = {'{s}'.format(s=ex.columns[a]):'ALT'})
    df.insert(5, "QUAL", 999)
    df.insert(6, "FILTER", ".")
    df.insert(7, "INFO", "DP=3552;ADF=0,2101;ADR=0,1364;AD=0,3465;VDB=1;SGB=-30.9728;MQSB=1;MQ0F=0;AC=45;AN=45;DP4=0,0,2101,1364;MQ=60")
    df.insert(8, "FORMAT", "GT:PL:DP:SP:ADF:ADR:AD")
    df.insert(9, '{s}'.format(s=ex.columns[a]), "1:255,0:87:0:0,71:0,16:0,87")
    df['REF'] = df['REF'].str.upper()
    df['ALT'] = df['ALT'].str.upper()
    df = df[df['REF'] != df['ALT']]
    df.to_csv('{s}.vcf'.format(s=ex.columns[a]), index=False,sep = "\t")
def create_headers(z):
    file_in = 'vcf_header.txt'
    file_out = '{s}_vcf_header.txt'.format(s=z)
    f_in = open(file_in, 'r')
    f_out = open(file_out, 'w') 
    for line in f_in:
        f_out.write(line)
    f_out.close()
    f_in.close()
def append_headers(z):
    subprocess.call(['bash', '-c', 'cat {s}_vcf_header.txt {s}.vcf >> {s}.processed.vcf'.format(s=z)])



file1 = glob.glob('*.tab')
ex = pd.read_csv(file1[0], sep="\t",skiprows=1)
ex = ex[['#Position','Ref',"13_11594","14_MBovis","15_11643","161_MBovis","17_11662","17_MBovis","182_MBovis","19_11957","19_MBovis","22_12200","23_MBovis","24_MBovis",
"25_MBovis","26_12883","26_MBovis","27_MBovis","28_12935","29_MBovis","30_MBovis","31_12952","35_MBovis","36_MBovis","37_MBovis","38_MBovis",
"39_MBovis","3_10110","41_2165","41_MBovis","42_MBovis","43_MBovis","44_MBovis","45_MBovis","47_MBovis","48_2919","48_MBovis","49_MBovis",
"50_MBovis","51_3292","51_MBovis","52_3698","54_MBovis","55_4348","56_MBovis","59_6110","59_MBovis","5_10284","7_10423"]]


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

files = glob.glob("*.vcf")
prefixes = [re.sub(r'.vcf','',vcfs) for vcfs in files]
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

