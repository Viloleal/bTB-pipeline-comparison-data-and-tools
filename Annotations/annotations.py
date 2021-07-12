#!/usr/bin/env python
# coding: utf-8
#This script extracts annotations from annotated vcf_files and shows them in a simplified manner


import vcf
import pandas as pd
import numpy as np
import sys
import gffpandas.gffpandas as gffpd
import glob
import os
import argparse
import re
import subprocess

parser = argparse.ArgumentParser(description="VCF annotation analysis")
parser.add_argument("sufix", metavar = "-P", type=str, help ="sufix to be used")
args = parser.parse_args()



def read(f):
    vcf_reader = vcf.Reader(open(f,"r"))
    sites = [record.POS for record in vcf_reader]
    sites = pd.Series(sites)
    vcf_reader = vcf.Reader(open(f,"r"))
    alt = [record.ALT for record in vcf_reader]
    alt= pd.Series(alt)
    vcf_reader = vcf.Reader(open(f,"r"))
    annot = [record.INFO['ANN'] for record in vcf_reader]
    annot = pd.Series(annot)
    return sites,alt,annot

def annot_genes(f):
    annotation = gffpd.read_gff3(f)
    gene =annotation.filter_feature_of_type(['gene'])
    gene = gene.attributes_to_columns()
    gene = gene[['ID','type','Name']]
    gene = gene.rename(columns={'Name':'name'})
    gene['product'] = 'none'
    CDS =annotation.filter_feature_of_type(['CDS'])
    CDS = CDS.attributes_to_columns()
    CDS = CDS[['locus_tag','type','gene','product']]
    CDS = CDS.rename(columns={'locus_tag':'ID','gene':'name'}) 
    CDS['ID']=CDS['ID'].str.replace('BQ','gene-BQ') #replace with gene so can be compared with vcf annotation
    return gene,CDS
    
def vcf_annotate(f):
#this function merges ID and gene product but also incorporates any other ID that is added by snpeff, such as "null" or "intergenic_regions".
    #sites,annot = read('{s}/{s}_{p}_ann.vcf'.format(s=sample, p = sufix))
    
    sites,alt,annot = read('{s}_{p}_ann.vcf'.format(s=sample, p = sufix))
    a = pd.DataFrame(columns = ['sites'],
                     data=sites)
    a = a.astype({'sites':str})
    x = pd.DataFrame(columns = ['ALT'],data=alt)
    x['ALT'] = x['ALT'].apply(lambda x: ', '.join(map(str, x)))
    a['sites'] = a['sites'] + x['ALT']
    b = pd.DataFrame(columns=['annot'],
                data=annot)
    d = a
    b['annot'] = b['annot'].apply(lambda x: ', '.join(map(str, x)))
    b = b.annot.str.split("|",expand=True,)
    c = b[[4]]
    a[['ID']]=c
    header_list = ['A','product','C','D','E','F','G','H','I','J','K','L','M','N','O','P']
    b = b.iloc[:,:16]
    b.columns = header_list
    b = b['product']
    d = d.join(b)
    d.insert(2, 'name', 'None')
    e = d[d['product'] == 'intergenic_region']
    f = d[~d['ID'].str.contains('gene', na=False)]
    g = pd.concat([e,f])
    
    return a,g    


#directories = list(filter(os.path.isdir, os.listdir()))
files = glob.glob("*ann.vcf")
directories = [re.sub(r'_\w+_ann.vcf','',vcf) for vcf in files]
file_count = 0
sufix = ref = args.sufix
length = len(directories)
df = pd.DataFrame()
while file_count < length:
    if file_count == length :
        pass
        file_count +=1
    else:
        sample = directories[file_count]
        print(sample)
        print('File count:{s}'.format(s=file_count))
        a,g = vcf_annotate(directories)
        gene, CDS = annot_genes('genes.gff')
        merged = pd.merge(a,
                          gene[['ID','name']],
                          on='ID')
        merged = pd.merge(merged,
                          CDS[['ID','product']],
                          on='ID', how='left')
        total = pd.concat([merged,g])
        total.reset_index(drop=True, inplace=True)
        total.to_csv('{s}_{p}_annotation.csv'.format(s=sample,p=sufix),header=True,index=True,index_label='{s}'.format(s=sample))
        total1 = total[['sites','ID','name','product']]	
        df = df.append(total1)
        file_count = file_count + 1
df.reset_index(drop=True, inplace=True)
df = df.sort_values(by=['sites'])
df = df.drop_duplicates()
subprocess.call(['bash', '-c', 'cat *.csv >> all_samples.csv'])
df.to_csv('total_sites.csv', header=True,index=False)

