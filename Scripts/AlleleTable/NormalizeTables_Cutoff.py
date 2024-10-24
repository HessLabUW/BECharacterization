###############################################################################
# Created by Gaelen Hess on 03/10/2023
# This script filters out alleles with one read and deletes renomralizes data
# This script is used for files generated by UMI_Filter.py and not ones generated by UMI_Filter_ProcComp.py
################################################################################
#

import argparse
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests
import os
import csv
import time

# this function renormalizes the allele frequency
def norm_per_reads(row,total_read_table):
  total_read_percent=total_read_table[(total_read_table['rep']==row['rep']) & (total_read_table['sample']==row['sample']) & (total_read_table['guideID']==row['guideID'])]['%Reads'].item()
  
  norm_percent=row['%Reads']*100/total_read_percent
  
  return norm_percent

parser = argparse.ArgumentParser(description='Extracts Allele Frequency Tables and concatenates them')

# Non-optional arguments: The files containing counts, as well as an output
parser.add_argument('name', help='Path and folder for output file.', type=str)

parser.add_argument('infile', help='Path for AlleleTable', type=str)

# Saves all input to object args
args = parser.parse_args()

infile=args.infile
outfile=os.path.join(args.name,infile.split("/")[-1].split(".tsv")[0]+"_Norm.tsv")

# Read in concatenated allele table
master_table=pd.read_csv(infile, sep="\t", header=0)
print ("Data Table Loaded")

# Remove all alleles that only have 1 read supporting their existence
master_table_cntfilt=master_table[(master_table['#Reads']>1)]

# Calculate the number of reads left after removing the single read alleles
total_per_reads_table=master_table_cntfilt.groupby(['rep', 'sample', 'guideID'])['%Reads'].sum().reset_index()

# Renormalize frequency of all the alleles
master_table_cntfilt['Norm%Reads']=master_table_cntfilt.apply(lambda x: norm_per_reads(x,total_per_reads_table), axis=1)
print("Norm Reads Calculated")

# Output new normalized allele frequency table
master_table_cntfilt.to_csv(outfile, sep="\t", header=True, index=False)