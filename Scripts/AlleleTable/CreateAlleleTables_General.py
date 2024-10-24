###############################################################################
# Created by Gaelen Hess on 03/10/2023
# This script concatenates all of the allele tables generated for each locus by CRISPRessoPooled
################################################################################
# Import neccessary modules

import argparse
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests
import os
import csv
import time
import glob
import zipfile

parser = argparse.ArgumentParser('Extracts Allele Frequency Tables and concatenates them')

# Non-optional arguments: The files containing counts, as well as an output
parser.add_argument('name', help='Path and Name for output file.', type=str)

parser.add_argument('datafolder', help='Path for crispresso run folder', type=str)

parser.add_argument('tempfolder', help ='Temp folder to unzip tables', type=str)

# Saves all input to object args
args = parser.parse_args()

ALLELE_TABLE_FOLDER=args.datafolder 
OUTPUT_SEED=args.name 
temp_seed=args.tempfolder

# Read in information on amplicons and whether they are on-target or the sgRNA targets the top or bottom strand
ampinfo_path="/project8/HessLabCHGPM/Users/ghess3/BEChar/FinalData/FastqProcessing/LibraryDesignFiles/AmpInfo_StrandTargetTable.csv"

ampinfo=pd.read_csv(ampinfo_path, header=0)
ampinfo=ampinfo.filter(items=['TargetSite','Strand','Targeted'],axis=1)
ampinfo.rename({'TargetSite' : 'guideID','Strand' : 'Strand', 'Targeted' : 'Targeted'}, axis=1, inplace=True)

alleletablelist=[]

# Iterate through all of the folders in the crispresso run folder
for name in glob.glob(os.path.join(ALLELE_TABLE_FOLDER,"*/*/Alleles_frequency_table.zip")):
  
# Unzip files and read in the data
  with zipfile.ZipFile(name, 'r') as zip:
    zip.extract("Alleles_frequency_table.txt",temp_seed)
  alleletable=pd.read_csv(os.path.join(temp_seed,"Alleles_frequency_table.txt"), sep="\t", header=0)
  
# Define where the file came from
  filebreakup=name.split(os.path.join(ALLELE_TABLE_FOLDER,"CRISPRessoPooled_on_"))[1].split("/Alleles_frequency_table.zip")[0]
  filebreakuplist=filebreakup.split("/CRISPResso_on_")
  alleletable['rep']=filebreakuplist[0].split("_")[2]
  alleletable['sample']=filebreakuplist[0].split("_")[1]
  alleletable['guideID']=filebreakuplist[1]
  alleletablelist.append(alleletable)
  
# Delete unzipped file
  os.remove(os.path.join(temp_seed,"Alleles_frequency_table.txt"))
  print("Added "+name)
totaltable=pd.concat(alleletablelist)

# Output the concatenated table
mergedtotaltable=totaltable.merge(ampinfo, how='left', on='guideID')
mergedtotaltable.to_csv(OUTPUT_SEED+"_AlleleTable.tsv", sep="\t", header=True, index=False)