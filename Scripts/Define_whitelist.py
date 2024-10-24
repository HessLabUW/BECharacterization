###############################################################################
# Modified by Gaelen Hess on 5/30/2023
# This script generates list of acceptable sgRNA-target site-UMI triplets that pass using differeing filtering methods.
################################################################################
# Import neccessary modules

# Requires bowtie, bowtie2, and samtools

import pandas as pd
import numpy as np
import sys
import os
import argparse
import subprocess
import shlex
from Bio.SeqIO.QualityIO import FastqGeneralIterator

sgRNA_Locus_SeqFile="/project8/HessLabCHGPM/Users/ghess3/BEChar/FinalData/FastqProcessing/LibraryDesignFiles/CXSensor_hu6_Small.csv"

sgRNA_index="/project8/HessLabCHGPM/Users/ghess3/BEChar/FinalData/FastqProcessing/Indices/sg200_sgRNA"
Locus_index="/project8/HessLabCHGPM/Users/ghess3/BEChar/FinalData/FastqProcessing/Indices/sg200_Targets"


# this function goes through a fastq file and determines whether the read matches a sgRNA-target site-UMI triplet among the list
def fastqparser(input_fastq_path, id_list, output):
  os.system("gunzip "+input_fastq_path)
  input_seed=input_fastq_path.split(".gz")[0]
  wanted = set (line.strip("\n").split(None,1)[0] for line in open(id_list))
  
  count=0
  handle = open(output,"w")
  for title, seq, qual in FastqGeneralIterator(open(input_seed)):
    if title.split(None,1)[0] in wanted:
      handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
      count +=1
  handle.close()
  print("Found and wrote "+ str(count)+" records")
  os.system("gzip "+input_seed)

# Initiate Argument Parser
parser = argparse.ArgumentParser('Filter Fastqs to generate Fastqs for CRISPResso')

parser.add_argument('parameter_file', help='Path to csv containing parameters of files', type=str)

parser.add_argument('trimmed_folder',help='Path to folder where UMI extracted and UMI trimmed files are located', type=str)

parser.add_argument('output_name', help='Path to folder with output', type=str)

args = parser.parse_args()

# Read in sequence of the oligos to determine the expected sequence of the target site
Locus_SeqInfo=pd.read_csv(sgRNA_Locus_SeqFile, header=None, names=['LocusID','FullOligoSeq'])
Locus_SeqInfo['LocusSeq']=Locus_SeqInfo.apply(lambda x: x['FullOligoSeq'][63:167],axis=1)

# Read in Parameter File Data, and only consider Parent samples
FileInfo=pd.read_csv(args.parameter_file, header=0)
FilteredFileInfo=FileInfo[FileInfo['Parent']==1]

for index, row in FilteredFileInfo.iterrows():
  
# Align sgRNA and Locus fastqs
  sgRNA_File=os.path.join(args.trimmed_folder,row['Output']+"_sgRNA.fastq.gz")
  Locus_File=os.path.join(args.trimmed_folder,row['Output']+"_Locus.fastq.gz")
  Temp_output_seed=os.path.join(args.output_name,row['Output'])
  subprocess.call(shlex.split(" ".join(["bowtie","-v", "1","--best -p", "20", "-x", sgRNA_index, "-q", sgRNA_File, Temp_output_seed+"_sgRNA.map", "--un", Temp_output_seed+"_sgRNA.unmapped"])))
  subprocess.call(shlex.split(" ".join(["bowtie2","--local -p", "20", "-x", Locus_index, "-q", Locus_File,"--no-hd -S", Temp_output_seed+"_Locus.map", "--un", Temp_output_seed+"_Locus.unmapped"])))
  print("sgRNA and Locus alligned for "+row['Output']+" sample")
  
# Read in alignment file and removed reads that are not aligned
  sgRNA_map_info=pd.read_csv(Temp_output_seed+"_sgRNA.map", header=None, sep='\t', names=['preReadID','strand','PairID','Sequence','Quality','Extra1','Extra2','Extra3'])
  sgRNA_map_info['ReadID']=sgRNA_map_info.apply(lambda x: x['preReadID'].split(" 2:N")[0],axis=1)
  sgRNA_map_info=sgRNA_map_info.filter(items=['ReadID','PairID'],axis=1)
  
  cols=['ReadID','Align_FLAG','LocusID',3,4,5,6,7,8,'Locus_Sequence',10,11,12,13,14,15,16,17,18,19]
  Locus_map_info=pd.read_csv(Temp_output_seed+"_Locus.map", header=None, sep='\t', names=cols)
  Locus_map_info=Locus_map_info.filter(items=['ReadID','Align_FLAG','LocusID','Locus_Sequence'], axis=1)
  
  AlignedLocus=Locus_map_info[Locus_map_info['LocusID']!='*']
  sgRNA_Aligned=AlignedLocus.merge(sgRNA_map_info, how='left', left_on='ReadID',right_on='ReadID')
  
# Output the list of IDs and fastq that pass the Aligned filter.
  sgRNA_Aligned['ReadID'].to_csv(Temp_output_seed+"_IDList_Aligned.csv",header=None, index=False)
  fastqparser(Locus_File, Temp_output_seed+"_IDList_Aligned.csv", Temp_output_seed+"_Aligned.fastq")
  subprocess.call(shlex.split(" ".join(["gzip",Temp_output_seed+"_Aligned.fastq"])))
  
# Apply Match filter by requiring Locus and sgRNA to be aligned th paired targets followed by outputing the list of triplets and fastq based on this filter
  sgRNA_Locus_Filter=sgRNA_Aligned[sgRNA_Aligned['LocusID']==sgRNA_Aligned['PairID']]
  sgRNA_Locus_Filter['ReadID'].to_csv(Temp_output_seed+"_IDList_Matched.csv",header=None, index=False)
  fastqparser(Locus_File, Temp_output_seed+"_IDList_Matched.csv",Temp_output_seed+"_Matched.fastq")
  subprocess.call(shlex.split(" ".join(["gzip",Temp_output_seed+"_Matched.fastq"])))
 
# Aligned the matched fastq and sorted bam files to be compatible with umi_tools
  subprocess.call(shlex.split(" ".join(["bowtie2","--local -p", "20", "-x", Locus_index, "-q", Temp_output_seed+"_Matched.fastq.gz"," | samtools view -Sb -o",Temp_output_seed+"_Matched.bam"])))
  subprocess.call(shlex.split(" ".join(["samtools","sort", Temp_output_seed+"_Matched.bam", "-o", Temp_output_seed+"_sorted_Matched.bam"])))
  subprocess.call(shlex.split(" ".join(["samtools","index", Temp_output_seed+"_sorted_Matched.bam"])))

# Use umi_tools to generate a list of all sgRNA-target site-UMI triplets observed in the sample
  subprocess.call(shlex.split(" ".join(["umi_tools","group","-I",Temp_output_seed+"_sorted_Matched.bam","--group-out="+Temp_output_seed+"_Matched_UMIGroups.tsv"])))
  UMI_Matched_read_table=pd.read_csv(Temp_output_seed+"_Matched_UMIGroups.tsv", sep="\t", header=0)

  UMI_Matched_read_table=UMI_Matched_read_table.filter(items=['read_id','final_umi'],axis=1)
  UMI_Matched=sgRNA_Locus_Filter.merge(UMI_Matched_read_table, how='left', left_on='ReadID', right_on='read_id')

# Iterate over the triplets identified and determine whether they pass the Matched-UMI-Indel or Matched-UMI-Perfect system
  mperfect_whitelist=[]
  mindel_whitelist=[]
  for name, group in UMI_Matched.groupby(['LocusID','final_umi']):
    group['Read_Length']=group.apply(lambda x: len(x['Locus_Sequence']),axis=1)
    totreads=group.shape[0]
    group_lenfilter=group[group['Read_Length']==104]
    noindelreads=group_lenfilter.shape[0]
    trueSequence=Locus_SeqInfo.loc[Locus_SeqInfo['LocusID']==name[0],'LocusSeq'].iloc[0]
    group_perfFilter=group[group['Locus_Sequence']==trueSequence]
    perfectreads=group_perfFilter.shape[0]
    if float(noindelreads)/float(totreads)>0.9: # the > 90% no indel filter required
      Locus_UMI=pd.DataFrame([[name[0],name[1]]],columns=['Locus','UMI'])
      mindel_whitelist.append(Locus_UMI)
    if float(perfectreads)/float(totreads)>0.66: # the > 66% perfect match filter required
      Locus_UMI=pd.DataFrame([[name[0],name[1]]],columns=['Locus','UMI'])
      mperfect_whitelist.append(Locus_UMI)
  print(str(len(mindel_whitelist))+" No Indel UMIs")
  print(str(len(mperfect_whitelist))+" Perfect UMIs")

# Output list of sgRNA-target site-UMI triplets that pass each filter
  mindel_whitelistfull=pd.concat(mindel_whitelist)
  mindel_whitelistfull.to_csv(Temp_output_seed+"_Matched_IndelWhitelist.csv",header=True, index=False)
  mperfect_whitelistfull=pd.concat(mperfect_whitelist)
  mperfect_whitelistfull.to_csv(Temp_output_seed+"_Matched_PerfectWhitelist.csv",header=True, index=False)