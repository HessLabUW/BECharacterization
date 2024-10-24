###############################################################################
# Created by Gaelen Hess on 10/17/2023
# This script filters the fastqs based on whitelist or other filters generated
# This script does all possible filtering methods compared to UMI_Filter
################################################################################
# Requires bowtie, bowtie2, umi_tools, and samtools

# Import neccessary modules

import pandas as pd
import numpy as np
import sys
import os
import argparse
import subprocess
import shlex
from Bio.SeqIO.QualityIO import FastqGeneralIterator

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

# this function checks whether the locus-umi id is prsent in the whitelist
def check_whitelist(row, whitelist_table):
  nrows=whitelist_table[(whitelist_table['Locus']== row['contig']) & (whitelist_table['UMI']==row['final_umi'])].shape[0]
  return(nrows)

# this function reads in whitelist file from before and filters the umi list for the current sample to only includes those that match the whitelist
# the list of filtered umis is written to a file as output
def whitelist_filter(whitelist_file, read_umi_list, output_seed):
  whitelist=pd.read_csv(whitelist_file, header=0)
  final_outfile=output_seed+"_IDList_"+whitelist_file.split("_")[-2]+"_"+whitelist_file.split("_")[-1].split("Whitelist")[0]+".csv"
  read_umi_list['Whitelist']=read_umi_list.apply(lambda x: check_whitelist(x, whitelist), axis=1)
  filtered_read_umi_list=read_umi_list[read_umi_list['Whitelist']!=0]
  filtered_read_umi_list['read_id'].to_csv(final_outfile, header=None, index=False)
  return final_outfile
  
#Initiate Argument Parser
parser = argparse.ArgumentParser('Filter Fastqs to generate Fastqs for CRISPResso')

parser.add_argument('parameter_file', help='Path to csv containing parameters of files', type=str)

parser.add_argument('trimmed_folder',help='Path to folder where UMI extracted and UMI trimmed files are located', type=str)

parser.add_argument('whitelist_folder',help='Path to folder where UMI whitelists are located', type=str)

parser.add_argument('output_name', help='Path to folder with output', type=str)

parser.add_argument('-p','--num_processors', help='Number of processors to be used', type=int, default=1)

args = parser.parse_args()

#Read in Parameter File Data
FileInfo=pd.read_csv(args.parameter_file, header=0)

# Iterate through sample files
for index, row in FileInfo.iterrows():

# Align sgRNA and Locus fastq files
  sgRNA_File=os.path.join(args.trimmed_folder,row['Output']+"_sgRNA.fastq.gz")
  Locus_File=os.path.join(args.trimmed_folder,row['Output']+"_Locus.fastq.gz")
  Temp_output_seed=os.path.join(args.output_name,row['Output'])
  subprocess.call(shlex.split(" ".join(["bowtie","-v", "1","--best -p", str(args.num_processors), "-x", sgRNA_index, "-q", sgRNA_File, Temp_output_seed+"_sgRNA.map", "--un", Temp_output_seed+"_sgRNA.unmapped"])))
  subprocess.call(shlex.split(" ".join(["bowtie2","--local -p", str(args.num_processors), "-x", Locus_index, "-q", Locus_File,"--no-hd -S", Temp_output_seed+"_Locus.map", "--un", Temp_output_seed+"_Locus.unmapped"])))
  print("Initial Alignment of "+row['Output']+" complete.")
  
# Filters out reads that are not aligned
  sgRNA_map_info=pd.read_csv(Temp_output_seed+"_sgRNA.map", header=None, sep='\t', names=['preReadID','strand','PairID','Sequence','Quality','Extra1',6,7])
  sgRNA_map_info['ReadID']=sgRNA_map_info.apply(lambda x: x['preReadID'].split(" 2:N")[0],axis=1)
  cols=['ReadID','Align_FLAG','LocusID',3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
  Locus_map_info=pd.read_csv(Temp_output_seed+"_Locus.map", header=None, sep='\t', names=cols)
  AlignedLocus=Locus_map_info[Locus_map_info['LocusID']!='*']
  
# Output Aligned_fastq
  sgRNA_Aligned=AlignedLocus.merge(sgRNA_map_info, how='left', left_on='ReadID',right_on='ReadID')
  sgRNA_Aligned['ReadID'].to_csv(Temp_output_seed+"_IDList_Aligned.csv",header=None, index=False)
  fastqparser(Locus_File, Temp_output_seed+"_IDList_Aligned.csv", Temp_output_seed+"_Aligned.fastq")
  subprocess.call(shlex.split(" ".join(["gzip",Temp_output_seed+"_Aligned.fastq"])))
  print(Temp_output_seed+"_Aligned.fastq.gz file generated")
  
# Applies Matched filter between Locus (target site) and sgRNA IDs needing to match
# Outputs Matched fastq.gz file here.
  sgRNA_Locus_Filter=sgRNA_Aligned[sgRNA_Aligned['LocusID']==sgRNA_Aligned['PairID']]
  sgRNA_Locus_Filter['ReadID'].to_csv(Temp_output_seed+"_IDList_Matched.csv",header=None, index=False)
  fastqparser(Locus_File, Temp_output_seed+"_IDList_Matched.csv",Temp_output_seed+"_Matched.fastq")
  subprocess.call(shlex.split(" ".join(["gzip",Temp_output_seed+"_Matched.fastq"])))
  print(Temp_output_seed+"_Matched.fastq.gz file generated")

  whitelist_seed=os.path.join(args.whitelist_folder,FileInfo[(FileInfo['Rep']==row['Rep']) & (FileInfo['Parent']==1)]['Output'].item())
  mch_indel_wl=whitelist_seed+"_Matched_IndelWhitelist.csv"
  mch_Perf_wl=whitelist_seed+"_Matched_PerfectWhitelist.csv"

# Aligned the matched fastq and sorted bam files to be compatible with umi_tools  
  subprocess.call(shlex.split(" ".join(["bowtie2","--local -p", str(args.num_processors), "-x", Locus_index, "-q", Temp_output_seed+"_Matched.fastq.gz"," | samtools view -Sb -o",Temp_output_seed+"_Matched.bam"])))
  subprocess.call(shlex.split(" ".join(["samtools","sort", Temp_output_seed+"_Matched.bam", "-o", Temp_output_seed+"_sorted_Matched.bam"])))
  subprocess.call(shlex.split(" ".join(["samtools","index", Temp_output_seed+"_sorted_Matched.bam"])))

# Use umi_tools to generate a list of all sgRNA-target site-UMI triplets observed in the sample
  subprocess.call(shlex.split(" ".join(["umi_tools","group","-I",Temp_output_seed+"_sorted_Matched.bam","--group-out="+Temp_output_seed+"_Matched_UMIGroups.tsv"])))
  UMI_Matched_read_table=pd.read_csv(Temp_output_seed+"_Matched_UMIGroups.tsv", sep="\t", header=0)

# filter list of Locus-UMIs found by umi_tools for those that match the whitelist followed by filtering fastq based on these matches and the fastq is produced as output
# This section does the filtering for the Indel Whitelist
  mch_indel_wl_file=whitelist_filter(mch_indel_wl, UMI_Matched_read_table,Temp_output_seed)
  fastqparser(Locus_File, mch_indel_wl_file, Temp_output_seed+"_UMI_Matched_Indel.fastq")
  subprocess.call(shlex.split(" ".join(["gzip",Temp_output_seed+"_UMI_Matched_Indel.fastq"])))
  print(Temp_output_seed+"_UMI_Matched_Indel.fastq.gz file generated")
  
# Repeat filtering step above but for Perfect Whitelist
  mch_Perf_wl_file=whitelist_filter(mch_Perf_wl, UMI_Matched_read_table,Temp_output_seed)
  fastqparser(Locus_File, mch_Perf_wl_file, Temp_output_seed+"_UMI_Matched_Perf.fastq")
  subprocess.call(shlex.split(" ".join(["gzip",Temp_output_seed+"_UMI_Matched_Perf.fastq"])))
  print(Temp_output_seed+"_UMI_Matched_Perf.fastq.gz file generated")
