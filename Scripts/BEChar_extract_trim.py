###############################################################################
# Modified by Gaelen Hess on 5/30/2023
# This script processes initial reads to generate trimmed fastqs with umis in the read coordiante
################################################################################
# Import neccessary modules

#Requires cutadapt and umi_tools

import pandas as pd
import numpy as np
import sys
import os
import argparse
import subprocess
import shlex
from Bio.SeqIO.QualityIO import FastqGeneralIterator

#Initiate Argument Parser
parser = argparse.ArgumentParser('Extract UMIs and Trim Read1 and Read 2 files')

parser.add_argument('parameter_file', help='Path to csv containing parameters of files', type=str)

parser.add_argument('output_name', help='Path to folder with output', type=str)

args = parser.parse_args()

#Read in Parameter File Data
FileInfo=pd.read_csv(args.parameter_file, header=0)

for index, row in FileInfo.iterrows():
  R1_file=row['Filename']+"R1.fastq.gz"
  R2_file=row['Filename']+"R2.fastq.gz"
  Out_R1_1=os.path.join(args.output_name,"Temp_R1.fastq.gz")
  Out_R2_1=os.path.join(args.output_name,"Temp_R2.fastq.gz")
  Out_R1_file=os.path.join(args.output_name,row['Output']+"_Locus.fastq.gz")
  Out_R2_file=os.path.join(args.output_name,row['Output']+"_sgRNA.fastq.gz")

###############################################################################
# This is an umi_tools command we are attempting to run
# umi_tools extract --extract-method=regex --bc-pattern=".{95,}(?P<discard_1>GGATC)(?P<umi_1>.{7})(?P<discard_2>.*)" --stdin Example_R1.fastq.gz --stdout Output_R1.fastq.gz --read2-in=Example_R2.fastq.gz --read2-out=Output_R2.fastq.gz
###############################################################################

  subprocess.call(shlex.split(" ".join(["umi_tools","extract", "--extract-method=regex",
                            "--bc-pattern=\".{95,}(?P<discard_1>GGATC)(?P<umi_1>.{7})(?P<discard_2>.*)\"", "--stdin", R1_file, "--stdout", Out_R1_1,
                            "--read2-in="+R2_file, "--read2-out="+Out_R2_1])))

###############################################################################
# This is the cutadapt command we are trying to run.
# cutadapt -l 20 -o OutputTrim_R2.fastq.gz Example_R2.fastq.gz
###############################################################################                            
                            
  subprocess.call(shlex.split(" ".join(["cutadapt","-l","20","-o",Out_R2_file,Out_R2_1])))
  
###############################################################################
# This is the cutadapt command we are trying to run.
# cutadapt -e 0.12 -g GCGCCgcttagcgccctc -o OutputTrim_R1.fastq.gz Example_R1.fastq.gz
###############################################################################
  
  subprocess.call(shlex.split(" ".join(["cutadapt","-e","0.12","-g","GCGCCgcttagcgccctc","-o",Out_R1_file, Out_R1_1])))
  
# Clean up the temporary files 
  subprocess.call(shlex.split(" ".join(["rm",Out_R1_1])))
  subprocess.call(shlex.split(" ".join(["rm",Out_R2_1])))
