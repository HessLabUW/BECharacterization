# BECharacterization

These scripts process fastq files to create allele tables for analyzing base editing for many sgRNA-target site pairs.

# Dependencies

The scripts use python3.7 and R version 4.4.1

The scripts were developed with the following dependencies:\
bowtie: version 1.3.1\
bowtie2: version 2.4.5\
samtools: version 0.1.19\
umi_tools: version 1.1.6\
cutadapt: version 4.9\
CRISPResso2: version 2.2.11

# Step 1: Create indices for sgRNA and target site library
*This step requires bowtie and bowtie2 to be installed*

If you have generated an index for this library previously, you don't have to recreate it.  Using "python /PATH_TO_SCRIPTS/makeCounts.py -h" will show a list of previously made libraries. You can also open the screen_type_index.txt file in the Indices folder to see the list.

The oligo_file is a .csv file that contains a Oligo_ID in the first column and the sequence of the oligo ordered from a chip. The sequence encodes the full oligo including parts that include handles for PCR and cloning. These parts need to be trimmed off the oligo in these scripts, so we make an alignment against just the sgRNA or target site sequence.

#### make_fa_sgRNA_Indices.py
>usage: Make index for alignment [-h] [-s STRIM] [-e ETRIM] [-n NUMS] [-o] [-t] [-b BOWTIE] oligo_file short_name full_name

###### positional arguments:
  >*oligo_file*            Input oligo file; csv format\
  >*short_name*            The screen type for reference\
  >*full_name*             Name for output files

###### optional arguments:
  >*-h, --help*            show this help message and exit\
  >*-s STRIM, --strim STRIM*
                        Trim bases from start; default is 0\
  >*-e ETRIM, --etrim ETRIM*
                        Trim bases from end; default is 0\
  >*-n NUMS, --nums NUMS*  Number of oligos to output (Default is 2)\
  >*-o, --override*        Flag to override existing indexes\
  >*-t, --test*            Flag to not run bowtie \
  >*-b BOWTIE, --bowtie BOWTIE*
                        Location of Bowtie aligner; default is home drive

Run this script first without the *-t* flag and it will show you what the actual sequences you are using to build the indices. You can then determine the porper trims using the *-s* and *-e* flags.  Once you are sure you have the right oligos.  Run with the *-t* flag which will generate the bowtie index files in the Indices folder as well as update the screen_type_index.txt

###### Output:

This script generates the bowtie files in the Indices folder and updates the screen_type_index.txt in the Indices folder.

#### make_fa_TargetSite_Indices_b2.py
>usage: Make index for alignment [-h] [-s STRIM] [-e ETRIM] [-n NUMS] [-o] [-t] [-b BOWTIE] oligo_file short_name full_name

###### positional arguments:
  >*oligo_file*            Input oligo file; csv format\
  >*short_name*            The screen type for reference\
  >*full_name*             Name for output files

###### optional arguments:
  >*-h, --help*            show this help message and exit\
  >*-s STRIM, --strim STRIM*
                        Trim bases from start; default is 0\
  >*-e ETRIM, --etrim ETRIM*
                        Trim bases from end; default is 0\
  >*-n NUMS, --nums NUMS*  Number of oligos to output (Default is 2)\
  >*-o, --override*        Flag to override existing indexes\
  >*-t, --test*            Flag to not run bowtie \
  >*-b BOWTIE2, --bowtie2 BOWTIE2*
                        Location of Bowtie2 aligner; default is home drive

Run this script first without the *-t* flag and it will show you what the actual sequences you are using to build the indices. You can then determine the porper trims using the *-s* and *-e* flags.  Once you are sure you have the right oligos.  Run with the *-t* flag which will generate the bowtie2 index files in the Indices folder as well as update the screen_type_index.txt

###### Output:

This script generates the bowtie2 files in the Indices folder and updates the screen_type_index.txt in the Indices folder.

# Step 2: Extract UMIs and trim sgRNA and target site reads
*This step requires umi_tools and cutadapt to be installed. This script also requires biopython package (version 1.84)*

The script requires an parameter file, which is a csv with the following columns:
>*Filename* Path and Seed for fastq file that is being used (doesn't need the R1.fastq.gz or R2.fastq.gz)\
>*Rep* Designates whether the file is from replicate 1 or 2\
>*Parent* Designates whether sample is a Parent sample (field is 1) or not (field is 0)\
>*Output* Seed name for fastq that are generated.

#### BEChar_extract_trim.py
>usage: Extract UMIs and Trim Read1 and Read 2 files [-h] parameter_file output_name

###### positional arguments:
  >*parameter_file*        Path to csv containing parameters of files\
  >*output_name*           Path to folder with output

###### optional arguments:
  >*-h, --help*            show this help message and exit

###### Output:

This script creates two fastq files per sample. The "_Locus.fastq.gz" file has the UMI extracted from Read 1 (target site) and common sequence is truncated. The "_sgRNA.fastq.gz" file has the UMI inserted in to the Read ID, and the last base is truncated off of Read 2. These files are output in the designated output folder.

# Step 3: Define list of acceptable sgRNA-target site-UMI triplets based on filtering methods
*This step requires bowtie, bowtie2, umi_tools and samtools to be installed. This script also requires biopython package (version 1.84)*

The script requires the parameter file used in Step 2.

#### Define_whitelist.py
>usage: Filter Fastqs to generate Fastqs for CRISPResso [-h] parameter_file trimmed_folder output_name

###### positional arguments:
  >*parameter_file*        Path to csv containing parameters of files\
  >*trimmed_folder*        Path to folder where UMI extracted and UMI trimmed files are located\
  >*output_name*           Path to folder with output

###### optional arguments:
  >*-h, --help*            show this help message and exit

###### Output:

For each parent/replicate sample in the parameter file, this script produces the following files:
1. "_sgRNA.map", "_sgRNA.unmapped","_Locus.map", "_Locus.unmapped" - files from the Locus and sgRNA alignments (4 files)
2. "_Aligned.fastq.gz" - fastq where unaligned reads are removed
3. "_Matched.fastq.gz" - fastq of matched and aligned reads only
4. "_IDList_Aligned.csv" and "_IDList_Matched.csv" - list of Read IDs that were maintained in the fastq (2 files)
5. "_Matched.bam", "_sorted_Matched.bam", and "_sorted_Matched.bam.bai" - alignment files needed as inputs for umitools (3 files)
6. "_Matched_UMIGroups.tsv" - List of sgRNA-target site-UMI triplets identified in sample
7. "_Matched_IndelWhitelist.csv" and "_Matched_PerfectWhitelist.csv" - These files contain all of the Locus/UMI pairs that pass the Indel or Perfect cutoff (2 files)

# Step 4: Filter fastqs generated by other samples
*This step requires bowtie, bowtie2, umi_tools and samtools to be installed. This script also requires biopython package (version 1.84)*

These scripts require the parameter file used in Step 2.

#### UMI_Filter.py
>usage: Filter Fastqs to generate Fastqs for CRISPResso [-h] [-p NUM_PROCESSORS] parameter_file trimmed_folder whitelist_folder output_name

###### positional arguments:
  >*parameter_file*        Path to csv containing parameters of files\
  >*trimmed_folder*        Path to folder where UMI extracted and UMI trimmed files are located\
  >*whitelist_folder*      Path to folder where UMI whitelists are located
  >*output_name*           Path to folder with output

###### optional arguments:
  >*-h, --help*            show this help message and exit\
  >*-p NUM_PROCESSORS, --num_processors NUM_PROCESSORS*
                        Number of processors to be used

###### Output:

For each parent/replicate sample in the parameter file, this script produces the following files:
1. "_sgRNA.map", "_sgRNA.unmapped","_Locus.map", "_Locus.unmapped" - files from the Locus and sgRNA alignments (4 files)
2. "_UMI_Matched_Indel.fastq.gz" - fastq file produced after all filtering
3. "_Matched.fastq.gz" - fastq of matched and aligned reads only
4. "_IDList_Matched.csv" and "_IDList_Matched_Indel.csv" - list of Read IDs that were maintained in the Matched and UMI_Matched_Indel fastq (2 files)
5. "_Matched.bam", "_sorted_Matched.bam", and "_sorted_Matched.bam.bai" - alignment files needed as inputs for umitools (3 files)
6. "_Matched_UMIGroups.tsv" - List of sgRNA-target site-UMI triplets identified in sample

#### UMI_Filter_ProcComp.py
>usage: Filter Fastqs to generate Fastqs for CRISPResso [-h] [-p NUM_PROCESSORS] parameter_file trimmed_folder whitelist_folder output_name

###### positional arguments:
  >*parameter_file*        Path to csv containing parameters of files\
  >*trimmed_folder*        Path to folder where UMI extracted and UMI trimmed files are located\
  >*whitelist_folder*      Path to folder where UMI whitelists are located
  >*output_name*           Path to folder with output

###### optional arguments:
  >*-h, --help*            show this help message and exit\
  >*-p NUM_PROCESSORS, --num_processors NUM_PROCESSORS*
                        Number of processors to be used

###### Output:

This script differs from UMI_Filter.py in that it produces the fastqs for all 4 filtering methods rather than just the indel version.

For each parent/replicate sample in the parameter file, this script produces the following files:
1. "_sgRNA.map", "_sgRNA.unmapped","_Locus.map", "_Locus.unmapped" - files from the Locus and sgRNA alignments (4 files)
2. "_UMI_Matched_Indel.fastq.gz" and "_UMI_Matched_Perf.fastq.gz" - fastq file produced after filtering on Indels or Perfect whitelists (2 files)
3. "_Aligned.fastq.gz" and "_Matched.fastq.gz" - fastq of aligned or matched and aligned reads (2 files)
4. "_IDList_Aligned.csv", "_IDList_Matched.csv", "_IDList_Matched_Indel.csv", and "_IDList_Matched_Perfect.csv" - list of Read IDs that were maintained in the Aligned, Matched, UMI_Matched_Indel, and UMI_Matched_Perf fastqs (4 files)
5. "_Matched.bam", "_sorted_Matched.bam", and "_sorted_Matched.bam.bai" - alignment files needed as inputs for umitools (3 files)
6. "_Matched_UMIGroups.tsv" - List of sgRNA-target site-UMI triplets identified in sample

# Step 5: Process fastqs through CRISPResso2

This is an example command run for each sample.

CRISPRessoPooled -r1 Example.fastq.gz -f AmpInfo_CRISPResso.txt --name Example_Output_Name  --exclude_bp_from_left 5 --exclude_bp_from_right 0 --min_reads_to_use_region 10 --base_editor_output --suppress_plots -w 0 --write_detailed_allele_table --output_folder Example_Output_Folder > Example_Log.log 2>&1

This uses an AmpInfo_CRISPResso.txt. This is a tab separated file with the following columns:
1. Target Site ID
2. Target Site Sequence
3. sgRNA sequence (none are listed for the scrambled sgRNAs)

# Step 6: Concatenate and normalize CRISPResso allele tables

These scripts require the parameter file used in Step 2.

#### UMI_Filter.py
>usage: Filter Fastqs to generate Fastqs for CRISPResso [-h] [-p NUM_PROCESSORS] parameter_file trimmed_folder whitelist_folder output_name

###### positional arguments:
  >*parameter_file*        Path to csv containing parameters of files\
  >*trimmed_folder*        Path to folder where UMI extracted and UMI trimmed files are located\
  >*whitelist_folder*      Path to folder where UMI whitelists are located
  >*output_name*           Path to folder with output

###### optional arguments:
  >*-h, --help*            show this help message and exit\
  >*-p NUM_PROCESSORS, --num_processors NUM_PROCESSORS*
                        Number of processors to be used

###### Output:

For each parent/replicate sample in the parameter file, this script produces the following files:
1. "_sgRNA.map", "_sgRNA.unmapped","_Locus.map", "_Locus.unmapped" - files from the Locus and sgRNA alignments (4 files)
2. "_UMI_Matched_Indel.fastq.gz" - fastq file produced after all filtering
3. "_Matched.fastq.gz" - fastq of matched and aligned reads only
4. "_IDList_Matched.csv" and "_IDList_Matched_Indel.csv" - list of Read IDs that were maintained in the Matched and UMI_Matched_Indel fastq (2 files)
5. "_Matched.bam", "_sorted_Matched.bam", and "_sorted_Matched.bam.bai" - alignment files needed as inputs for umitools (3 files)
6. "_Matched_UMIGroups.tsv" - List of sgRNA-target site-UMI triplets identified in sample