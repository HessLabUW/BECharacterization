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

#### CreateAlleleTables_General.py
>usage: Extracts Allele Frequency Tables and concatenates them [-h] name datafolder tempfolder

###### positional arguments:
  >*name*        Path and Name for output file\
  >*datafolder*  Path for crispresso run folder\
  >*tempfolder*  Temp folder to unzip tables

###### optional arguments:
  >*-h, --help*            show this help message and exit

###### Output:

This file generates a tsv file that concatenates the allele frequency table generated by CRISPresso. The files contains the following columns:
1. #Reads
2. Aligned_Sequence
3. Reference_Sequence
4. n_inserted
5. n_deleted 
6. n_mutated 
7. Reference_Name 
8. Read_Status
9. Aligned_Reference_Names
10. Aligned_Reference_Scores
11. ref_positions
12. all_insertion_positions
13. all_insertion_left_positions
14. insertion_positions
15. insertion_coordinates
16. insertion_sizes
17. all_deletion_positions
18. deletion_positions
19. deletion_coordinates
20. deletion_sizes
21. all_substitution_positions
22. substitution_positions
23. substitution_values
24. %Reads
25. rep - Replicate of sample
26. sample - Which editor was used
27. guideID - target site id
28. Strand - does guide target top or bottom strand
29. Targeted - On Target or Scrambled

After concatenating tables, we will normalize the frequencies after removing alleles with only a single read

#### NormalizeTables_Cutoff.py
usage: NormalizeTables_Cutoff.py [-h] name infile

###### positional arguments:
  >*name*        Path and folder for output file.\
  >*infile*      Path for AlleleTable

###### optional arguments:
  >*-h, --help*            show this help message and exit

###### Output:

This file generates a new tsv with "_Norm.tsv" added to the end. It keeps all of the same columns from the concatenated table and adds an addition 'Norm%Reads' column.

#### NormalizeTables_Initial_Cutoff.py
usage: NormalizeTables_Initial_Cutoff.py [-h] name infile

###### positional arguments:
  >*name*        Path and folder for output file.\
  >*infile*      Path for AlleleTable

###### optional arguments:
  >*-h, --help*            show this help message and exit

###### Output:

This script differs from NormalizeTables_Cutoff.py in that it is used for comparing filtering methods and generates an additional column for this purpose. This file generates a new tsv with "_Norm.tsv" added to the end. It keeps all of the same columns from the concatenated table and adds 'Norm%Reads', 'treatment' and 'process_group' columns.

# Step 7: Normalize allele frequency to parent samples and make base substitution frequency table

#### MakeAllTables_Full.r and MakeAllTables_Lenti_Full.r

These scripts have no inputs and just directly read in the _Norm.tsv table generated in Step 6. They are the same scripts they just point to different input tables.

###### Output:

This file generates two tsv files _propDiff.tsv and _Subs.tsv. The first has the allele frequencies after normalization to parent sample. The Subs has the frequencies for each base substitution detected in the alleles.

_propDiff.tsv has the following columns:
1. guideID
2. deletion_positions
3. deletion_sizes
4. insertion_positions
5. insertion_sizes
6. substitution_positions
7. substitution_values
8. n_deleted
9. n_inserted
10. n_mutated
11. Reference_Sequence
12. Read_Status
13. Targeted
14. Strand
15. sample
16. NormReads_Rep1
17. NormReads_Rep2
18. NumReads_Rep1
19. NumReads_Rep2
20. PercReads_Rep1
21. PercReads_Rep2
22. sample_control
23. NormReadsCon_Rep1
24. NormReadsCon_Rep2
25. NumReadsCon_Rep1
26. NumReadsCon_Rep2
27. PercReadsCon_Rep1
28. PercReadsCon_Rep2
29. Rep1_RawDiff - frequency after normalization to Parent
30. Rep2_RawDiff - frequency after normalization to Parent
31. propRawDiff - mean frequency after normalization to Parent
32. propDiff - propRawDiff but frequencies below 0.5% are set to 0.

_Subs.tsv have the same columns as _propDiff with one additional column (#33) and changes to two other columns:
6. substitution_positions - This is a list in the _propDiff, but this is a single value for a specific substition in this table
7. substitution_values - This is a list in the _propDiff, but this is a single value for a specific substition in this table
33. refBase - This is the reference (starting) base at the position that is being mutated
