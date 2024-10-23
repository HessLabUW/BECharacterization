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

The oligo_file is a .csv file that contains a Oligo_ID in the first column and the sequence of the oligo ordered from a chip. The sequence encodes the full oligo including parts that include handles for PCR and cloning. These parts need to be trimmed off the oligo in these scripts, so we make an alignment against just the peptide sequence.

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

Run this script first without the *-t* flag and it will show you what the actual sequences you are using to build the indice..  You can then determine the porper trims using the *-s* and *-e* flags.  Once you are sure you have the right oligos.  Run with the *-t* flag which will generate the bowtie index files in the Indices folder as well as update the screen_type_index.txt
