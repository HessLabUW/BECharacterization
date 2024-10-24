###############################################################################
# Created by Gaelen Hess on 11/10/2023
# This generates allele tables normalized by the parent samples
# This script is for the lenti dataset, but the only difference is the starting table
################################################################################

library(tidyverse)
library(data.table)
library(fs)
library(magrittr)
library(ggseqlogo)
library(ggpubr)
library(patchwork)

# Filter out columns that are not needed in analysis
alleleTable.format = fread("/project8/HessLabCHGPM/Users/ghess3/BEChar/FinalData/FastqProcessing/Data/AlleleTables/Lenti/Lenti_AlleleTable_Norm.tsv",stringsAsFactors = F, select=c(25, 26, 27, 30, 18, 20, 14, 16, 22, 23, 5, 4, 6, 3, 8, 1,24, 29, 28))

# Merge rows that have alleles that are similar within the established quantification window into a single allele and sum frequencies for all of the alleles that match
p = alleleTable.format %>% 
        select(rep,
               sample,
               guideID,
               `Norm%Reads`,
               deletion_positions,
               deletion_sizes,
               insertion_positions,
               insertion_sizes,
               substitution_positions,
               substitution_values,
               n_deleted,
               n_inserted,
               n_mutated, Reference_Sequence, Read_Status, '#Reads', '%Reads', Targeted, Strand) %>%
        group_by(rep,sample, guideID,
                 deletion_positions,deletion_sizes,
                 insertion_positions,insertion_sizes,
                 substitution_positions,substitution_values,
                 n_deleted,n_inserted,n_mutated,Reference_Sequence, Read_Status, Targeted, Strand) %>%
        summarize_at(c('Norm%Reads','#Reads','%Reads'), sum, na.rm = TRUE) %>%
        pivot_wider(id_cols= c(sample, guideID,
                 deletion_positions,deletion_sizes,
                 insertion_positions,insertion_sizes,
                 substitution_positions,substitution_values,
                 n_deleted,n_inserted,n_mutated,Reference_Sequence, Read_Status, Targeted, Strand), names_from = rep, values_from=c('Norm%Reads','#Reads','%Reads'), values_fill = 0) %>%
        rename ( NormReads_Rep1 = 'Norm%Reads_Rep1',NormReads_Rep2 = 'Norm%Reads_Rep2', NumReads_Rep1 = '#Reads_Rep1', NumReads_Rep2 = '#Reads_Rep2',PercReads_Rep1 = '%Reads_Rep1', PercReads_Rep2 = '%Reads_Rep2')

# Make a table for Parent/untrans samples and treated samples
p.untrans = subset(p, sample == "Untr")
p.trans= subset(p, sample != "Untr")

# Merge tables based on whether the allele in the trans table is in the untrans table, so they can be normalized
p.merge = merge(p.trans, p.untrans, by = names(p.untrans)[2:15], all.x = T)

names(p.merge)[15:28] = c("sample","NormReads_Rep1","NormReads_Rep2","NumReads_Rep1",
                          "NumReads_Rep2","PercReads_Rep1","PercReads_Rep2","sample_control",
                          "NormReadsCon_Rep1","NormReadsCon_Rep2","NumReadsCon_Rep1","NumReadsCon_Rep2","PercReadsCon_Rep1","PercReadsCon_Rep2")

# if no matching allele is found in the trans or untrans, then set the frequency to 0
p.merge$NormReads_Rep1[is.na(p.merge$NormReads_Rep1)] = 0
p.merge$NormReads_Rep2[is.na(p.merge$NormReads_Rep2)] = 0
p.merge$NormReadsCon_Rep1[is.na(p.merge$NormReadsCon_Rep1)] = 0
p.merge$NormReadsCon_Rep2[is.na(p.merge$NormReadsCon_Rep2)] = 0

# Calculate difference between transfected and untransfected samples
p.merge$Rep1_RawDiff = p.merge$NormReads_Rep1 - p.merge$NormReadsCon_Rep1
p.merge$Rep2_RawDiff = p.merge$NormReads_Rep2 - p.merge$NormReadsCon_Rep2

# Calculate mean across two replicates and set differences less than 0.5% to 0
p.merge$propRawDiff = rowMeans(p.merge[,c("Rep1_RawDiff","Rep2_RawDiff")],na.rm = T)
p.merge$propDiff = ifelse(p.merge$propRawDiff>0.5, p.merge$propRawDiff,0)

# Some alleles with many edits break the line in this value and this corrects this so it doesn't mess things up in downstream processing
p.merge$substitution_values = gsub("\n",",",p.merge$substitution_values, fixed = T)

# Identifies loci/target sites that have enough reads to be analyzed (greater than 200)
useful_loci = p.merge %>%
  group_by (sample, guideID) %>%
  summarize_at(c('NumReads_Rep1', 'NumReads_Rep2'), sum, na.rm = TRUE) %>%
  filter ((NumReads_Rep1+NumReads_Rep2 > 200))

# Filter table based on whether the loci had enough reads
p.merge2=semi_join(p.merge,useful_loci, by=join_by(sample, guideID))

# Outputs normalized table
write.table(p.merge2, "/project8/HessLabCHGPM/Users/ghess3/BEChar/FinalData/FastqProcessing/Data/AlleleTables/Lenti/K562_sg200_Lenti_AlleleTable_propDiff.tsv", col.names = T, row.names = F, sep = "\t", quote= F)

print("Norm File Written")

#Making the substitution files.

# Filter table for alleles that contain base substitutions
subs = p.merge2
subs = subset(subs, n_mutated >=1)

# process field of substitution positions/values to deal with each individual mutation
subs$substitution_positions = gsub("[","",subs$substitution_positions, fixed = T)
subs$substitution_positions = gsub("]","",subs$substitution_positions, fixed = T)

subs$substitution_values = gsub("[","",subs$substitution_values, fixed = T)
subs$substitution_values = gsub("]","",subs$substitution_values, fixed = T)
subs$substitution_values = gsub("\' \'",", ",subs$substitution_values, fixed = T)
subs$substitution_values = gsub("\'","",subs$substitution_values, fixed = T)

# break up each row by mutation
subs.sep = subs %>% separate_rows(substitution_positions, substitution_values, sep=", ", convert = TRUE)

# calculate the identity of the reference base by extracting from reference sequence
subs.sep$refBase = apply(subs.sep, 1, FUN = function(x) {substr(x[11], start = as.numeric(x[6])+1, stop = as.numeric(x[6])+1)})

# Output File
write.table(subs.sep,"/project8/HessLabCHGPM/Users/ghess3/BEChar/FinalData/FastqProcessing/Data/AlleleTables/Lenti/K562_sg200_Lenti_AlleleTable_Subs.tsv", col.names = T, row.names = F, sep = "\t", quote= F)

print("Subs File Written")