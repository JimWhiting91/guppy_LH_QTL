#!/bin/bash

# Blast candidate genes onto the Sharma et al Transcriptome

MASTER=~/Exeter/qtl_crosses
cd $MASTER

GENOME=~/Exeter/Genomes/Poecilia_reticulata.Guppy_female_1.0_MT.dna.toplevel.fa
TRANSCRIPTOME=~/Exeter/Genomes/trin_cuff_v14_cdhit90.fa
CANDIDATES=data/chr19_candidate_genes.txt
RUN_NAME=chr19_candidate_genes

# Loop over to make fasta
rm -f outputs/$RUN_NAME.fa
for gene in {1..40}
do

# Get details
CHR=$(cat $CANDIDATES | cut -f2 | sed "${gene}q;d")
START=$(cat $CANDIDATES | cut -f3 | sed "${gene}q;d")
END=$(cat $CANDIDATES | cut -f4 | sed "${gene}q;d")

# Pipe sequence to transcriptome
samtools faidx $GENOME ${CHR}:${START}-${END} >> outputs/$RUN_NAME.fa

done

# Also align with minimap instead of blast...
minimap2 -t6 $TRANSCRIPTOME outputs/$RUN_NAME.fa > outputs/$RUN_NAME.paf
