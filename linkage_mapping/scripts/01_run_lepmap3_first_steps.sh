#!/bin/bash
#PBS -d .
#PBS -q pq
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=16
#PBS -A Research_Project-T110748
#PBS -N run_lepmap3
#PBS -e logs/run_lepmap3.err.txt
#PBS -o logs/run_lepmap3.out.txt
#PBS -V
#PBS -m e -M j.whiting2@exeter.ac.uk

###################################################
# Run of LepMap3 up to final step...

module load Java/1.8.0_162

MASTER=/gpfs/ts0/home/jw962/qtl_crosses/
VCF=$MASTER/data/qtl_crosses_FINAL_inds.sorted.filtered.vcf.gz
PEDIGREE=$MASTER/data/pedigrees/all_crosses_pedigree_dummyParents.txt
MASK=123
OUT=qtl_cross_full_pedigree_v10_inform${MASK}
CPU=4

cd $MASTER

############################################
# ParentCall2
# Run the first part of the pipeline, ignoring (dummy) parent order, set as half-sibs, and remove non-informative
zcat $VCF | java -cp ~/software/lepmap3/bin ParentCall2 data=$PEDIGREE \
  vcfFile=- \
  XLimit=2 \
  halfSibs=1 | gzip > outputs/${OUT}.call.gz

FILE=outputs/${OUT}.call.gz


####################################################################
# Filtering2 - Filter out data that is missing in roughly 10% of a family aside from that keep defaults...
zcat $FILE | java -cp ~/software/lepmap3/bin Filtering2 data=- \
MAFLimit=0.1 \
missingLimit=18 | gzip > outputs/${OUT}_filtered.call.gz
FILE=outputs/${OUT}_filtered.call.gz

####################################################################
# SeparateChromosomes2
zcat $FILE | java -cp ~/software/lepmap3/bin SeparateChromosomes2 data=- \
  distortionLod=1 \
  grandparentPhase=1 \
  informativeMask=$MASK \
  lodLimit = 20 \
  sizeLimit = 20 \
  numThreads=$CPU > outputs/${OUT}_lod20.map

# Observe size distribution
sort outputs/${OUT}_lod20.map | uniq -c | sort -n

# Need to split up LG1 and maybe LG2
#  zcat $FILE | java -cp ~/software/lepmap3/bin SeparateChromosomes2 data=- \
#      distortionLod=1 \
#    grandparentPhase=1 \
#    informativeMask=2 \
#    lodLimit = 30 \
#    sizeLimit = 20 \
#    map = outputs/${OUT}_lod20.map \
#    lg = 23 \
#    numThreads=$CPU > outputs/${OUT}_lod30_lg1.map
# sort outputs/${OUT}_lod30_lg1.map | uniq -c | sort -n

#
# # And again...
# zcat $FILE | java -cp ~/software/lepmap3/bin SeparateChromosomes2 data=- \
#   grandparentPhase=1 \
#   informativeMask=${MASK} \
#   lodLimit = 35 \
#   sizeLimit = 20 \
#   map = outputs/${OUT}_lod30_lg1.map \
#   lg = 1 \
#   numThreads=$CPU > outputs/${OUT}_lod30_final.map
# sort outputs/${OUT}_lod30_final.map | uniq -c | sort -r -n

################################################
# We now need to Join Singles to the map...
zcat $FILE | java -cp ~/software/lepmap3/bin JoinSingles2All data=- \
distortionLod=1 \
map=outputs/${OUT}_lod20.map \
numThreads = $CPU \
lodLimit=5 iterate=1 \
> outputs/${OUT}_lod20_final_js_iterated.map
#informativeMask=${MASK} \

cut -f 1 outputs/${OUT}_lod20_final_js_iterated.map |sort|uniq -c|sort -n
