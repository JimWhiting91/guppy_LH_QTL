#!/bin/bash
#PBS -d .
#PBS -q pq
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=16
#PBS -A Research_Project-T110748
#PBS -N run_lepmap3_ordermarkers
#PBS -e logs/run_lepmap3_ordermarkers.err.txt
#PBS -o logs/run_lepmap3_ordermarkers.out.txt
#PBS -V
#PBS -m e -M j.whiting2@exeter.ac.uk
#PBS -t 1-23

###################################################
# Run of LepMap3 final step...
# This runs separately for each chromosome as an arrayjob

module load Java/1.8.0_162 R/3.5.1-foss-2018b

MASTER=/gpfs/ts0/home/jw962/qtl_crosses/
VCF=$MASTER/data/qtl_crosses_FINAL_inds.sorted.filtered.vcf.gz
PEDIGREE=$MASTER/data/pedigrees/all_crosses_pedigree_dummyParents.txt
MASK=123
OUT=qtl_cross_full_pedigree_v10_inform${MASK}
CPU=16

cd $MASTER

FILE=outputs/${OUT}_filtered.call.gz
MAP=outputs/${OUT}_lod20_final_js_iterated.map

# Now OrderMarkers and run multiple iterations to take the best likelihood
zcat $FILE | java -cp ~/software/lepmap3/bin OrderMarkers2 data=- \
  map = $MAP \
  numThreads = $CPU \
  outputPhasedData=0 \
  numMergeIterations=50 \
  chromosome = $MOAB_JOBARRAYINDEX \
  usePhysical=1 0.01 \
  computeLODScores=outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_LODscores_STEP1.txt \
  sexAveraged=1 > outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_STEP1.SA.txt

# Filter markers based on LOD
Rscript R/find_erroneous_markers.R outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_LODscores_STEP1.txt outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_to_remove
grep -v "#" outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_STEP1.SA.txt | cut -f1 > outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_markers
to_remove=""
for i in $(cat outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_to_remove)
do
  marker=$(sed "${i}q;d" outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_markers)
  to_remove+="$marker "
done

# Refine (no recombination for males)
  zcat $FILE | java -cp ~/software/lepmap3/bin OrderMarkers2 data=- \
    evaluateOrder=outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_STEP1.SA.txt \
    improveOrder=0 \
    numThreads = 16 \
    grandparentPhase=1 \
    outputPhasedData=1 \
    removeMarkers=$to_remove \
    calculateIntervals=outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_calcIntervals_BEST.txt \
    computeLODScores=outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_LODscores_BEST.txt \
    sexAveraged=1 > outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_BEST.SA.txt

# And convert marker co-ords back to genome co-ords...
zcat $FILE | cut -f 1,2 |awk '(NR>=7)' > outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_snps.txt
awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_snps.txt outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_BEST.SA.txt > outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_BEST.SA.mapped.txt
sed 's/+//g' outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_BEST.SA.mapped.txt > outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_BEST.SA.mapped.clean.txt
sed -i 's/-//g' outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_BEST.SA.mapped.clean.txt
sed -i 's/_//g' outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_BEST.SA.mapped.clean.txt

# Tidy up
rm -f outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_snps.txt
rm -f outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_markers
rm -f outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_STEP1.SA.txt
rm -f outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_LODscores_STEP1.txt
rm -f outputs/${OUT}_order_LG${MOAB_JOBARRAYINDEX}_to_remove
