#!/bin/bash

# Script makes yaml inputs for each cross...
MASTER=/Users/jimwhiting/Exeter/qtl_crosses/
cd $MASTER

# Set em up
males=data/qtl_cross_males.yaml
females=data/qtl_cross_females.yaml

crosses=(QF2YHM2 QF4YHM8 YHF5QM6 YHF6QM7)

# Now edit each yaml
for cross in "${crosses[@]}"
do

# Make new yamls
cp $males data/qtl_cross_${cross}_males.yaml
cp $females data/qtl_cross_${cross}_females.yaml

# Need to edit the geno and the pheno
sed -i '' "s/male_phenotypes_qtl2_format/male_phenotypes_qtl2_format_${cross}/g" data/qtl_cross_${cross}_males.yaml
sed -i '' "s/female_phenotypes_qtl2_format/female_phenotypes_qtl2_format_${cross}/g" data/qtl_cross_${cross}_females.yaml

sed -i '' "s/lepmap_phased_genos_MALE_allchr/lepmap_phased_genos_MALE_allchr_${cross}/g" data/qtl_cross_${cross}_males.yaml
sed -i '' "s/lepmap_phased_genos_FEMALE_allchr/lepmap_phased_genos_FEMALE_allchr_${cross}/g" data/qtl_cross_${cross}_females.yaml

sed -i '' "s/male_qtl2_covar/male_qtl2_covar_${cross}/g" data/qtl_cross_${cross}_males.yaml
sed -i '' "s/female_qtl2_covar/female_qtl2_covar_${cross}/g" data/qtl_cross_${cross}_females.yaml

done
