# Guppy Life History QTL
Associated scripts and data to accompany:
"On the genetic architecture of rapidly adapting and convergent life history traits in guppies"

This study investigated the genetic architecture of guppy life history traits that are known to evolve rapidly and repeatedly across natural populations under high- and low-predation in northern Trinidad. The repository includes scripts for linkage-mapping, phenotype analyses, GCTA estimates of heritability and QTL scans.

### Author Contact Information
james.whiting@ucalgary.ca
b.fraser@exeter.ac.uk

### Usage and license information
If you use or are inspired by code in this repository please cite the following work or contact me about how to cite.
Whiting et al. (2021) *bioRxiv* [doi_link](https://doi.org/10.1101/2021.03.18.435980)

---
# Environment setup
The majority of scripts in this repository can be run locally in an R environment (v4.0.3 was used).

Scripts for linkage mapping were designed to be run on a PBS-based HPC.

---
# Software
  * SNP calling was done using Stacks 2 (parameters available in manuscript) (https://catchenlab.life.illinois.edu/stacks/)
  * Linkage mapping is performed with LepMAP3 (https://sourceforge.net/p/lep-map3/wiki)
  * Heritability analyses are performed with GCTA (https://cnsgenomics.com/software/gcta)
  * QTL scans were done with r/qtl (https://rqtl.org/) and r/qtl2 (https://kbroman.org/qtl2/)

---
# Data files
 * Raw phenotype information for males and females can be found in: `/data/qtl_female_phenotypes_clean_with_litters.csv` and `qtl_male_phenotypes_clean.csv`.
 * Litter size was added later as it became available according to: `R/merge_litter_size_with_clean_females.R`.
 * SNP data = `data/qtl_crosses_FINAL_inds.sorted.filtered.vcf.gz`
 * r/qtl2 cross information are available as `data/*.yaml`
 * Full genotypes and linkage map data (including trimming guide) is available in `data/v10/*`
 * Rearing covariates (temp and DOB) = `data/*_rearing_covariates.csv`
 * Power analyses matrices = `data/*_qtl_power_matrix.csv`
 * A list of sex-biased transcripts from *Sharma et al. 2014 BMC Genomics* (doi=10.1186/1471-2164-15-400) are in `data/sharma_sex_biased_transcripts.csv`.

---
# Script information
### Data preparation/cleaning
 * `R/merge_litter_size_with_clean_females.R` - Adds on litter size information into the cleaned phenotype csv.
 * `R/merge_scaffolds_to_chroms_based_on_linkage_map.R` - Takes information from the linkage map to merge scaffolds onto the start/ends of chromosomes in the VCF for GCTA analyses.
 * `R/identify_sequencing_outliers_and_make_pedigrees_from_metadata.R` - Used to identify individuals to remove based on mislabelling and produce pedigree for linkage mapping.
 * `R/prepare_qtl2_inputs_geno_pheno_with_informative_only_with_litters.R` - Takes the linkage map genotypes and phenotype data and builds all input formats for r/qtl2 analyses (see below). Also determines which markers are fully/family-informative.
 * `R/prepare_rqtl_inputs_geno_pheno_with_informative_only.R` - Modification of the above but produces inputs in the format for r/qtl.
 * `R/plot_lepmap_vs_STAR_v10_trimming.R` - Compares linkage map to the male guppy genome and produces trimming guidelines.
 * `scripts/write_cross_specific_yaml.sh` - Loop to produce cross-specific .yaml files

### Linkage mapping
All scripts (and pedigree data) for linkage mapping are stored in the `linkage_mapping/` subdirectory.
 * `linkage_mapping/scripts/01_run_lepmap3_first_steps.sh` - Performs first few modules of LepMAP3, including: ParentCall2, Filtering2, SeparateChromosomes2, and JoinSingles2All
 * `linkage_mapping/scripts/02_run_lepmap3_OrderMarkers_NoMaleRecomb.sh` - Finalises marker order with OrderMarkers2
 * `linkage_mapping/scripts/03_rerun_map2genotypes.sh` - Final step for fetching genotypes from map files.

### Phenotype analysis (GLMs) and general figures
 * `R/phenotype_analysis.R` - Phenotype GLMs testing effects of Family, Temp, and DOB on phenotypes.
 * `R/figure1_map.R` - Produces map in figure 1 with rivers highlighted.

### GCTA analyses of heritability
 * `R/GCTA_estimate_heritability_of_phenotypes_with_litters_FINAL.R` - Main GCTA analyses, produces all GRM and does genome-wide and chr-specific tests of heritability for all phenotypes for males and females separately.
 * `R/GCTA_estimate_heritability_of_phenotypes_chr_partitioning.R` - A modification of the above, tailored specifically to chr-specific models used for the likelihood-ratio tests.
 * `R/GCTA_estimate_heritability_downsample_males.R` - Modification of the main analysis, but randomly downsamples males to the same number/power as females and repeats analyses.
 * `R/HC_correction_polygenic_test.R` - Performs accounting for Heteroscedasticity and Censoring according to *Kemppainen and Husby 2018, Evolution Letters* (https://doi.org/10.1002/evl3.88) and runs tests of associations between chromosome size and heritability.
 * `R/GCTA_power_analyses.R` - Runs power analyses for GCTA analyses.

### QTL scans
 * `R/qtl2_analysis_single_locus_final_with_litters.R` - Runs single-locus QTL scans using the R package r/qtl2 for all phenotypes and the sex locus.
 * `R/rqtl_analysis_multi-loci_final_with_litters.R` - Explores additional multi-locus QTL models with r/qtl.
 * `R/annotate_qtl_regions.R` - Annotates the QTL candidates by aligning regions to the guppy genome on Ensembl and pulling gene information from BioMart.
 * `scripts/blast_candidates_to_transcriptome.sh` - Uses minimap to align candidates to Sharma transcriptome.
 * `R/compare_qtl_candidates_to_transcriptome.R` - Compares QTL candidate genes to sex-biased transcripts in the Sharma transcriptome.
 * `R/dnds_candidate_alignment.R` - dN/dS analyses for the *ythdc1* gene candidate across poeciliids.
