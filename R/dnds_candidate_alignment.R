######################################################
# Candidate analysis for ythdc1
lib <- c("ape","seqinr","vcfR","hierfstat","adegenet","data.table","ggplot2")
lapply(lib,library,character.only=T)

### dN/dS analysis of ythdc1 alignments
# Read in alignment
aln <- read.alignment("outputs/ythdc1/ythdc1_CDS_all_poecilid_only.fas",format="fasta")
aln$nam <- c("pret",
             "pform",
             "plat",
             "pmex",
             "xiph")

# Convert
aln_dnabin <- as.DNAbin(aln)

# Calculate dnds
ythdc1_dnds <- dnds(aln_dnabin)
ythdc1_res <- as.matrix(ythdc1_dnds)

# Save this as a table
colnames(ythdc1_res) <- c("P reticulata","P formosa","P latipinna","P mexicana","Xiphophorus maculatus")
rownames(ythdc1_res) <- c("P reticulata","P formosa","P latipinna","P mexicana","Xiphophorus maculatus")
ythdc1_res <- round(ythdc1_res,3)
write.csv(ythdc1_res,
          "tables/TableSX_dnds_ythdc1_candidate_poecillids.csv",
          quote=F)
