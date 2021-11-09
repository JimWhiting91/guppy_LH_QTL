#########################################################################
# Quick script to merge scaffolds and make new VCF for chr partitioning
lib<-c("dplyr","data.table","vcfR","parallel","adegenet","qtl2","ggplot2","viridis","adegenet")
lapply(lib,library,character.only=T)

# Get old vcf
vcf_path <- "data/qtl_crosses_FINAL_inds.sorted.filtered.vcf.gz"
old_vcf <- read.vcfR(vcf_path)

# Get all the chr and remove autosomes
scafs <- unique(old_vcf@fix[,1])
scafs <- scafs[!(scafs %in% paste0("chr",1:23))]

# Get linkage map info
linkage_map <- read.table("tables/TableSX_linkage_map_chromosome_position_of_unplaced_scaffolds.txt",skip = 1)
colnames(linkage_map) <- c("chr","scaf","min_cm","max_cm","relative_pos")

# Get scaf size info
fai <- read.table("~/Exeter/Genomes/STAR.chromosomes.release.fasta.fai")

# Get metadata
metadata <- old_vcf@meta
  
# Go by chromosome and make new "per chrom" vcfs
chr_vcfs <- lapply(1:23,function(x){
  
  scafs_to_add <- linkage_map[linkage_map$chr == paste0("chr",x),]
  
  tmp_chr <- old_vcf[old_vcf@fix[,1] == paste0("chr",x),]
  new_chr_gt <- tmp_chr@gt
  new_chr_fix <- tmp_chr@fix
  
  # Add the starter
  start_scafs <- scafs_to_add[scafs_to_add$relative_pos == "Start","scaf"]
  to_add <- 0
  new_fix <- matrix(ncol=8,nrow=0)
  new_gt <- matrix(ncol=ncol(new_chr_gt),nrow=0)
  for(scaf in start_scafs){
    scaf_tmp <- old_vcf[old_vcf@fix[,1] == scaf,]
    new_pos <- as.integer(scaf_tmp@fix[,2]) + to_add
    to_add <- to_add + fai[fai[,1] == scaf,2]
    scaf_fix <- scaf_tmp@fix
    scaf_fix[,2] <- as.character(new_pos)
    new_fix <- rbind(new_fix,scaf_fix)
    new_gt <- rbind(new_gt,scaf_tmp@gt)
  }
  
  # Merge these now with our main chr
  new_chr_fix[,2] <- as.character(as.integer(new_chr_fix[,2])+to_add)
  new_chr_gt2 <- rbind(new_gt,new_chr_gt)
  new_chr_fix2 <- rbind(new_fix,new_chr_fix)
  
  # Now add the ends
  to_add <- max(as.integer(new_chr_fix2[,2]))
  end_scafs <- scafs_to_add[scafs_to_add$relative_pos == "End","scaf"]
  for(scaf in end_scafs){
    scaf_tmp <- old_vcf[old_vcf@fix[,1] == scaf,]
    new_pos <- as.integer(scaf_tmp@fix[,2]) + to_add
    to_add <- to_add + fai[fai[,1] == scaf,2]
    scaf_fix <- scaf_tmp@fix
    scaf_fix[,2] <- as.character(new_pos)
    new_chr_fix2 <- rbind(new_chr_fix2,scaf_fix)
    new_chr_gt2 <- rbind(new_chr_gt2,scaf_tmp@gt)
  }
  
# Replace
  tmp_chr@gt <- new_chr_gt2
  tmp_chr@fix <- new_chr_fix2
  tmp_chr@fix[,1] <- paste0("chr",x)
  return(tmp_chr)
})

new_vcf <- chr_vcfs[[1]]
for(i in 2:23){
  new_vcf@gt <- rbind(new_vcf@gt,chr_vcfs[[i]]@gt)
  new_vcf@fix <- rbind(new_vcf@fix,chr_vcfs[[i]]@fix)
}

# Remove allocated scaffolds from metadata
scaffolds_placed <- unique(linkage_map$scaf)
to_remove <- grep("F_",metadata,value=T)
new_metadata <- metadata[!(metadata %in% to_remove)]
new_vcf@meta <- new_metadata

# Write the new VCF
write.vcf(new_vcf,"data/qtl_crosses_FINAL_inds.sorted.filtered.scafs_placed_with_chroms.vcf.gz")
