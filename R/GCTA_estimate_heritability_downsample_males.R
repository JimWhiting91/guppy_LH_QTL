# Estimate heritability of individuals based on GCTA approach
lib<-c("pbmcapply","pbapply","dplyr","data.table","vcfR","parallel","adegenet","qtl2","ggplot2","viridis","adegenet")
lapply(lib,library,character.only=T)

####################################################################
# Extract only the autosomes
vcf <- "data/qtl_crosses_FINAL_inds.sorted.filtered.scafs_placed_with_chroms.vcf.gz"
autosomes <- paste0("chr",1:23)
autosomes <- paste(autosomes,collapse = ",")

# Fetch the temperature and date of birth covariance data
female_rearing <- read.csv("data/female_rearing_covariates.csv",header=T)
male_rearing <- read.csv("data/male_rearing_covariates.csv",header=T)
metadata <- read.table("data/qtl_cross_ID_metadata_FINAL.tsv",header=T)

# Add in IDs
female_rearing$vcf_id <- NA
for(i in 1:nrow(female_rearing)){
  if(length(metadata[metadata$TRUE_ID == female_rearing$ID[i],"VCF_ID"]) != 0){
    female_rearing$vcf_id[i] <- metadata[metadata$TRUE_ID == female_rearing$ID[i],"VCF_ID"]
  }
}
male_rearing$vcf_id <- NA
for(i in 1:nrow(male_rearing)){
  if(length(metadata[metadata$TRUE_ID == male_rearing$ID[i],"VCF_ID"]) != 0){
    male_rearing$vcf_id[i] <- metadata[metadata$TRUE_ID == male_rearing$ID[i],"VCF_ID"]
  }
}

rearing_list <- list(female_rearing,male_rearing)
names(rearing_list) <- c("female","male")

# Make our output
dir.create("outputs/GCTA_scafs_merged_to_chrom_with_litters_male_downsample",showWarnings = F)

# Split the analysis into males and females
sex="male"

# Which indidividuals are we looking at?
cross <- read_cross2(paste0("data/qtl_cross_",sex,"s_informative.yaml"))
pedigree <- cross$covar
phenos <- data.frame(cross$pheno)

# Log-transform Days
phenos[,"Days"] <- log(phenos[,"Days"])

# Start to permute over random downsamples of the dataset to see whether we retain effects...
permN <- 100
femaleN <- 200
pblapply(1:permN,function(perm){
  set.seed(perm)
  
  # Set perm-specific output
  output <- paste0("outputs/GCTA_scafs_merged_to_chrom_with_litters_male_downsample/",sex,"_qtl_crosses_downsample_perm_",perm,"_heritability")
  
  # Sample phenos
  phenos_perm <- phenos[sample(1:nrow(phenos),femaleN),]
  
  # Make vector of ind names for VCF sampling
  all_inds <- rownames(phenos_perm)
  inds <- paste(all_inds,collapse=",")
  
  # Make autosome only vcf
  system(paste0("bcftools view -r ",autosomes," -s ",inds," ",vcf, " > ",output,"_autosomes.vcf"))
  autosome_vcf <-  paste0(output,"_autosomes.vcf")
  
  # Convert to Plink format
  system(paste0("~/bin/plink --double-id --allow-extra-chr --vcf ",autosome_vcf," --make-bed --out ",output))
  
  # Now make edits to the family files and make a phenotypes file
  # Family = column 1 of .fam
  fam <- read.table(paste0(output,".fam"))
  
  # Loop over to fill
  for(i in 1:nrow(fam)){
    fam$V1[i] <- pedigree[rownames(pedigree) == fam$V2[i],"cross"]
  }
  
  # And overwite
  write.table(fam,
              paste0(output,".fam"),row.names = F,sep = "\t",quote = F,col.names = F)
  
  # Phenotypes file...
  out_pheno <- phenos_perm
  pheno_ped <- data.frame(ind=rownames(phenos_perm))
  for(i in 1:nrow(pheno_ped)){
    pheno_ped$family[i] <- pedigree[rownames(pedigree) == pheno_ped$ind[i],"cross"]
  }
  out_pheno <- cbind(pheno_ped,out_pheno)
  out_pheno <- out_pheno[,c("family","ind",colnames(out_pheno)[3:ncol(out_pheno)])]
  
  # Write it
  write.table(out_pheno,
              paste0(output,".phen"),row.names = F,quote = F,sep="\t",col.names = F)
  
  # Estimate the Genetic Relatedness Matrix...
  system(paste0("~/bin/gcta64 --bfile ",output," --maf 0.01 --autosome-num 23 --autosome --make-grm --out ",output,"_GRM --thread-num 6"))
  
  # And Per Chromosome...
  for(i in 1:23){
    system(paste0("~/bin/gcta64 --bfile ",output," --maf 0.01 --autosome-num 23 --chr ",i," --make-grm --out ",output,"_GRM_chr",i," --thread-num 6"))
  }
  
  # Estimate additional GRM for family data...
  # Creating an additional GRM from the GRM above (setting the off-diagonals that are < 0.05 to 0)
  system(paste0("~/bin/gcta64 --grm ",output,"_GRM --make-bK 0.05 --out ",output,"_GRM_bK"))
  
  # And Per Chromosome...
  for(i in 1:23){
    system(paste0("~/bin/gcta64 --grm ",output,"_GRM_chr",i," --make-bK 0.05 --out ",output,"_GRM_chr",i,"_bK"))
  }
  
  # Get genome-wide PCs, we take 10 PCs
  system(paste0("~/bin/gcta64 --grm ",output,"_GRM --pca 10 --out ",output))
  
  # Now read these back in, confirms we need 3 PCs to separate out 4 families
  #pc_data_eigenvals <- read.table(paste0("outputs/male_qtl_crosses_heritability.eigenval"))
  pc_data <- read.table(paste0(output,".eigenvec"))
  pca_mat <- pc_data[,1:5]
  
  # Merge PCA with rearing conditions for full covar
  rearing_mat <- matrix(nrow=nrow(pca_mat),ncol=2)
  for(i in 1:nrow(pca_mat)){
    if(nrow(na.omit(rearing_list[[sex]][rearing_list[[sex]]$vcf_id == pca_mat[i,2],])) != 0){
      rearing_mat[i,1] <- na.omit(rearing_list[[sex]][rearing_list[[sex]]$vcf_id == pca_mat[i,2],"mean_temp"])
      rearing_mat[i,2] <- na.omit(rearing_list[[sex]][rearing_list[[sex]]$vcf_id == pca_mat[i,2],"DOB"])
    }
  }
  
  # Merge
  pca_rearing_mat <- cbind(pca_mat,rearing_mat)
  
  # Male share covar matrices
  # Age at maturity
  write.table(pca_rearing_mat,
              paste0(output,"_pca_rearing_qcovar_",1,".txt"),
              row.names = F,quote = F,sep = "\t",col.names = F)
  # Size at maturity
  write.table(pca_rearing_mat,
              paste0(output,"_pca_rearing_qcovar_",2,".txt"),
              row.names = F,quote = F,sep = "\t",col.names = F)
  
  # Save...
  write.table(pca_rearing_mat,
              paste0(output,"_pca_rearing_qcovar.txt"),
              row.names = F,quote = F,sep = "\t",col.names = F)
  
  # Make the combined grm + bK input
  grm_combined <- data.frame(X=c(paste0(output,"_GRM"),paste0(output,"_GRM_bK")))
  write.table(grm_combined,paste0(output,"_GRM_mgrm.txt"),
              row.names = F,quote = F,sep = "\t",col.names = F)
  
  # And per chromosome
  for(i in 1:23){
    grm_combined <- data.frame(X=c(paste0(output,"_GRM_chr",i),paste0(output,"_GRM_chr",i,"_bK")))
    write.table(grm_combined,paste0(output,"_GRM_chr",i,"_mgrm.txt"),
                row.names = F,quote = F,sep = "\t",col.names = F)
  }
  
  # Estimate variance explained for each phenotype
  for(i in 1:2){
    system(paste0("~/bin/gcta64 --mgrm ",output,"_GRM_mgrm.txt --pheno ",output,".phen --mpheno ",i," --reml --qcovar ",output,"_pca_rearing_qcovar_",i,".txt --out ",output,"_pheno_",i,"_RES --thread-num 6"))
    # And Per Chromosome...
    for(j in 1:23){
      system(paste0("~/bin/gcta64 --grm ",output,"_GRM_chr",j," --pheno ",output,".phen --mpheno ",i," --reml --qcovar ",output,"_pca_rearing_qcovar_",i,".txt --out ",output,"_pheno_",i,"_chr",j,"_RES --thread-num 6"))
    }
  }
  
  # End of permutation_set
})


##### Read back in results and summarise #####

# Phenos
pheno_vec <- c("Age at maturity (M)",
               "Size at maturity (M)")

# Read back in a summarise results from all perms and average through them
male_res <- data.frame(rbindlist(lapply(1:2,function(i){
  
  all_perm_res <- data.frame(rbindlist(pbmclapply(1:permN,function(perm){
    
    # Read in hsq file
    if(file.exists(paste0("outputs/GCTA_scafs_merged_to_chrom_with_litters_male_downsample/male_qtl_crosses_downsample_perm_",perm,"_heritability_pheno_",i,"_RES.hsq"))){
      tmp <- read.table(paste0("outputs/GCTA_scafs_merged_to_chrom_with_litters_male_downsample/male_qtl_crosses_downsample_perm_",perm,"_heritability_pheno_",i,"_RES.hsq"),fill=T,header=T)
      
      # Make output based on bK method
      out <- matrix(ncol=8,nrow=1)
      out[1,1] <- pheno_vec[i]
      out[1,2] <- as.numeric(tmp$Variance[1])
      out[1,3] <- as.numeric(tmp$Variance[2])
      out[1,4] <- as.numeric(tmp$Variance[3])
      out[1,5] <- as.numeric(tmp$Variance[4])
      out[1,6] <- as.numeric(tmp$Variance[5])
      out[1,7] <- as.numeric(tmp$Source[8])
      out[1,8] <- tmp$Variance[9]
      
      out <- data.frame(out)
      colnames(out) <- c("Phenotype","V(G1)","V(G2)","V(e)","V(P)","V(G1)/V(P)","V(G1+G2)/V(P)","logL")
      out[,2:ncol(out)] <- as.numeric(out[,2:ncol(out)])
      
      # Save 
      return(data.frame(out))
    } else {
      out <- data.frame(matrix(ncol=8,nrow=1,NA))
      colnames(out) <- c("Phenotype","V(G1)","V(G2)","V(e)","V(P)","V(G1)/V(P)","V(G1+G2)/V(P)","logL")
      return(data.frame(out))
    }
    
  },mc.cores=6)))
  
  # Now average through cols...
  perm_avg <- data.frame(Phenotype=all_perm_res$Phenotype[1],
                         t(colMeans(na.omit(all_perm_res[,2:ncol(all_perm_res)]))))
  
})))


# Merge
greml_results <- male_res
colnames(greml_results) <- c("Phenotype","V(G1)","V(G2)","V(e)","V(P)","V(G1)/V(P)","V(G1+G2)/V(P)","logL")

# Round and save these...
for(i in 2:ncol(greml_results)){
  for(j in 1:nrow(greml_results)){
    greml_results[j,i] <- round(greml_results[j,i],3)
  }
}
write.table(greml_results,
            "tables/TableSX_male_downsampled_GREML_heritability_results.txt",
            row.names = F, quote = F,sep = "\t")
