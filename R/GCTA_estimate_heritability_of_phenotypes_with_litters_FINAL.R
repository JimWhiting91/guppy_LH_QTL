# Estimate heritability of individuals based on GCTA approach
lib<-c("dplyr","data.table","vcfR","parallel","adegenet","qtl2","ggplot2","viridis","adegenet")
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
dir.create("outputs/GCTA_scafs_merged_to_chrom_with_litters",showWarnings = F)

# Split the analysis into males and females
lapply(c("male","female"),function(sex){
  
  # Set sex-specific output...
  output <- paste0("outputs/GCTA_scafs_merged_to_chrom_with_litters/",sex,"_qtl_crosses_heritability")
  
  # Which indidividuals are we looking at?
  cross <- read_cross2(paste0("data/qtl_cross_",sex,"s_informative.yaml"))
  pedigree <- cross$covar
  phenos <- data.frame(cross$pheno)
  
  # Log-transform our female variables
  if(sex=="female"){
    
    # Filter out bad interbroods...
    phenos[phenos$interbrood < 15 | 
             phenos$interbrood > 40, c("interbrood",
                                       "first_brood_size",
                                       "dry_offspring_weight_1",
                                       "dry_offspring_weight_2")]<-NA
    
    # Remove dry offspring weight 2 and 3
    phenos <- phenos[,1:5]
    
    # Get residual first brood size
    size_mod <- lm(first_brood_size~size_at_first_brood,data=phenos)
    for(i in 1:nrow(phenos)){
      phenos$predicted_first_brood_size[i] <- predict(size_mod)[rownames(phenos)[i]]
    }
    phenos$first_brood_size_resid <- phenos$first_brood_size - as.numeric(phenos$predicted_first_brood_size)

    # Log-transform all values, test, check hists
    to_transform <- c("age_first_brood","size_at_first_brood","interbrood","dry_offspring_weight_1")
    for(i in to_transform){
      phenos[,i] <- log(phenos[,i])
    }
    
    phenos <- phenos[,c(c("age_first_brood","size_at_first_brood","first_brood_size_resid","interbrood","dry_offspring_weight_1"))]
  } else {
    phenos[,"Days"] <- log(phenos[,"Days"])
  }
  
  inds <- paste(rownames(phenos),collapse=",")
  
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
  out_pheno <- phenos
  pheno_ped <- data.frame(ind=rownames(phenos))
  for(i in 1:nrow(pheno_ped)){
    pheno_ped$family[i] <- pedigree[rownames(pedigree) == pheno_ped$ind[i],"cross"]
  }
  out_pheno <- cbind(pheno_ped,out_pheno)
  out_pheno <- out_pheno[,c("family","ind",colnames(out_pheno)[3:ncol(out_pheno)])]
  
  # Write it
  write.table(out_pheno,
              paste0(output,".phen"),row.names = F,quote = F,sep="\t",col.names = F)
  
  # Estimate the Genetic Relatedness Matrix...
  system(paste0("~/bin/gcta64 --make-grm-gz --bfile ",output," --maf 0.01 --autosome-num 23 --autosome --make-grm --out ",output,"_GRM --thread-num 6"))
  
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
  
  # Break these down into female specific covar matrices, males can share
  if(sex == "female"){
    female_phenos <- colnames(phenos)
    
    # Age first brood
    write.table(pca_rearing_mat[,-7],
                paste0(output,"_pca_rearing_qcovar_",1,".txt"),
                row.names = F,quote = F,sep = "\t",col.names = F)
    # Size at first brood
    write.table(pca_rearing_mat[,-c(6,7)],
                paste0(output,"_pca_rearing_qcovar_",2,".txt"),
                row.names = F,quote = F,sep = "\t",col.names = F)
    # First brood size
    write.table(pca_rearing_mat[,-c(6,7)],
                paste0(output,"_pca_rearing_qcovar_",3,".txt"),
                row.names = F,quote = F,sep = "\t",col.names = F)
    # Interbrood
    write.table(pca_rearing_mat,
                paste0(output,"_pca_rearing_qcovar_",4,".txt"),
                row.names = F,quote = F,sep = "\t",col.names = F)
    # Dry Offspring Weight 1
    write.table(pca_rearing_mat[,-c(6,7)],
                paste0(output,"_pca_rearing_qcovar_",5,".txt"),
                row.names = F,quote = F,sep = "\t",col.names = F)
  } else {
    # Age at maturity
    write.table(pca_rearing_mat,
                paste0(output,"_pca_rearing_qcovar_",1,".txt"),
                row.names = F,quote = F,sep = "\t",col.names = F)
    # Size at maturity
    write.table(pca_rearing_mat,
                paste0(output,"_pca_rearing_qcovar_",2,".txt"),
                row.names = F,quote = F,sep = "\t",col.names = F)
  }
  
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
  if(sex == "female"){
    for(i in 1:5){
      system(paste0("~/bin/gcta64 --mgrm ",output,"_GRM_mgrm.txt --pheno ",output,".phen --mpheno ",i," --reml --qcovar ",output,"_pca_rearing_qcovar_",i,".txt --out ",output,"_pheno_",i,"_RES --thread-num 6"))
      # And Per Chromosome...
      for(j in 1:23){
        system(paste0("~/bin/gcta64 --grm ",output,"_GRM_chr",j," --pheno ",output,".phen --mpheno ",i," --reml --qcovar ",output,"_pca_rearing_qcovar_",i,".txt --out ",output,"_pheno_",i,"_chr",j,"_RES --thread-num 6"))
      }
    }
  }  else {
    for(i in 1:2){
      system(paste0("~/bin/gcta64 --mgrm ",output,"_GRM_mgrm.txt --pheno ",output,".phen --mpheno ",i," --reml --qcovar ",output,"_pca_rearing_qcovar_",i,".txt --out ",output,"_pheno_",i,"_RES --thread-num 6"))
      # And Per Chromosome...
      for(j in 1:23){
        system(paste0("~/bin/gcta64 --grm ",output,"_GRM_chr",j," --pheno ",output,".phen --mpheno ",i," --reml --qcovar ",output,"_pca_rearing_qcovar_",i,".txt --out ",output,"_pheno_",i,"_chr",j,"_RES --thread-num 6"))
      }
    }
  }
})

##### Read back in results and summarise #####

# Phenos
pheno_vec <- c("Age at first brood (F)",
               "Size at first brood (F)",
               "First brood size (F)",
               "Interbrood period (F)",
               "First brood offspring weight (F)",
               "Age at maturity (M)",
               "Size at maturity (M)")

# Read back in a summarise results
female_res <- data.frame(rbindlist(lapply(1:5,function(i){

  # Read in hsq file
  tmp <- read.table(paste0("outputs/GCTA_scafs_merged_to_chrom_with_litters/female_qtl_crosses_heritability_pheno_",i,"_RES.hsq"),fill=T,header=T)
  
  # Make output based on bK method
  out <- matrix(ncol=8,nrow=1)
  out[1,1] <- pheno_vec[i]
  out[1,2] <- paste0(round(as.numeric(tmp$Variance[1]),2)," ± ",round(as.numeric(tmp$SE[1]),2))
  out[1,3] <- paste0(round(as.numeric(tmp$Variance[2]),2)," ± ",round(as.numeric(tmp$SE[2]),2))
  out[1,4] <- paste0(round(as.numeric(tmp$Variance[3]),2)," ± ",round(as.numeric(tmp$SE[3]),2))
  out[1,5] <- paste0(round(as.numeric(tmp$Variance[4]),2)," ± ",round(as.numeric(tmp$SE[4]),2))
  out[1,6] <- paste0(round(as.numeric(tmp$Variance[5]),2)," ± ",round(as.numeric(tmp$SE[5]),2))
  out[1,7] <- paste0(round(as.numeric(tmp$Source[8]),3)," ± ",round(as.numeric(tmp$Variance[8]),2))
  out[1,8] <- tmp$Variance[9]
  
  # Save 
  return(data.frame(out))
})))

# Read back in a summarise results
male_res <- data.frame(rbindlist(lapply(1:2,function(i){
  # Read in hsq file
  tmp <- read.table(paste0("outputs/GCTA_scafs_merged_to_chrom_with_litters/male_qtl_crosses_heritability_pheno_",i,"_RES.hsq"),fill=T,header=T)
  
  # Make output based on bK method
  out <- matrix(ncol=8,nrow=1)
  out[1,1] <- pheno_vec[i+5]
  out[1,2] <- paste0(round(as.numeric(tmp$Variance[1]),3)," ± ",round(as.numeric(tmp$SE[1]),2))
  out[1,3] <- paste0(round(as.numeric(tmp$Variance[2]),3)," ± ",round(as.numeric(tmp$SE[2]),2))
  out[1,4] <- paste0(round(as.numeric(tmp$Variance[3]),3)," ± ",round(as.numeric(tmp$SE[3]),2))
  out[1,5] <- paste0(round(as.numeric(tmp$Variance[4]),3)," ± ",round(as.numeric(tmp$SE[4]),2))
  out[1,6] <- paste0(round(as.numeric(tmp$Variance[5]),3)," ± ",round(as.numeric(tmp$SE[5]),2))
  out[1,7] <- paste0(round(as.numeric(tmp$Source[8]),3)," ± ",round(as.numeric(tmp$Variance[8]),2))
  out[1,8] <- tmp$Variance[9]
  
  # Save 
  return(data.frame(out))
})))

# Merge
greml_results <- rbind(female_res,male_res)
colnames(greml_results) <- c("Phenotype","V(G1)","V(G2)","V(e)","V(P)","V(G1)/V(P)","V(G1+G2)/V(P)","logL")

# Save these
write.table(greml_results,
            "tables/TableX_GREML_heritability_results.txt",
            row.names = F, quote = F,sep = "\t")
write.csv(greml_results,
            "tables/TableX_GREML_heritability_results.csv",
            row.names = F, quote = F)

###############################
# Also read back in the per-chromosome estimates...
female_res <- data.frame(rbindlist(lapply(1:5,function(i){
  
  # Get res per chromosome
  chr_merge <- data.frame(rbindlist(lapply(1:23,function(chr){
    print(chr)
    
    # Read in hsq file
    tmp <- read.table(paste0("outputs/GCTA_scafs_merged_to_chrom_with_litters/female_qtl_crosses_heritability_pheno_",i,"_chr",chr,"_RES.hsq"),fill=T,header=T)
    
    # Make output
    out <- matrix(ncol=7,nrow=1)
    out[1,1] <- pheno_vec[i]
    out[1,2] <- tmp$Variance[1]
    out[1,3] <- tmp$Variance[2]
    out[1,4] <- tmp$Variance[3]
    out[1,5] <- tmp$Variance[4]
    out[1,6] <- tmp$Variance[9]
    out[1,7] <- paste0("chr",chr)
    return(data.frame(out))
  })))
  
  # Save 
  return(chr_merge)
})))

male_res <- data.frame(rbindlist(lapply(1:2,function(i){
  
  # Get res per chromosome
  chr_merge <- data.frame(rbindlist(lapply(1:23,function(chr){
    
    # Read in hsq file
    tmp <- read.table(paste0("outputs/GCTA_scafs_merged_to_chrom_with_litters/male_qtl_crosses_heritability_pheno_",i,"_chr",chr,"_RES.hsq"),fill=T,header=T)
    
    # Make output
    out <- matrix(ncol=7,nrow=1)
    out[1,1] <- pheno_vec[i+5]
    out[1,2] <- tmp$Variance[1]
    out[1,3] <- tmp$Variance[2]
    out[1,4] <- tmp$Variance[3]
    out[1,5] <- tmp$Variance[4]
    out[1,6] <- tmp$Variance[9]
    out[1,7] <- paste0("chr",chr)
    return(data.frame(out))
  })))
  
  # Save 
  return(chr_merge)
})))

# Per chr merge
greml_results_chr <- rbind(female_res,male_res)
colnames(greml_results_chr) <- c("Phenotype","V(G)","V(e)","V(P)","V(G)/V(P)","p_val","chr")

# Plot
greml_results_chr$chr_F <- factor(greml_results_chr$chr,levels=paste0("chr",1:23))
greml_results_chr$Phenotype_F <- factor(greml_results_chr$Phenotype,levels=c("Age at first brood (F)",
                                                                             "Size at first brood (F)",
                                                                             "First brood size (F)",
                                                                             "Interbrood period (F)",
                                                                             "First brood offspring weight (F)",
                                                                            # "Second brood offspring weight (F)",
                                                                             "Age at maturity (M)",
                                                                             "Size at maturity (M)"))
greml_results_chr$`V(G)/V(P)` <- as.numeric(greml_results_chr$`V(G)/V(P)`)
greml_results_chr$tile_fill <- greml_results_chr$`V(G)/V(P)`

# Correct pvals
phenos <- unique(greml_results_chr$Phenotype)
greml_results_chr$fdr <- NA
for(pheno in phenos){
  greml_results_chr[greml_results_chr$Phenotype == pheno,"fdr"] <- p.adjust(greml_results_chr[greml_results_chr$Phenotype == pheno,"p_val"],method="fdr")
}

greml_results_chr$fdr_label <- NA
greml_results_chr[greml_results_chr$fdr <= 0.05,"fdr_label"] <- "*"
greml_results_chr[greml_results_chr$fdr <= 0.01,"fdr_label"] <- "**"
greml_results_chr[greml_results_chr$fdr <= 0.001,"fdr_label"] <- "***"

# Plot per chr heritability
per_chr_heritability <- ggplot(greml_results_chr,aes(x=chr_F,y=Phenotype_F,fill=tile_fill))+
  geom_tile()+
  scale_fill_viridis(option = "A")+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=14,angle=45,hjust=1),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=16),
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        panel.grid = element_blank())+
  labs(y="Phenotype",fill=expression(V[G]/V[P]))+
  geom_text(aes(x=chr_F,y=Phenotype_F,label=fdr_label),size=5,colour="black")

pdf("figs/FigureX_perchr_heritability_GCTA_scafs_merged_to_chroms.pdf",width=12,height=3)
per_chr_heritability
dev.off()

saveRDS(per_chr_heritability,
        "figs/FigureX_perchr_heritability_GCTA_scafs_merged_to_chroms.rds")

na.omit(greml_results_chr)

# In each case, let's sum the proportion of variance....
greml_results_chr %>% group_by(Phenotype) %>% summarise(sum_var=sum(`V(G)/V(P)`))

# Save per chromosome heritabilies somewhere
greml_results_chr_out <- greml_results_chr[,c(1,7,2:6,11)]
for(i in 3:8){
  greml_results_chr_out[,i] <- round(as.numeric(greml_results_chr_out[,i]),3)
}
write.csv(greml_results_chr_out[order(greml_results_chr_out$fdr),],
            "tables/TableSX_single_chr_estimates_h2c.csv",
            row.names = F,quote = F)

