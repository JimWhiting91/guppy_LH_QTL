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
dir.create("outputs/GCTA_chr_partitioning_with_litters")

# Split the analysis into males and females
lapply(c("male","female"),function(sex){
  
  # Set sex-specific output...
  output <- paste0("outputs/GCTA_chr_partitioning_with_litters/",sex,"_qtl_crosses_heritability")
  
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
    phenos$resid_first_brood_size <- phenos$first_brood_size-phenos$predicted_first_brood_size
    
    # Remove dry offspring weight 2 and 3
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
  
  # Estimate the Genetic Relatedness Matrix for all autosomes...
  system(paste0("~/bin/gcta64 --bfile ",output," --maf 0.01 --autosome-num 23 --autosome --make-grm --out ",output,"_GRM --thread-num 6"))
  
  # And Per Chromosome, estimate a GRM for all chr except focal and one for focal only
  for(i in 1:23){
    message(paste0("Starting Chr",i))
    ####################
    # list chr except for focal
    chrs <- paste0("chr",1:23)
    chrs <- chrs[chrs != paste0("chr",i)]
    chrs_filter <- paste(chrs,collapse = ",")
    
    # Make vcf for all chr minus focal
    system(paste0("bcftools view -r ",chrs_filter," -s ",inds," ",vcf, " > ",output,"_NoFocalChr_autosomes.vcf"))
    nofocal_vcf <-  paste0(output,"_NoFocalChr_autosomes.vcf")
    # Convert to Plink format
    system(paste0("~/bin/plink --double-id --allow-extra-chr --vcf ",nofocal_vcf," --make-bed --out ",output,"_NoFocalChr"))
    
    # Now make edits to the family files and make a phenotypes file
    # Family = column 1 of .fam
    fam <- read.table(paste0(output,"_NoFocalChr.fam"))
    
    # Loop over to fill
    for(j in 1:nrow(fam)){
      fam$V1[j] <- pedigree[rownames(pedigree) == fam$V2[j],"cross"]
    }
    
    # And overwite
    write.table(fam,
                paste0(output,"_NoFocalChr.fam"),row.names = F,sep = "\t",quote = F,col.names = F)
    
    
    # Make GRM
    system(paste0("~/bin/gcta64 --bfile ",output,"_NoFocalChr --maf 0.01 --autosome-num 22 --autosome --make-grm --out ",output,"_GRM_NoFocalChr_chr",i," --thread-num 6"))
    # Tidy
    system(paste0("rm -f ",nofocal_vcf))
    
    ####################
    # Focal chr
    system(paste0("~/bin/gcta64 --bfile ",output," --maf 0.01 --autosome-num 23 --chr ",i," --make-grm --out ",output,"_GRM_chr",i," --thread-num 6"))
  }
  
  # Estimate additional GRM for family data...
  # Creating an additional GRM from the GRM above (setting the off-diagonals that are < 0.05 to 0)
  system(paste0("~/bin/gcta64 --grm ",output,"_GRM --make-bK 0.05 --out ",output,"_GRM_bK"))
  
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
  }
  
  # Save...
  write.table(pca_rearing_mat,
              paste0(output,"_pca_rearing_qcovar.txt"),
              row.names = F,quote = F,sep = "\t",col.names = F)
  
  # Make the combined grm + bK input
  grm_combined <- data.frame(X=c(paste0(output,"_GRM"),paste0(output,"_GRM_bK")))
  write.table(grm_combined,paste0(output,"_GRM_mgrm.txt"),
              row.names = F,quote = F,sep = "\t",col.names = F)
  
  # And per chromosome, here we want to test 4 models, so we need to make 4 inputs reflecting the 4 model types
  for(i in 1:23){
    mod2_grm_combined <- data.frame(X=c(paste0(output,"_GRM_NoFocalChr_chr",i),paste0(output,"_GRM_chr",i)))
    mod4_grm_combined <- data.frame(X=c(paste0(output,"_GRM"),paste0(output,"_GRM_chr",i)))
    
    write.table(mod2_grm_combined,paste0(output,"_GRM_chr",i,"_mgrm_model2.txt"),
                row.names = F,quote = F,sep = "\t",col.names = F)
    write.table(mod4_grm_combined,paste0(output,"_GRM_chr",i,"_mgrm_model4.txt"),
                row.names = F,quote = F,sep = "\t",col.names = F)
  }
  
  # Estimate variance explained for each phenotype
  if(sex == "female"){
    for(i in 1:5){
      system(paste0("~/bin/gcta64 --mgrm ",output,"_GRM_mgrm.txt --pheno ",output,".phen --mpheno ",i," --reml --qcovar ",output,"_pca_rearing_qcovar_",i,".txt --out ",output,"_pheno_",i,"_RES --thread-num 6"))
      # And Per Chromosome
      # We run 4 models here as described above
      for(j in 1:23){
        # Model 1 = GRM excluding chr
        system(paste0("~/bin/gcta64 --grm ",output,"_GRM_NoFocalChr_chr",j," --pheno ",output,".phen --mpheno ",i," --reml --qcovar ",output,"_pca_rearing_qcovar_",i,".txt --out ",output,"_pheno_",i,"_chr",j,"_model1_RES --thread-num 6"))
        # Model 2 = GRM excluding chr and focal only
        system(paste0("~/bin/gcta64 --mgrm ",output,"_GRM_chr",j,"_mgrm_model2.txt --pheno ",output,".phen --mpheno ",i," --reml --qcovar ",output,"_pca_rearing_qcovar_",i,".txt --out ",output,"_pheno_",i,"_chr",j,"_model2_RES --thread-num 6"))
        # Model 3 = GRM all chr
        system(paste0("~/bin/gcta64 --grm ",output,"_GRM --pheno ",output,".phen --mpheno ",i," --reml --qcovar ",output,"_pca_rearing_qcovar_",i,".txt --out ",output,"_pheno_",i,"_chr",j,"_model3_RES --thread-num 6"))
        # Model 4 = GRM all chr + GRM focal chr
        system(paste0("~/bin/gcta64 --mgrm ",output,"_GRM_chr",j,"_mgrm_model4.txt --pheno ",output,".phen --mpheno ",i," --reml --qcovar ",output,"_pca_rearing_qcovar_",i,".txt --out ",output,"_pheno_",i,"_chr",j,"_model4_RES --thread-num 6"))
      }
    }
  }  else {
    for(i in 1:3){
      system(paste0("~/bin/gcta64 --mgrm ",output,"_GRM_mgrm.txt --pheno ",output,".phen --mpheno ",i," --reml --qcovar ",output,"_pca_rearing_qcovar.txt --out ",output,"_pheno_",i,"_RES --thread-num 6"))
      # And Per Chromosome
      # We run 4 models here as described above
      for(j in 1:23){
        # Model 1 = GRM excluding chr
        system(paste0("~/bin/gcta64 --grm ",output,"_GRM_NoFocalChr_chr",j," --pheno ",output,".phen --mpheno ",i," --reml --qcovar ",output,"_pca_rearing_qcovar.txt --out ",output,"_pheno_",i,"_chr",j,"_model1_RES --thread-num 6"))
        # Model 2 = GRM excluding chr and focal only
        system(paste0("~/bin/gcta64 --mgrm ",output,"_GRM_chr",j,"_mgrm_model2.txt --pheno ",output,".phen --mpheno ",i," --reml --qcovar ",output,"_pca_rearing_qcovar.txt --out ",output,"_pheno_",i,"_chr",j,"_model2_RES --thread-num 6"))
        # Model 3 = GRM all chr
        system(paste0("~/bin/gcta64 --grm ",output,"_GRM --pheno ",output,".phen --mpheno ",i," --reml --qcovar ",output,"_pca_rearing_qcovar.txt --out ",output,"_pheno_",i,"_chr",j,"_model3_RES --thread-num 6"))
        # Model 4 = GRM all chr + GRM focal chr
        system(paste0("~/bin/gcta64 --mgrm ",output,"_GRM_chr",j,"_mgrm_model4.txt --pheno ",output,".phen --mpheno ",i," --reml --qcovar ",output,"_pca_rearing_qcovar.txt --out ",output,"_pheno_",i,"_chr",j,"_model4_RES --thread-num 6"))
      }
    }
  }
})

# Phenos
pheno_vec <- c("Age at first brood (F)",
               "Size at first brood (F)",
               "First brood size (F)",
               "Interbrood period (F)",
               "First brood offspring weight (F)",
               "Age at maturity (M)",
               "Size at maturity (M)")

###############################
# Also read back in the per-chromosome estimates...
female_res <- data.frame(rbindlist(lapply(1:5,function(i){
  
  # Get res per chromosome
  chr_merge <- data.frame(rbindlist(lapply(1:23,function(chr){
    #print(chr)
    
    # List inputs
    inputs <- list.files("outputs/GCTA_chr_partitioning_with_litters/")
    inputs <- grep(paste0("female_qtl_crosses_heritability_pheno_",i,"_chr",chr,"_model"),inputs,value=T)
    inputs <- grep(".hsq",inputs,value = T)
    
    # Read in hsq file for all models
    model_res <- lapply(inputs,function(mod){
      # print(mod)
      read.table(paste0("outputs/GCTA_chr_partitioning_with_litters/",mod),fill=T,header=T)
    })
    
    # Make output
    out <- matrix(ncol=10,nrow=1)
    out[1,1] <- pheno_vec[i] # Phenotype
    out[1,2] <- paste0("chr",chr) # Chromosome
    
    if(length(inputs) == 4){
      out[1,3] <- model_res[[1]]$Variance[5] # logL model 1
      out[1,4] <- model_res[[2]]$Variance[9] # logL model 2
      out[1,5] <- model_res[[3]]$Variance[5] # logL model 3
      out[1,6] <- model_res[[4]]$Variance[9] # logL model 4
      
      # Calculate LRT for comps of 2v1 and 4v2
      out[1,7] <- -2*(as.numeric(out[1,3]) - as.numeric(out[1,4]))
      out[1,9] <- -2*(as.numeric(out[1,5]) - as.numeric(out[1,6]))
      
      # Calculate p values
      out[1,8] <- pchisq(as.numeric(out[1,7]),df=1,lower.tail=FALSE)
      out[1,10] <- pchisq(as.numeric(out[1,9]),df=1,lower.tail=FALSE)
      
      # Tidy up
      for(val in 7:10){
        out[1,val] <- round(as.numeric(out[1,val]),3)
      }
    } else {
      out[1,3:10] <- NA
    }
    
    return(data.frame(out))
  })))
  
  # Save 
  return(chr_merge)
})))

male_res <- data.frame(rbindlist(lapply(1:2,function(i){
  
  # Get res per chromosome
  chr_merge <- data.frame(rbindlist(lapply(1:23,function(chr){
    
    # List inputs
    inputs <- list.files("outputs/GCTA_chr_partitioning_with_litters/")
    inputs <- grep(paste0("male_qtl_crosses_heritability_pheno_",i,"_chr",chr,"_model"),inputs,value=T)
    inputs <- grep("female",inputs,invert=T,value=T)
    inputs <- grep(".hsq",inputs,value = T)
    
    # Read in hsq file for all models
    model_res <- lapply(inputs,function(mod){
      # print(mod)
      read.table(paste0("outputs/GCTA_chr_partitioning_with_litters/",mod),fill=T,header=T)
    })
    
    # Make output
    out <- matrix(ncol=10,nrow=1)
    out[1,1] <- pheno_vec[i+5] # Phenotype
    out[1,2] <- paste0("chr",chr) # Chromosome
    
    if(length(inputs) == 4){
      
      out[1,3] <- model_res[[1]]$Variance[5] # logL model 1
      out[1,4] <- model_res[[2]]$Variance[9] # logL model 2
      out[1,5] <- model_res[[3]]$Variance[5] # logL model 3
      out[1,6] <- model_res[[4]]$Variance[9] # logL model 4
      
      # Calculate LRT for comps of 2v1 and 4v2
      out[1,7] <- -2*(as.numeric(out[1,3]) - as.numeric(out[1,4]))
      out[1,9] <- -2*(as.numeric(out[1,5]) - as.numeric(out[1,6]))
      
      # Calculate p values
      out[1,8] <- pchisq(as.numeric(out[1,7]),df=1,lower.tail=FALSE)
      out[1,10] <- pchisq(as.numeric(out[1,9]),df=1,lower.tail=FALSE)
      
      # Tidy up
      for(val in 7:10){
        out[1,val] <- round(as.numeric(out[1,val]),3)
      }
    } else {
      out[1,3:10] <- NA
    }
    
    return(data.frame(out))
  })))
  
  # Save 
  return(chr_merge)
})))

# Per chr merge
greml_results_chr <- rbind(female_res,male_res)
colnames(greml_results_chr) <- c("Phenotype","Chromosome","Mod1 logL","Mod2 logL","Mod3 logL","Mod4 logL","LRT1","LRT1 p","LRT2","LRT2 p")

# FDR correct within phenotypes
greml_results_chr$`LRT1 fdr` <- NA
greml_results_chr$`LRT2 fdr` <- NA
for(pheno in pheno_vec){
  greml_results_chr[greml_results_chr$Phenotype == pheno,"LRT1 fdr"] <- p.adjust(greml_results_chr[greml_results_chr$Phenotype == pheno,"LRT1 p"],method = "fdr")
  greml_results_chr[greml_results_chr$Phenotype == pheno,"LRT2 fdr"] <- p.adjust(greml_results_chr[greml_results_chr$Phenotype == pheno,"LRT2 p"],method = "fdr")
}

greml_results_chr$`LRT1 fdr` <- round(greml_results_chr$`LRT1 fdr`,3)
greml_results_chr$`LRT2 fdr` <- round(greml_results_chr$`LRT2 fdr`,3)

# Where is strongest effects of chr?
head(greml_results_chr[order(greml_results_chr$`LRT1 p`),],10)
head(greml_results_chr[order(greml_results_chr$`LRT2 p`),],10)

# Output to table 
write.table(na.omit(greml_results_chr[order(greml_results_chr$`LRT1 p`),]),
            "tables/TableSX_chr_dropping_GRM_heritability.txt",
            row.names = F,quote = F,sep = "\t")

# How do these correlate with single chr estimates of h2c
single_chr_estimates <- read.csv("tables/TableSX_single_chr_estimates_h2c.csv",header=T)
single_chr_estimates$pheno_chr <- paste0(single_chr_estimates$Phenotype,"_",single_chr_estimates$chr)
greml_results_chr$pheno_chr <- paste0(greml_results_chr$Phenotype,"_",greml_results_chr$Chromosome)

# Filter single chr to get those that are also in na.omit greml_results_chr
greml_results_chr_na <- na.omit(greml_results_chr)
single_chr_estimates <- single_chr_estimates[single_chr_estimates$pheno_chr %in% greml_results_chr_na$pheno_chr,]

# SOrt each
greml_results_chr_na <- greml_results_chr_na[order(greml_results_chr_na$pheno_chr),]
single_chr_estimates <- single_chr_estimates[order(single_chr_estimates$pheno_chr),]

# Check sort
single_chr_estimates$Phenotype == greml_results_chr_na$Phenotype
single_chr_estimates$chr == greml_results_chr_na$Chromosome
single_chr_estimates$pheno_chr == greml_results_chr_na$pheno_chr

# Correlate
# pvals 
p_val_corr_res <- cor.test(as.numeric(greml_results_chr_na$`LRT1 p`),as.numeric(single_chr_estimates$p_val),method = "spearman")
to_plot <- data.frame(single_chr=as.numeric(single_chr_estimates$p_val),
                      model_chr=as.numeric(greml_results_chr_na$`LRT1 p`))

# Plot these
pdf("figs/FigureSX_correlation_between_h2c_approaches.pdf",width=6,height=6)
ggplot(to_plot,aes(x=log10(single_chr),y=log10(model_chr)))+
  geom_point()+
  theme_bw()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        title = element_text(size=20))+
  labs(x=expression(log[10](Single~Chr~GRM~p-value)),
       y=expression(log[10](LRT~Approach~p-value)))+
  ggtitle(paste0("Spearman's rho = ",round(p_val_corr_res$estimate,3)))+
  geom_abline(linetype="dotted")
dev.off()

