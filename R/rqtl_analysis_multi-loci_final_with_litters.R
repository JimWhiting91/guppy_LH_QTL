# rQTL Analysis for one and two-loci models -------------------------------
lib<-c("qtl","ggplot2")
lapply(lib,library,character.only=T)

# Which sex?
male_cross<-"outputs/rqtl_format_male_informative_with_phenos.csv"
female_cross<-"outputs/rqtl_format_female_informative_with_phenos.csv"

# Fetch rearing covariates and amend IDs
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

# Count the crossovers...
# Read in the cross
females<-read.cross(format="csv",
                    file=female_cross)
# Read in the cross
males<-read.cross(format="csv",
                    file=male_cross)

crossovers <- data.frame(sex=c(rep("female",nrow(females$pheno)),
                               rep("male",nrow(males$pheno))),
                         xo=c(countXO(females),countXO(males)))
ggplot(crossovers,aes(x=xo))+
  geom_density()+
  facet_wrap(~sex,ncol=1)

#  Female QTL Analysis ----------------------------------------------------
summary(females)

# Jitter the map a bit
females<-jittermap(females)
females2 <- qtl2::read_cross2('data/qtl_cross_females_informative.yaml')

# Make family and temperature covars...
# Set up family and rearing conditions as  covariate...
family <- factor(females2$covar$cross)

# use model.matrix to create a numeric matrix of 0/1 covariates
family_matrix <- model.matrix( ~ family)

# omit the first column; drop=FALSE ensures it stays as a matrix
family_matrix <- family_matrix[,-1, drop=FALSE]
rownames(family_matrix) <- rownames(females2$covar)

# and rearing...
rear_mat <- matrix(nrow = nrow(family_matrix),ncol=2)
rownames(rear_mat) <- rownames(family_matrix)
colnames(rear_mat) <- c("DOB","mean_temp")
for(i in 1:nrow(rear_mat)){
  if(length(na.omit(female_rearing[female_rearing$vcf_id == rownames(rear_mat)[i],"DOB"])) != 0){
    rear_mat[i,1] <- na.omit(female_rearing[female_rearing$vcf_id == rownames(rear_mat)[i],"DOB"])
  }
  if(length(na.omit(female_rearing[female_rearing$vcf_id == rownames(rear_mat)[i],"mean_temp"])) != 0){
    rear_mat[i,2] <- na.omit(female_rearing[female_rearing$vcf_id == rownames(rear_mat)[i],"mean_temp"])
  }}

# Merge
family_rear_covar <- cbind(family_matrix,rear_mat)

# Transform and edit female phenos
female_phenos <- females$pheno

# Filter out bad interbroods...
female_phenos[female_phenos$interbrood < 15 | 
                female_phenos$interbrood > 40, c("interbrood",
                                   "first_brood_size",
                                   "dry_offspring_weight_1")]<-NA

# Get residual first brood size
size_mod <- lm(first_brood_size~size_at_first_brood,data=female_phenos)
for(i in 1:nrow(phenos)){
  female_phenos$predicted_first_brood_size[i] <- predict(size_mod)[rownames(female_phenos)[i]]
}

female_phenos$first_brood_size_resid <- female_phenos$first_brood_size - as.numeric(female_phenos$predicted_first_brood_size)
#hist(phenos$first_brood_size_resid)

# Log-transform all values, test, check hists
to_transform <- c("age_first_brood","size_at_first_brood","interbrood","dry_offspring_weight_1")
for(i in to_transform){
  female_phenos[,i] <- log(female_phenos[,i])
  #print(shapiro.test(phenos[,i]))
}

# Set up final output
female_phenos <- female_phenos[,c("age_first_brood","size_at_first_brood","first_brood_size_resid","interbrood","dry_offspring_weight_1")]
female_phenos <- cbind(females$pheno[,1],female_phenos)
colnames(female_phenos)[1] <- colnames(females$pheno)[1]
females$pheno <- female_phenos

# Simulate genotype probabilities
females <- calc.genoprob(females, step=2, error.prob=0.01)
females <- sim.geno(females, step=2, n.draws = 16, error.prob=0.01)

# Step analysis first to test for any potential multi-QTLs.
# Additive
if(!(file.exists("outputs/female_stepqtl_additive_models.rds"))){
female_step_out <- mclapply(2:6,function(i){
tmp_step <- stepwiseqtl(females, additive.only=TRUE, max.qtl=6,pheno.col = i,method = "imp",covar = family_rear_covar)
return(tmp_step)
},mc.cores=5)
saveRDS(female_step_out,"outputs/female_stepqtl_additive_models.rds")
}
female_step_out_add <- readRDS("outputs/female_stepqtl_additive_models.rds")

# Check each model
female_step_out_add[[1]]
female_step_out_add[[2]]
female_step_out_add[[3]]
female_step_out_add[[4]]
female_step_out_add[[5]]

##################################################################################### 
# Run two-loci model across all phenos...
if(!(file.exists("outputs/female_scantwo_models.rds"))){
female_pheno_out <- mclapply(2:6,function(i){
  female_pheno_model <- scantwo(females, method="imp",pheno.col = i,addcovar = family_rear_covar)
  return(female_pheno_model)
},mc.cores=5)
saveRDS(female_pheno_out,"outputs/female_scantwo_models.rds")
}
female_pheno_out <- readRDS("outputs/female_scantwo_models.rds")


if(!(file.exists("outputs/female_scantwo_models_perms.rds"))){
  female_pheno_perms <- mclapply(2:6,function(i){
    female_pheno_model <- scantwo(females, method="imp",pheno.col = i,n.perm = 100,addcovar = family_rear_covar)
    return(female_pheno_model)
  },mc.cores=4)
  saveRDS(female_pheno_perms,"outputs/female_scantwo_models_perms.rds")
}
female_pheno_perms <- readRDS("outputs/female_scantwo_models_perms.rds")

# Summarise
for(i in 1:5){
summary(female_pheno_out[[i]], perms=female_pheno_perms[[i]], pvalues=TRUE,
        alphas=c(0.05, 0.05, 0, 0.05, 0.05))
}
# Plot
for(i in 1:5){
plot(female_pheno_out[[5]])
}


# Male QTL Analysis -------------------------------------------------------

# Read in the cross
males<-read.cross(format="csv",
                    file=male_cross)
summary(males)

# Jitter the map a bit
males<-jittermap(males)
males2 <- qtl2::read_cross2('data/qtl_cross_males_informative.yaml')

# Make family and temperature covars...
# Set up family and rearing conditions as  covariate...
family <- factor(males2$covar$cross)

# use model.matrix to create a numeric matrix of 0/1 covariates
family_matrix <- model.matrix( ~ family)

# omit the first column; drop=FALSE ensures it stays as a matrix
family_matrix <- family_matrix[,-1, drop=FALSE]
rownames(family_matrix) <- rownames(males2$covar)

# and rearing...
rear_mat <- matrix(nrow = nrow(family_matrix),ncol=2)
rownames(rear_mat) <- rownames(family_matrix)
colnames(rear_mat) <- c("DOB","mean_temp")
for(i in 1:nrow(rear_mat)){
  if(length(na.omit(male_rearing[male_rearing$vcf_id == rownames(rear_mat)[i],"DOB"])) != 0){
    rear_mat[i,1] <- na.omit(male_rearing[male_rearing$vcf_id == rownames(rear_mat)[i],"DOB"])
  }
  if(length(na.omit(male_rearing[male_rearing$vcf_id == rownames(rear_mat)[i],"mean_temp"])) != 0){
    rear_mat[i,2] <- na.omit(male_rearing[male_rearing$vcf_id == rownames(rear_mat)[i],"mean_temp"])
  }}

# Merge
family_rear_covar <- cbind(family_matrix,rear_mat)

# Simulate genotype probabilities
males <- calc.genoprob(males, step=2, error.prob=0.01)
males <- sim.geno(males, step=2, n.draws = 16, error.prob=0.01)

# Step analysis first to test for any potential multi-QTLs.
# Additive
if(!(file.exists("outputs/male_stepqtl_additive_models.rds"))){
  male_step_out <- mclapply(2:3,function(i){
    tmp_step <- stepwiseqtl(males, additive.only=TRUE, max.qtl=6,pheno.col = i,method = "imp",covar = family_rear_covar)
    return(tmp_step)
  },mc.cores=2)
  saveRDS(male_step_out,"outputs/male_stepqtl_additive_models.rds")
}
male_step_out_add <- readRDS("outputs/male_stepqtl_additive_models.rds")

# Observe scantwo model for male age
male_age_scantwo <- scantwo(males, method="imp",pheno.col = 1,addcovar = family_rear_covar)


##################################################################################### 
# Run two-loci model across all phenos...
if(!(file.exists("outputs/male_scantwo_models.rds"))){
  male_pheno_out <- mclapply(2:5,function(i){
    male_pheno_model <- scantwo(males, method="imp",pheno.col = i,addcovar = family_rear_covar)
    return(male_pheno_model)
  },mc.cores=4)
  saveRDS(male_pheno_out,"outputs/male_scantwo_models.rds")
}
