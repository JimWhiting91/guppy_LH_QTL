# Prepare all our inputs for qtl2 - Converting Lepmap3 genotypes to fully informative genotypes
lib<-c("data.table","vcfR","parallel")
lapply(lib,library,character.only=T)

####################################################################
# Genotypes, chr by chr (can combine later if needs be)
# Also trim the LG
trimming<-read.table("data/v10/trimming_guide.txt",header=T)
trimmed_map<-read.table("data/v10/trimmed_map.txt",header=T)

# Read in the VCF
vcf<-read.vcfR("data/qtl_crosses_FINAL_inds.sorted.filtered.vcf.gz")
chrBP<-paste0(vcf@fix[,1],"-",vcf@fix[,2])

####################################################################
# Find "fully informative" GP sites
all_gt<-extract.gt(vcf)
gt<-all_gt[,c(grep("Q",colnames(all_gt),value=T),
          grep("YH",colnames(all_gt),value=T))]

# Set up crosses
cross_info<-list(c("QF2YHM2","QF2_merged.sorted.rg","YHM2_rerun.merged.sorted.rg"),
                 c("QF4YHM8","QF4_merged.sorted.rg","YHM8_merged.sorted.rg"),
                 c("YHF5QM6","YHF5_merged.sorted.rg","QM6_merged.sorted.rg"),
                 c("YHF6QM7","YHF6_merged.sorted.rg","QM7_merged.sorted.rg"))

# Cross specific sites
cross_informative <- sapply(cross_info,function(x){
  # Filter
  gt_tmp<-gt[,x[2:3]]
  
  # Fully informative
  inform1<-which(gt_tmp[,1] == "0/0" & gt_tmp[,2] == "1/1")
  inform2<-which(gt_tmp[,1] == "1/1" & gt_tmp[,2] == "0/0")
  
  # Double-check informative sites...
  row_genos <- apply(gt_tmp[sort(c(inform1,inform2)),],1,paste,collapse="-")
  if(any(!(row_genos %in% c("0/0-1/1","1/1-0/0")))){
    message("Warning: Non-het sites detected in cross...")
  } else {
    message("All sites are informative :)")
  }
  
  return(chrBP[sort(c(inform1,inform2))])
})

# How many sites are informative for all families but not necessarily in pop-direction...
all_crosses<-table(unlist(cross_informative))
all_crosses<-names(all_crosses[all_crosses == 4])

# Common homozygote sites across all families, regardless of diretion?
all_informative<-table(unlist(cross_informative))
all_informative<-names(all_informative[all_informative == 4])

# Population specific sites
pop_informative<-na.omit(unlist(lapply(1:nrow(gt),function(x){

# Filter NAs
if(length(gt[x,][is.na(gt[x,])]) == 0)
{
# Test
pop_inform1 <- all(gt[x,1:4] == "0/0") & all(gt[x,5:8] == "1/1")
pop_inform2 <- all(gt[x,1:4] == "1/1") & all(gt[x,5:8] == "0/0")

# Return
if(pop_inform1 == TRUE | pop_inform2 == TRUE){
  return(chrBP[x])
} else {
return(NA)  
}
} else {
return(NA)  
}
})))

####################################################################
# This is the order of LG relative to CHR
chr_LG<-data.frame(LG=1:23,
                   chr=c(2,6,11,3,20,15,1,13,4,9,10,5,17,18,14,19,16,23,22,8,7,21,12))
chr_LG2<-chr_LG[order(chr_LG$chr),]

chr_genos<-mclapply(1:23,function(n){
  
  # First set the LG...
  x<-chr_LG2$LG[n]
  
  # Read in the data
  data_ = read.table(paste0("data/v10/qtl_cross_full_pedigree_v10_inform123_genotypes_LG",x,"_BEST.SA.txt"),header=F, sep="\t")
  mapped_info = read.table(paste0("data/v10/qtl_cross_full_pedigree_v10_inform123_order_LG",x,"_BEST.SA.mapped.clean.txt"),header=F, sep="\t")[,1:4]
  marker_info<-data_[,1:4]
  data_ = data_[,5:length(data_)]

  # Now we filter rows based on trimming and fetch the new trimmed positions
  trim_tmp<-trimming[trimming$LG == paste0("LG",x),]
  mapped_info<-mapped_info[trim_tmp$trim == "No",]
  pos<-trimmed_map[trimmed_map$LG == paste0("LG",x),"cM"]
  data_<-data_[trim_tmp$trim == "No",]

  # Tidy
  mapped_info$V2<-gsub("[*].*$","",mapped_info$V2)
  mapped_info$V1<-gsub("F0","F_0",mapped_info$V1)
  
  # Fetch the corresponding genotypes from the vcf
  to_keep<-unlist(sapply(1:nrow(mapped_info),function(x){
    if(length(which(vcf@fix[,1] == mapped_info$V1[x] &
            as.integer(vcf@fix[,2]) == mapped_info$V2[x])) > 1){
      return(NA)
    } else {
      return(which(vcf@fix[,1] == mapped_info$V1[x] &
                     as.integer(vcf@fix[,2]) == mapped_info$V2[x]))
    }
  }))
  
  # Amend pos
  pos<-pos[!(is.na(to_keep))]
  
  # Individuals based on lepmap filtered input...
  lep_inds<-as.character(read.table("data/lepmap_individuals_INPUT.txt")[,1])
  GP<-c(grep("Q",lep_inds,value=T),grep("YH",lep_inds,value=T))
  lep_inds<-lep_inds[!(lep_inds %in% GP)]
  individuals<-lep_inds
  
  ###############################################################
  # Sort by males and females
  ID_info<-read.table("data/qtl_cross_ID_metadata_FINAL.tsv",header=T)
  
  # Do males
  males<-na.omit(read.csv("data/qtl_male_phenotypes_clean.csv"))
  for(i in 1:nrow(males)){
    # print(i)
    if(length(as.character(ID_info[ID_info$TRUE_ID == males$ID[i],"VCF_ID"])) > 0)
    {
      males$vcf_id[i]<-as.character(ID_info[ID_info$TRUE_ID == males$ID[i],"VCF_ID"])
    }
    else {
      males$vcf_id[i]<-NA
    }
  }
  males<-males[males$vcf_id %in% individuals,]
  
  # Do females
  females<-na.omit(read.csv("data/qtl_female_phenotypes_clean.csv"))
  females<-na.omit(females)
  for(i in 1:nrow(females)){
    # print(i)
    if(length(as.character(ID_info[ID_info$TRUE_ID == females$ID[i],"VCF_ID"])) > 0)
    {
      females$vcf_id[i]<-as.character(ID_info[ID_info$TRUE_ID == females$ID[i],"VCF_ID"])
    }
    else {
      females$vcf_id[i]<-NA
    }
  }
  females<-females[females$vcf_id %in% individuals,]
  
  # Now filter genotypes as well
  colnames(data_)<-individuals
  data_<-data_[,individuals %in% c(males$vcf_id,females$vcf_id)]
  
  # And fetch
  male_inds<-individuals[individuals %in% males$vcf_id]
  female_inds<-individuals[individuals %in% females$vcf_id]
  individuals<-individuals[individuals %in% c(males$vcf_id,females$vcf_id)]
  
  # Set up family info
  family_info<-data.frame(ind=individuals)
  pedigree<-t(read.table("outputs/all_crosses_pedigree_dummyParents.txt",header=F))[,1:2]
  pedigree<-data.frame(pedigree[3:nrow(pedigree),])
  colnames(pedigree)<-c("family","ind")
  for(i in 1:nrow(family_info)){
    family_info$family[i]<-as.character(pedigree[as.character(pedigree$ind) == family_info$ind[i],"family"])
  }
  
  
  # Set logical for direction to translate genotypes
  QY<-grep("QM", family_info$family)
  YQ<-grep("YHM", family_info$family)
  genos<-as.matrix(data_[!(is.na(to_keep)),])
  
  # Rename geno rows
  rownames(genos)<-chrBP[na.omit(to_keep)]
  
  # Set them up differently
  genos[,QY][genos[,QY] == "1 1"]<-"LL"
  genos[,QY][genos[,QY] == "2 2"]<-"HH"
  genos[,QY][genos[,QY] == "1 2"]<-"HL"
  genos[,QY][genos[,QY] == "2 1"]<-"HL"

  genos[,YQ][genos[,YQ] == "2 2"]<-"LL"
  genos[,YQ][genos[,YQ] == "1 1"]<-"HH"
  genos[,YQ][genos[,YQ] == "2 1"]<-"HL"
  genos[,YQ][genos[,YQ] == "1 2"]<-"HL"
  
  ###############################################################
  
  male_genos<-t(genos[,colnames(genos) %in% male_inds])
  female_genos<-t(genos[,colnames(genos) %in% female_inds])
  
  return(list(male_genos,female_genos,pos,rownames(genos)))
},mc.cores=6)

# Extract males and females and cbind, also re-order...
male_genos<-data.frame(t(data.frame(rbindlist(lapply(1:23,function(x){return(data.frame(t(chr_genos[[x]][[1]])))})))))
female_genos<-data.frame(t(data.frame(rbindlist(lapply(1:23,function(x){return(data.frame(t(chr_genos[[x]][[2]])))})))))

marker_ids <- unlist(lapply(chr_genos,'[[',4))

# Sort out naming problenm
rownames(male_genos)<-gsub("X","",rownames(male_genos))
rownames(female_genos)<-gsub("X","",rownames(female_genos))

# Write the GM
genetic_map<-data.frame(rbindlist(lapply(1:23,function(x){
  real_chr<-paste0("chr",x)
  out<-data.frame(marker=chr_genos[[x]][[4]],
                  chr=real_chr,
                  pos=chr_genos[[x]][[3]])
  return(out)
})))

# How big is the map?
chr_lengths <- sapply(1:23,function(x){return(max(genetic_map[genetic_map$chr == paste0("chr",x),"pos"]))})
sum(chr_lengths)

# Write
write.csv(genetic_map,
          "data/lepmap_phased_genos_gmap.csv",row.names = F)

# Write to a file
colnames(male_genos)<-c(as.character(genetic_map$marker))
colnames(female_genos)<-c(as.character(genetic_map$marker))

write.csv(male_genos,
          "data/lepmap_phased_genos_MALE_allchr.csv")
write.csv(female_genos,
          "data/lepmap_phased_genos_FEMALE_allchr.csv")

####################################################################
# Also write only fully-informative markers
male_genos_informative<-male_genos[,colnames(male_genos) %in% pop_informative]
female_genos_informative<-female_genos[,colnames(female_genos) %in% pop_informative]
informative_map<-genetic_map[genetic_map$marker %in% pop_informative,]

# Write
write.csv(male_genos_informative,
          "data/lepmap_phased_genos_MALE_allchr_informative.csv")
write.csv(female_genos_informative,
          "data/lepmap_phased_genos_FEMALE_allchr_informative.csv")
write.csv(informative_map,
          "data/lepmap_phased_genos_gmap_informative.csv",row.names = F)

####################################################################
# Phenotypes
# Get our individuals
male_geno_ids<-as.character(read.csv("data/lepmap_phased_genos_MALE_allchr.csv")[,1])
female_geno_ids<-as.character(read.csv("data/lepmap_phased_genos_FEMALE_allchr.csv")[,1])
individuals<-c(male_geno_ids,female_geno_ids)

# Sort by males and females
ID_info<-read.table("data/qtl_cross_ID_metadata_FINAL.tsv",header=T)

# Do males
males<-na.omit(read.csv("data/qtl_male_phenotypes_clean.csv"))
for(i in 1:nrow(males)){
  # print(i)
  if(length(as.character(ID_info[ID_info$TRUE_ID == males$ID[i],"VCF_ID"])) > 0)
  {
    males$vcf_id[i]<-as.character(ID_info[ID_info$TRUE_ID == males$ID[i],"VCF_ID"])
  }
  else {
    males$vcf_id[i]<-NA
  }
}
males<-males[males$vcf_id %in% individuals,]

# Do females
females<-read.csv("data/qtl_female_phenotypes_clean_with_litters.csv")
#females<-na.omit(females)
for(i in 1:nrow(females)){
  #print(i)
  if(length(as.character(ID_info[ID_info$TRUE_ID == females$ID[i],"VCF_ID"])) > 0)
  {
    females$vcf_id[i]<-as.character(ID_info[ID_info$TRUE_ID == females$ID[i],"VCF_ID"])
  }
  else {
    females$vcf_id[i]<-NA
  }
}
females<-females[females$vcf_id %in% individuals,]

# Male phenos
male_phenos<-matrix(nrow=length(male_geno_ids),ncol=3)
for(i in 1:length(male_geno_ids)){
  male_phenos[i,1]<-as.numeric(na.omit(males[males$vcf_id == male_geno_ids[i],"Days"]))
  male_phenos[i,2]<-as.numeric(na.omit(males[males$vcf_id == male_geno_ids[i],c("Mature_length")]))
  male_phenos[i,3]<-as.numeric(na.omit(males[males$vcf_id == male_geno_ids[i],c("Mature_weight")]))
}

female_phenos<-matrix(nrow=length(female_geno_ids),ncol=6)
for(i in 1:length(female_geno_ids)){
  print(i)
  female_phenos[i,1]<-as.numeric(na.omit(females[females$vcf_id == female_geno_ids[i],c("age_first_brood")]))
  female_phenos[i,2]<-as.numeric(na.omit(females[females$vcf_id == female_geno_ids[i],c("size_at_first_brood")]))
  female_phenos[i,3]<-as.numeric(na.omit(females[females$vcf_id == female_geno_ids[i],c("first_brood_size")]))
  female_phenos[i,4]<-as.numeric(na.omit(females[females$vcf_id == female_geno_ids[i],c("interbrood")]))
  female_phenos[i,5]<-as.numeric(females[females$vcf_id == female_geno_ids[i],c("dry_offspring_weight_1")])
  female_phenos[i,6]<-as.numeric(females[females$vcf_id == female_geno_ids[i],c("dry_offspring_weight_2")])
}

# Tidy
rownames(male_phenos)<-male_geno_ids
rownames(female_phenos)<-female_geno_ids
colnames(male_phenos)<-c("Days","Mature_length","Mature_weight")
colnames(female_phenos)<-c("age_first_brood","size_at_first_brood","first_brood_size","interbrood",paste0("dry_offspring_weight_",1:2))

# Output in same order as genotypes...
write.csv(male_phenos,
          "data/male_phenotypes_qtl2_format.csv")
write.csv(female_phenos,
          "data/female_phenotypes_qtl2_format.csv")

##############################################
# Also make the covar and phenocovar files
male_covar<-data.frame(id=male_geno_ids,
                       sex="M")
female_covar<-data.frame(id=female_geno_ids,
                         sex="F")
# Get cross direction...
pedigree<-(read.table("outputs/all_crosses_pedigree_dummyParents.txt"))
pedigree<-data.frame(t(pedigree[,3:ncol(pedigree)]))

# Fetch cross
for(i in 1:nrow(male_covar)){
  male_covar$cross[i]<-as.character(pedigree[as.character(pedigree$X2) == as.character(male_covar$id[i]),"X1"])
  tmp<-substring(male_covar$cross[i], 1, 1)
  if(tmp == "Q"){
    male_covar$cross_direction[i]<-"(HxL)x(HxL)"
  } else {
    male_covar$cross_direction[i]<-"(LxH)x(LxH)"
  }
}
for(i in 1:nrow(female_covar)){
  female_covar$cross[i]<-as.character(pedigree[as.character(pedigree$X2) == as.character(female_covar$id[i]),"X1"])
  tmp<-substring(female_covar$cross[i], 1, 1)
  if(tmp == "Q"){
    female_covar$cross_direction[i]<-"(HxL)x(HxL)"
  } else {
    female_covar$cross_direction[i]<-"(LxH)x(LxH)"
  }
}

# Write them
write.csv(male_covar,
          "data/male_qtl2_covar.csv",row.names = F)
write.csv(female_covar,
          "data/female_qtl2_covar.csv",row.names = F)

################################################################
# We also want to produce cross-specific geno files
peds<-list.files("outputs/")
peds<-grep("dummyParents",peds,value = T)
peds<-grep("all",peds,value = T,invert = T)

names(cross_info)<-peds
names(cross_informative)<-peds

# Run over each one
lapply(peds,function(cross){
  
  # Fetch info
  inform_tmp<-cross_informative[[cross]]
  
  # Get the pedigrees
  pedigree<-read.table(paste0("outputs/",cross))
  pedigree2<-data.frame(t(pedigree[,3:ncol(pedigree)]))
  
  # Get cross and inds
  cross_id<-as.character(pedigree2$X1[1])
  cross_inds<-as.character(pedigree2$X2)
  
  # Filter the genos
  male_cross<-male_genos[rownames(male_genos) %in% cross_inds,colnames(male_genos) %in% inform_tmp]
  female_cross<-female_genos[rownames(female_genos) %in% cross_inds,colnames(female_genos) %in% inform_tmp]
  
  # Filter the phenos
  male_pheno_cross<-male_phenos[rownames(male_phenos) %in% cross_inds,]
  female_pheno_cross<-female_phenos[rownames(female_phenos) %in% cross_inds,]
  
  # Filter the phenocovar
  male_cross_covar<-male_covar[male_covar$id %in% cross_inds,]
  female_cross_covar<-female_covar[female_covar$id %in% cross_inds,]
  
  # Filter the genetic map
  cross_map<-genetic_map[genetic_map$marker %in% inform_tmp,]

  # Write the new genos + phenos
  write.csv(male_cross,
            paste0("data/lepmap_phased_genos_MALE_allchr_",cross_id,"_informative.csv"))
  write.csv(female_cross,
            paste0("data/lepmap_phased_genos_FEMALE_allchr_",cross_id,"_informative.csv"))
  write.csv(male_pheno_cross,
            paste0("data/male_phenotypes_qtl2_format_",cross_id,".csv"))
  write.csv(female_pheno_cross,
            paste0("data/female_phenotypes_qtl2_format_",cross_id,".csv"))
  write.csv(male_cross_covar,
            paste0("data/male_qtl2_covar_",cross_id,".csv"),
            row.names = F)
  write.csv(female_cross_covar,
            paste0("data/female_qtl2_covar_",cross_id,".csv"),
            row.names = F)
  write.csv(cross_map,
            paste0("data/lepmap_phased_genos_gmap_",cross_id,"_informative.csv"),row.names = F)
})

# Finally write the phenocovar csv
male_phenocovar <- data.frame(phenos=colnames(male_phenos),
                              type=colnames(male_phenos))
female_phenocovar <- data.frame(phenos=colnames(female_phenos),
                              type=colnames(female_phenos))
write.csv(male_phenocovar,
          "data/male_qtl2_phenocovar.csv",
          quote = F,row.names = F)
write.csv(female_phenocovar,
          "data/female_qtl2_phenocovar.csv",
          quote = F,row.names = F)
