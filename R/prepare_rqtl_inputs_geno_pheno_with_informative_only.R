# Prepare all our inputs for r/qtl but using only informative sites
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
cross_informative<-sapply(cross_info,function(x){
  # Filter
  gt_tmp<-gt[,x[2:3]]
  
  # Fully informative
  inform1<-which(gt_tmp[,1] == "0/0" & gt_tmp[,2] == "1/1")
  inform2<-which(gt_tmp[,1] == "1/1" & gt_tmp[,2] == "0/0")
  
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
    pop_inform1<-all(gt[x,1:4] == "0/0") & all(gt[x,5:8] == "1/1")
    pop_inform2<-all(gt[x,1:4] == "1/1") & all(gt[x,5:8] == "0/0")
    
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
  genos[,QY][genos[,QY] == "1 1"]<-"A"
  genos[,QY][genos[,QY] == "2 2"]<-"B"
  genos[,QY][genos[,QY] == "1 2"]<-"H"
  genos[,QY][genos[,QY] == "2 1"]<-"H"
  
  genos[,YQ][genos[,YQ] == "2 2"]<-"A"
  genos[,YQ][genos[,YQ] == "1 1"]<-"B"
  genos[,YQ][genos[,YQ] == "2 1"]<-"H"
  genos[,YQ][genos[,YQ] == "1 2"]<-"H"
  
  ###############################################################
  
  male_genos<-t(genos[,colnames(genos) %in% male_inds])
  female_genos<-t(genos[,colnames(genos) %in% female_inds])
  
  return(list(male_genos,female_genos,pos,rownames(genos)))
},mc.cores=6)

# Extract males and females and cbind, also re-order...
male_genos<-data.frame(t(data.frame(rbindlist(lapply(1:23,function(x){return(data.frame(t(chr_genos[[x]][[1]])))})))))
female_genos<-data.frame(t(data.frame(rbindlist(lapply(1:23,function(x){return(data.frame(t(chr_genos[[x]][[2]])))})))))

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

# Get each into r/qtl format
genetic_map$chr<-gsub("chr","",genetic_map$chr)
map2<-data.frame(t(genetic_map))

# Merge map with genos
male_merge<-rbind(map2,male_genos)
female_merge<-rbind(map2,female_genos)

####################################################################
# Also write only fully-informative markers
male_informative<-male_merge[,genetic_map$marker %in% pop_informative]
female_informative<-female_merge[,genetic_map$marker %in% pop_informative]

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
female_qtl2 <- read_cross2("data/qtl_cross_females_informative.yaml")
female_qtl2_phenos <- female_qtl2$pheno

# Male phenos
male_phenos<-matrix(nrow=length(male_geno_ids),ncol=3)
for(i in 1:length(male_geno_ids)){
  male_phenos[i,1]<-as.numeric(na.omit(males[males$vcf_id == male_geno_ids[i],"Days"]))
  male_phenos[i,2]<-as.numeric(na.omit(males[males$vcf_id == male_geno_ids[i],c("Mature_length")]))
  male_phenos[i,3]<-as.numeric(na.omit(males[males$vcf_id == male_geno_ids[i],c("Mature_weight")]))
}

female_phenos<-matrix(nrow=length(female_geno_ids),ncol=5)
for(i in 1:length(female_geno_ids)){
 # print(i)
  female_phenos[i,1]<-as.numeric(na.omit(female_qtl2_phenos[female_geno_ids[i],c("age_first_brood")]))
  female_phenos[i,2]<-as.numeric(na.omit(female_qtl2_phenos[female_geno_ids[i],c("size_at_first_brood")]))
  female_phenos[i,3]<-as.numeric(na.omit(female_qtl2_phenos[female_geno_ids[i],c("first_brood_size")]))
  female_phenos[i,4]<-as.numeric(na.omit(female_qtl2_phenos[female_geno_ids[i],c("interbrood")]))
  female_phenos[i,5]<-as.numeric(female_qtl2_phenos[female_geno_ids[i],c("dry_offspring_weight_1")])
}

# Tidy
rownames(male_phenos)<-male_geno_ids
rownames(female_phenos)<-female_geno_ids
colnames(male_phenos)<-c("Days","Mature_length","Mature_weight")
colnames(female_phenos)<-c("age_first_brood","size_at_first_brood","first_brood_size","interbrood","dry_offspring_weight_1")

# Merge phenotypes with informative genotypes
male_pheno_empty<-data.frame(Days=c("Days"," "," "),
                             Mature_length=c("Mature_length"," "," "),
                             Mature_weight=c("Mature_weight"," "," "))
female_pheno_empty<-data.frame(age_first_brood=c("age_first_brood"," "," "),
                               size_at_first_brood=c("size_at_first_brood"," "," "),
                               first_brood_size=c("first_brood_size"," "," "),
                               interbrood=c("interbrood"," "," "),
                               dry_offspring_weight_1=c("dry_offspring_weight_1"," "," "))

# And combine
male_informative_pheno<-cbind(rbind(male_pheno_empty,male_phenos),male_informative)
female_informative_pheno<-cbind(rbind(female_pheno_empty,female_phenos),female_informative)

# Tidy
rownames(male_informative_pheno)[1:3]<-c(" ","  ","   ")
rownames(female_informative_pheno)[1:3]<-c(" ","  ","   ")

# Write these files
write.table(male_informative_pheno,
            "outputs/rqtl_format_male_informative_with_phenos.csv",
            col.names = FALSE,quote = F,row.names = TRUE,sep=",")

write.table(female_informative_pheno,
            "outputs/rqtl_format_female_informative_with_phenos.csv",
            col.names = FALSE,quote = F,row.names = TRUE,sep=",")
