###################################################################################
# Make the pedigree file for lepmap3 analysis

# We need these
lib<-c("vcfR","data.table","parallel","adegenet")
lapply(lib,library,character.only=T)

# Fetch data
metadata<-read.table("data/qtl_cross_ID_metadata_FINAL.tsv",header=T)
metadata<-metadata[metadata$TO_KEEP != "LOSE",]
inds<-metadata$VCF_ID
GPs<-metadata[metadata$TO_KEEP == "GP",]

# Get males and females
male_phenos<-read.csv("data/qtl_male_phenotypes_clean.csv")
males<-na.omit(male_phenos$ID)
males<-as.character(metadata[metadata$TRUE_ID %in% males,"VCF_ID"])
females<-as.character(inds[!(inds %in% males)])

# Make pedigree 
pedigree<-data.frame(family=metadata$CROSS,
                     ind=inds)

pedigree$father<-0
pedigree$mother<-0
pedigree$sex<-0
pedigree$pheno<-0


# Get the father and mother and sex
for(i in 1:length(inds)){

    # Handle GPs
  if(inds[i] %in% GPs$VCF_ID){
    
    # Just get sex
    if(length(grep("F",as.character(inds[i]))) == 1){
      pedigree$sex[i]<-2
    } else {
      pedigree$sex[i]<-1
    }
  } else {
  
  # And the F2s
  GP_tmp<-metadata[metadata$CROSS ==  metadata[metadata$VCF_ID == inds[i],"CROSS"] &
                     metadata$TO_KEEP == "GP",]
  
  pedigree$father[i]<-grep("M",GP_tmp$VCF_ID,value = T)
  pedigree$mother[i]<-grep("F",GP_tmp$VCF_ID,value = T)
  
  if(inds[i] %in% males){
    pedigree$sex[i] <- 1
  } else {
    pedigree$sex[i] <- 2
  }

  # End
  }
}

# Add in some random dummy parents as well for each cross...
dummy_parents<-data.frame(rbindlist(lapply(unique(pedigree$family),function(fam){
  
  # Extract 
  tmp<-pedigree[pedigree$family==fam,]
  
  # Get the male
  male<-tmp$father[1]
  female<-tmp$mother[1]
  
  # Parents
  parents<-data.frame(family=fam,
                      ind=c(paste0(fam,c("_P1","_P2"))),
                      father=male,
                      mother=female,
                      sex=c(1,2),
                      pheno=0)

  # And reset the parents of all F2s
  tmp[tmp$father==male,"father"]<-paste0(fam,"_P1")
  tmp[tmp$mother==female,"mother"]<-paste0(fam,"_P2")
  
  # Rbind
  out<-rbind(tmp,parents)
  
  return(out)
})))

###################################################################################
# Now we want to remove individuals on the basis of PCA space
vcf<-read.vcfR("data/qtl_crosses_FINAL_inds.sorted.vcf.gz")

# Convert to genind
pca_gen<-vcfR2genind(vcf)
x.gen<-tab(pca_gen,freq=TRUE,NA.method="mean")

# Do PCA
pca <- dudi.pca(x.gen,center=T,scale = F,scannf = FALSE, nf=10)
scores<-data.frame(pca$li)
scores$ind<-rownames(pca_gen$tab)
scores$family=NA
for(i in 1:nrow(scores)){
  scores$family[i]<-as.character(pedigree[as.character(pedigree$ind)==scores$ind[i],"family"])
}

# Save PCA scores
write.csv(scores,
          "outputs/qtl_cross_individuals_preFilter_PCA_scores.csv",
          quote = F,row.names = F)

# All individuals
starting_individuals<-rownames(scores)

# For each family, identify "mislabelled"
family_vec<-unique(scores$family)
YHF5QM6_scores<-scores[scores$family == "YHF5QM6",]
YHF5QM6_scores[order(YHF5QM6_scores$Axis1),]
YHF5QM6_scores[order(YHF5QM6_scores$Axis2),]

YHF6QM7_scores<-scores[scores$family == "YHF6QM7",]
YHF6QM7_scores[order(-YHF6QM7_scores$Axis1),]

QF2YHM2_scores<-scores[scores$family == "QF2YHM2",]
head(QF2YHM2_scores[order(-QF2YHM2_scores$Axis2),1:2])

QF4YHM8_scores<-scores[scores$family == "QF4YHM8",]
head(QF4YHM8_scores[order(-QF4YHM8_scores$Axis2),1:2])
head(QF4YHM8_scores[order(QF4YHM8_scores$Axis1),1:2])

# Mislabelled individuals are, these need removing:
mislabelled <- c("merged_240.sorted.rg","merged_054.sorted.rg","merged_597.sorted.rg",
                 "09_557.sorted.rg","10_846.sorted.rg","merged_046.sorted.rg",
                 "1_063.sorted.rg",  "10_835.sorted.rg", "10_839.sorted.rg","11_979.sorted.rg", "11_973.sorted.rg", "11_983.sorted.rg", "11_980.sorted.rg","11_974.sorted.rg", "11_975.sorted.rg",
                 "1_101.sorted.rg","10_887.sorted.rg","merged_027.sorted.rg","3_314.sorted.rg","8_235.sorted.rg","3_272.sorted.rg")

# Set cross ID
s.class(pca$li, fac=as.factor(scores$family), col=funky(15))

# Now mark individuals that deviate too far from their cluster
cross_lists<-list(c("QF2YHM2","QF2_merged.sorted.rg","YHM2_rerun.merged.sorted.rg"),
                  c("QF4YHM8","QF4_merged.sorted.rg","YHM8_merged.sorted.rg"),
                  c("YHF5QM6","YHF5_merged.sorted.rg","QM6_merged.sorted.rg"),
                  c("YHF6QM7","YHF6_merged.sorted.rg","QM7_merged.sorted.rg"))

# Function for calculating outliers based on distance
calc_distance_outliers<-function(x){
 # print(x)
  
  cross<-x[[1]]
  gp1<-x[[2]]
  gp2<-x[[3]]
  
  # Get means
  tmp<-scores[scores$family == cross ,paste0("Axis",1:2)]
  tmp2<-rbind(tmp,colMeans(tmp))
  rownames(tmp2)[nrow(tmp2)]<-"AVG"
  
  dists<-as.matrix(dist(tmp2,method="euclidean"))
  colnames(dists)<-rownames(tmp2)
  rownames(dists)<-colnames(dists)

  # Trim
  dists2<-data.frame(dists[!(colnames(dists) %in% c(gp1,gp2)),"AVG"])
  dists2$cross<-cross
  colnames(dists2)<-c("AVG","Cross")
  
  cutoff<-median(dists2$AVG)+2*(sd(dists2$AVG))

  dists2$Outlier<-"No"
  dists2[dists2$AVG > cutoff,"Outlier"]<-"Yes"

  dists2$ID<-as.character(rownames(dists)[!(rownames(dists) %in% c(gp1,gp2))])
  
  # Add scores
  for(i in 1:nrow(dists2)){
    dists2$PC1[i]<-tmp2[rownames(tmp2)==dists2$ID[i],"Axis1"]
    dists2$PC2[i]<-tmp2[rownames(tmp2)==dists2$ID[i],"Axis2"]
  }
  # Return
  return(dists2)
}

#######################################################
# Run outliers Step 1
outlier_dist<-data.frame(rbindlist(lapply(cross_lists,calc_distance_outliers)))

# Check outliers
outlier_dist[outlier_dist$Outlier == "Yes",]

# Plot to evaluate
plot_dd<-scores[scores$ind %in% outlier_dist$ID,]

# Match order with inds
outlier_dist<-outlier_dist[match(plot_dd$ind, outlier_dist$ID),]
s.class(plot_dd[outlier_dist$Outlier == "No",], fac=factor(plot_dd$family[outlier_dist$Outlier == "No"],levels=unique(plot_dd$family)), col=funky(15))

# And remove these individuals and repeat
scores<-scores[scores$ind %in% plot_dd[outlier_dist$Outlier == "No","ind"],]

# Repeat
# Run outliers Step 2
outlier_dist<-data.frame(rbindlist(lapply(cross_lists,calc_distance_outliers)))

# Check outliers
outlier_dist[outlier_dist$Outlier == "Yes",]

# Plot to evaluate
plot_dd<-scores[scores$ind %in% outlier_dist$ID,]

# Match order with inds
outlier_dist<-outlier_dist[match(plot_dd$ind, outlier_dist$ID),]
s.class(plot_dd[outlier_dist$Outlier == "No",], fac=factor(plot_dd$family[outlier_dist$Outlier == "No"],levels=unique(plot_dd$family)), col=funky(15))

# Remove these individuals as well
to_keep<-c(plot_dd[outlier_dist$Outlier == "No","ind"],
           as.character(GPs$VCF_ID),
           grep("_P",dummy_parents$ind,value=T))

# And mark removed...
removed<-starting_individuals[!(starting_individuals %in% to_keep)]
hybrids<-removed[!(removed %in% mislabelled)]

####################################################################################
# And filter the pedigree for individuals that do not cluster with their family
pedigree<-pedigree[pedigree$ind %in% to_keep,]
dummy_parents<-dummy_parents[dummy_parents$ind %in% to_keep,]

# Report pedigree information
# How many per cross?
family_sex<-paste0(pedigree$family,"_",pedigree$sex)
table(family_sex)-1
# How many F2s?
sum(table(family_sex)-1)

# Subset by cross
lapply(unique(pedigree$family),function(fam){
  
  # Extract and transpose
  tmp<-t(pedigree[pedigree$family==fam,])
  tmp2<-t(dummy_parents[dummy_parents$family==fam,])
  
  # Add the extra columns
  tmp_out<-cbind(rep("CHR",6),
                 rep("POS",6),
                 tmp)
  tmp2_out<-cbind(rep("CHR",6),
                 rep("POS",6),
                 tmp2)
  
  # Write
  write.table(tmp_out,
              paste0("outputs/",fam,"_pedigree.txt"),
              sep="\t",row.names = F,quote = F,col.names = F)
  write.table(tmp2_out,
              paste0("outputs/",fam,"_pedigree_dummyParents.txt"),
              sep="\t",row.names = F,quote = F,col.names = F)
})

# Write the whole thing
# Add the extra columns
full_out<-cbind(rep("CHR",6),
               rep("POS",6),
               t(pedigree))

full_out_dummy<-cbind(rep("CHR",6),
                rep("POS",6),
                t(dummy_parents))

write.table(full_out,
            "outputs/all_crosses_pedigree.txt",
            sep="\t",row.names = F,quote = F,col.names = F)

write.table(full_out_dummy,
            "outputs/all_crosses_pedigree_dummyParents.txt",
            sep="\t",row.names = F,quote = F,col.names = F)


###################################################################################################################

