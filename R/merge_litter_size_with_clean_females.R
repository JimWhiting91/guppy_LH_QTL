##### Short script to merge litter size with cleaned phenotypes #####
lib<-c("data.table")
lapply(lib,library,character.only=T)

# Get the data
litter <- read.csv("raw_data/litter_data.csv",header=T)
colnames(litter) <- c("Line","ID","Date1","Size1","Wet1","Dry1",
                      "Date2","Size2","Wet2","Dry2",
                      "Date3","Size3","Wet3","Dry3")
litter <- litter[,c("Line","ID","Date1","Size1","Wet1","Dry1",
                    "Date2","Size2","Wet2","Dry2",
                    "Date3","Size3","Wet3","Dry3")]
litter<-data.frame(litter)

# Get mean outliers
plot(litter$Dry1~litter$Size1)
plot(litter$Dry2~litter$Size2)

# Get mean dry weight
inds <- litter$ID
litter_mat <- matrix(ncol=4,nrow=length(inds))
litter_mat[,1] <- inds
for(i in 1:3){
  
  # Add together any that need adding...
  dry_mean <- sapply(inds,function(id){
    tmp <- litter[litter$ID == id,]
    tmp <- tmp[complete.cases(tmp[ ,c(paste0("Dry",i),paste0("Size",i)),]),]
    
    # Remove bad ones
    if(i == 1){
      tmp[tmp$Dry1 > 0.01 & tmp$Size1[i] < 10,"Dry1"] <- NA
    } else if (i == 2){
      tmp[tmp$Dry2 > 0.04,"Dry2"] <- NA
    }
    
    sum(as.numeric(tmp[,paste0("Dry",i)]))/sum(as.numeric(tmp[,paste0("Size",i)]))
  })

  litter_mat[,i+1] <- dry_mean  
}

colnames(litter_mat) <- c("inds",paste0("Dry_Mean",1:3))
litter_mat <- data.frame(litter_mat)

# Check final distributions
hist(log(litter_mat$Dry_Mean1))
hist(log(litter_mat$Dry_Mean2))

# Visualise together
plot(log(litter_mat$Dry_Mean1)~log(litter_mat$Dry_Mean2))
cor.test(log(litter_mat$Dry_Mean1),log(litter_mat$Dry_Mean2))

# Get the clean data
clean_phenos <- read.csv("data/qtl_female_phenotypes_clean.csv")

# Visualise interbrood vs brood size
plot(clean_phenos$interbrood,clean_phenos$first_brood_size)

# Merge them
for(i in 1:nrow(clean_phenos)){
  print(i)
  clean_phenos$dry_offspring_weight_1[i] <- litter_mat[litter_mat$inds == as.character(clean_phenos$ID[i]),"Dry_Mean1"]
  clean_phenos$dry_offspring_weight_2[i] <- litter_mat[litter_mat$inds == as.character(clean_phenos$ID[i]),"Dry_Mean2"]
}

# Visualise with interbrood, are very high dry weight linked with interbrood?
plot(clean_phenos$age_first_brood~clean_phenos$dry_offspring_weight_1)
cor.test(clean_phenos$age_first_brood,clean_phenos$dry_offspring_weight_1)

plot(clean_phenos$size_at_first_brood~clean_phenos$dry_offspring_weight_1)
cor.test(clean_phenos$size_at_first_brood,clean_phenos$dry_offspring_weight_1)

plot(clean_phenos$first_brood_size~clean_phenos$dry_offspring_weight_1)
cor.test(clean_phenos$first_brood_size,clean_phenos$dry_offspring_weight_1)

plot(clean_phenos$interbrood~clean_phenos$dry_offspring_weight_1)
cor.test(clean_phenos$interbrood,clean_phenos$dry_offspring_weight_1)

# Check PCA
female_pca <- prcomp(na.omit(clean_phenos[,3:ncol(clean_phenos)]),scale. = T,center = T)
female_pca$rotation

# Write this
write.csv(clean_phenos,
          "data/qtl_female_phenotypes_clean_with_litters.csv",quote = F,row.names = F)
