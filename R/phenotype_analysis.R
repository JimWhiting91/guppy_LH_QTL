######## Analysis of QTL Phenotypes
lib <- c("DHARMa","data.table","ggplot2","qtl2","effects","see","reshape2","multcompView","wesanderson","Hmisc","cowplot","rsq")
lapply(lib,library,character.only=T)

#########################################################
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


#### Males #####
males <- read_cross2("data/qtl_cross_males_informative.yaml")
pheno_dd <- data.frame(inds=rownames(males$covar),
                       family=males$covar$cross)
pheno_dd <- cbind(pheno_dd,males$pheno[,1:2])

# Add rearing
pheno_dd$temp <- NA
pheno_dd$DOB <- NA
for(i in 1:nrow(pheno_dd)){
  if(length(na.omit(male_rearing[male_rearing$vcf_id == pheno_dd$inds[i],"mean_temp"]))==1){
    pheno_dd$temp[i] <- na.omit(male_rearing[male_rearing$vcf_id == pheno_dd$inds[i],"mean_temp"])
  } else {
    pheno_dd$temp[i] <- NA 
  }
  if(length(na.omit(male_rearing[male_rearing$vcf_id == pheno_dd$inds[i],"DOB"]))==1){
    pheno_dd$DOB[i] <- na.omit(male_rearing[male_rearing$vcf_id == pheno_dd$inds[i],"DOB"])
  } else {
    pheno_dd$DOB[i] <- NA 
  }
}

# Correlations between phenotypes and covar
male_pheno_covar <- rcorr(as.matrix(pheno_dd[,c("Days","Mature_length","temp","DOB")]),type = "spearman")
male_pheno_covar$r
male_pheno_covar$P
male_covar_cor_dd <- melt(male_pheno_covar$r)
male_covar_cor_dd$p <- melt(male_pheno_covar$P)$value
male_covar_cor_dd <- na.omit(male_covar_cor_dd)
# Filter for comps of interest
male_covar_cor_dd<- male_covar_cor_dd[male_covar_cor_dd$Var1 %in% c("DOB","temp") & male_covar_cor_dd$Var1 %in% c("Days","Mature_length") |
                                        male_covar_cor_dd$Var2 %in% c("DOB","temp") & male_covar_cor_dd$Var1 %in% c("Days","Mature_length"),]
male_covar_cor_dd$fdr <- p.adjust(male_covar_cor_dd$p,method = "fdr")
male_covar_cor_dd$signif <- NA
male_covar_cor_dd[male_covar_cor_dd$fdr < 0.05,"signif"] <- "*"
male_covar_cor_dd[male_covar_cor_dd$fdr < 0.01,"signif"] <- "**"
male_covar_cor_dd[male_covar_cor_dd$fdr < 0.001,"signif"] <- "***"
for(col in c("value","p","fdr")){
  male_covar_cor_dd[,col] <- round(male_covar_cor_dd[,col],3)
}
male_covar_cor_dd

####################################################################################
# Compare Days between families...
hist(log(pheno_dd$Days))
shapiro.test(pheno_dd$Days)
shapiro.test(log(pheno_dd$Days))
days_glm <- glm(Days~family+temp+DOB,data=pheno_dd,family=Gamma(link="log"))
days_glm_step <- step(days_glm)

# Check model...
shapiro.test(days_glm_step$residuals)
plot(days_glm_step)
check_gamma_model <- simulateResiduals(fittedModel = days_glm_step, n = 500)
plot(check_gamma_model)

# Get % variance explained only by family and other effects
mod_rsq <- rsq.partial(days_glm_step,adj = T,type="v")
round(sum(mod_rsq$partial.rsq),3)
mod_rsq$variable
round(mod_rsq$partial.rsq,3)

# Analyse
summary(days_glm_step)
anova(days_glm_step,test="F")
drop1(days_glm_step,test="F")
days_aov <- aov(days_glm_step,test="F")
days_tuk <- TukeyHSD(x=days_aov, 'family', conf.level=0.95)
plot(days_tuk , las=1 , col="brown")
plot(allEffects(days_glm_step))

# Make significance groups...
days.levels <- days_tuk[["family"]][,4]
days.labels <- data.frame(multcompLetters(days.levels)['Letters'])
days.labels$variable <- "Days to Maturity"
for(i in 1:nrow(days.labels)){
  days.labels$y_pos[i] <- max(pheno_dd[pheno_dd$family == rownames(days.labels)[i],"Days"])+0.1*(max(pheno_dd[pheno_dd$family == rownames(days.labels)[i],"Days"]))
}
days.labels$family <- rownames(days.labels)

####################################################################
# Compare length between families...
tmp <- na.omit(pheno_dd[,c("Mature_length","temp","family","DOB")])
size_glm <- glm(Mature_length~family+temp+DOB,data=tmp,family="gaussian")
size_glm_step <- step(size_glm)

# Check model...
shapiro.test(size_glm_step$residuals) # All good
#plot(size_glm_step)
check_normal_model <- simulateResiduals(fittedModel = size_glm_step, n = 500)
plot(check_normal_model)

summary(size_glm_step)
anova(size_glm_step,test="F")
drop1(size_glm_step,test="F")
size_aov <- aov(size_glm_step,test="F")
size_tuk <- TukeyHSD(x=size_aov, 'family', conf.level=0.95)
plot(size_tuk , las=1 , col="brown")
plot(allEffects(size_glm_step))

# Get % variance explained only by family and other effects
mod_rsq <- rsq.partial(size_glm_step,adj = T,type="v")
round(sum(mod_rsq$partial.rsq),3)
mod_rsq$variable
round(mod_rsq$partial.rsq,3)

# Make significance groups...
size.levels <- size_tuk[["family"]][,4]
size.labels <- data.frame(multcompLetters(size.levels)['Letters'])
size.labels$variable <- "Mature Length (mm)"
for(i in 1:nrow(size.labels)){
  size.labels$y_pos[i] <- max(pheno_dd[pheno_dd$family == rownames(size.labels)[i],"Mature_length"])+0.1*(max(pheno_dd[pheno_dd$family == rownames(size.labels)[i],"Mature_length"]))
}
size.labels$family <- rownames(size.labels)

# Rbind
tukey_labs <- rbind(days.labels,size.labels)

# Plot each...
plot_dd <- melt(pheno_dd,
                measure.vars=c("Days","Mature_length"),
                id.vars=c("family"))
plot_dd$variable <- gsub("Days","Days to Maturity",plot_dd$variable)
plot_dd$variable <- gsub("Mature_length","Mature Length (mm)",plot_dd$variable)

male_pheno <- ggplot(plot_dd,aes(family,value,fill=family))+
  theme_bw()+
  geom_violin(alpha=0.5,draw_quantiles = c(0.25,0.5,0.75),show.legend = F)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        axis.text = element_text(size=16),
        strip.background = element_blank(),
        strip.text = element_text(size=16),
        axis.text.x = element_text(size=16,angle=45,hjust=1))+
  labs(y="Phenotype")+
  facet_wrap(~variable,ncol=1,nrow=2,scales = "free_y",strip.position = "right")+
  geom_text(data=tukey_labs,aes(x=family,y=y_pos,label=Letters),size=6)+
  scale_fill_manual(values=wes_palette("Rushmore1")[c(1,3:5)])

# Are male traits correlated?
cor.test(log10(pheno_dd$Days),pheno_dd$Mature_length,method = "pearson")

# Compare length between families...
male_glm <- glm(Mature_length~Days+family+temp+DOB,data=na.omit(pheno_dd),family="gaussian")
male_glm_step <- step(male_glm)
drop1(male_glm_step,test="F")
plot(allEffects(male_glm_step))

################################################################################################
#### Females #####
females <- read_cross2("data/qtl_cross_females_informative.yaml")

# Adjust phenos
pheno_dd <- data.frame(inds=rownames(females$covar),
                       family=females$covar$cross)
pheno_dd <- cbind(pheno_dd,females$pheno)

# Filter out bad interbroods...
pheno_dd[pheno_dd$interbrood < 15 | 
           pheno_dd$interbrood > 40, c("interbrood",
                                       "first_brood_size",
                                       "dry_offspring_weight_1",
                                       "dry_offspring_weight_2")]<-NA

# Get residual first brood size
size_mod <- lm(first_brood_size~size_at_first_brood,data=pheno_dd)
for(i in 1:nrow(pheno_dd)){
  pheno_dd$predicted_first_brood_size[i] <- predict(size_mod)[rownames(pheno_dd)[i]]
}
pheno_dd$first_brood_size_resid <- pheno_dd$first_brood_size-pheno_dd$predicted_first_brood_size
pheno_dd <- pheno_dd[,c("inds","family","age_first_brood","size_at_first_brood","first_brood_size_resid","interbrood","dry_offspring_weight_1")]

# First run PCA to examine tightly associated variables
# PCA
female_pca <- prcomp(na.omit(pheno_dd[,3:ncol(pheno_dd)]),scale=T,center=T)
pca_sum <- summary(female_pca)
pca_sum <- pca_sum$importance
pca_out <- female_pca$rotation
colnames(pca_out) <- paste0(colnames(pca_out),": ",round(pca_sum[2,],3)*100,"%")

# Plot the females
scores <- data.frame(female_pca$x)
scores$family <- NA
for(i in 1:nrow(scores)){
  scores$family[i] <- females$covar[rownames(scores)[i],"cross"]
}

# Save
write.table(pca_out,
            "tables/TableSX_female_PCA_covariance.txt",
            quote = F,sep = "\t")

# Add rearing
pheno_dd$temp <- NA
pheno_dd$DOB <- NA
for(i in 1:nrow(pheno_dd)){
  if(length(na.omit(female_rearing[female_rearing$vcf_id == pheno_dd$inds[i],"mean_temp"]))==1){
    pheno_dd$temp[i] <- na.omit(female_rearing[female_rearing$vcf_id == pheno_dd$inds[i],"mean_temp"])
  } else {
    pheno_dd$temp[i] <- NA 
  }
  if(length(na.omit(female_rearing[female_rearing$vcf_id == pheno_dd$inds[i],"DOB"]))==1){
    pheno_dd$DOB[i] <- na.omit(female_rearing[female_rearing$vcf_id == pheno_dd$inds[i],"DOB"])
  } else {
    pheno_dd$DOB[i] <- NA 
  }
}

# Correlations between phenotypes and covar
female_phenotypes_of_interest <- c("age_first_brood","size_at_first_brood","first_brood_size_resid","interbrood","dry_offspring_weight_1")
female_pheno_covar <- rcorr(as.matrix(pheno_dd[,c(female_phenotypes_of_interest,"temp","DOB")]),type = "spearman")
female_pheno_covar$r
female_pheno_covar$P
female_covar_cor_dd <- melt(female_pheno_covar$r)
female_covar_cor_dd$p <- melt(female_pheno_covar$P)$value
female_covar_cor_dd <- na.omit(female_covar_cor_dd)

# Filter for comps of interest
female_covar_cor_dd<- female_covar_cor_dd[female_covar_cor_dd$Var1 %in% c("DOB","temp") & female_covar_cor_dd$Var1 %in% female_phenotypes_of_interest |
                                            female_covar_cor_dd$Var2 %in% c("DOB","temp") & female_covar_cor_dd$Var1 %in% female_phenotypes_of_interest,]
female_covar_cor_dd$fdr <- p.adjust(female_covar_cor_dd$p,method = "fdr")
female_covar_cor_dd$signif <- NA
female_covar_cor_dd[female_covar_cor_dd$fdr < 0.05,"signif"] <- "*"
female_covar_cor_dd[female_covar_cor_dd$fdr < 0.01,"signif"] <- "**"
female_covar_cor_dd[female_covar_cor_dd$fdr < 0.001,"signif"] <- "***"
for(col in c("value","p","fdr")){
  female_covar_cor_dd[,col] <- round(female_covar_cor_dd[,col],3)
}
female_covar_cor_dd

# Merge with males and save
female_covar_cor_dd$sex <- "female"
male_covar_cor_dd$sex <- "male"
all_covar_cor_dd <- rbind(male_covar_cor_dd,female_covar_cor_dd)
colnames(all_covar_cor_dd) <- c("Phenotype","Rearing Variable","rho","p","fdr","significance","sex")
all_covar_cor_dd <- all_covar_cor_dd[order(all_covar_cor_dd$`Rearing Variable`),]
all_covar_cor_dd$Phenotype <- rep(c("Age at maturity (M)","Size at maturity (M)","Age at first brood (F)","Size at first brood (F)","First brood size (F)","Interbrood Period (F)","Offspring weight (F)"),2)
write.csv(all_covar_cor_dd,
          "tables/TableSX_rearing_effect_correlations.csv",
          row.names = F,quote = F)


# Compare age between families...
hist(log10(pheno_dd$age_first_brood))
tmp <- na.omit(pheno_dd[,c("age_first_brood","DOB","temp","family")])
age_glm <- glm(age_first_brood~family+DOB+temp,
               data=tmp,
               family = Gamma(link="log"))
age_glm_step <- step(age_glm)

# Test model
hist(age_glm_step$residuals)
shapiro.test(age_glm_step$residuals)
#plot(age_glm_step)
check_model <- simulateResiduals(fittedModel = age_glm_step, n = 500)
plot(check_model)

summary(age_glm_step)
anova(age_glm_step,test="F")
drop1(age_glm_step,test="F")

# Get % variance explained only by family and other effects
mod_rsq <- rsq.partial(age_glm_step,adj = T,type="v")
round(sum(mod_rsq$partial.rsq),3)
mod_rsq$variable
round(mod_rsq$partial.rsq,3)

# Just plot a correlation
cor.test(tmp$age_first_brood,tmp$temp,method="spearman")

# Compare length between families...
hist(pheno_dd$size_at_first_brood)
tmp <- na.omit(pheno_dd[,c("size_at_first_brood","family","temp","DOB")])
fem_size_glm <- glm(size_at_first_brood~family+temp+DOB,
                    data=tmp,
                    family=Gamma(link="log"))
fem_size_glm_step <- step(fem_size_glm)

# Assess
shapiro.test(fem_size_glm_step$residuals)
#plot(fem_size_glm)
check_model <- simulateResiduals(fittedModel = fem_size_glm_step, n = 500)
plot(check_model)

summary(fem_size_glm_step)
anova(fem_size_glm_step,test="F")
drop1(fem_size_glm_step,test="F")

# Get % variance explained only by family and other effects
mod_rsq <- rsq.partial(fem_size_glm_step,adj = T,type="v")
round(sum(mod_rsq$partial.rsq),3)
mod_rsq$variable
round(mod_rsq$partial.rsq,3)

fem_size_aov <- aov(fem_size_glm_step)
fem_size_tuk <- TukeyHSD(x=fem_size_aov, 'family', conf.level=0.95)
plot(fem_size_tuk , las=1 , col="brown")
plot(allEffects(fem_size_glm_step))

# Make significance groups...
fem_size.levels <- fem_size_tuk[["family"]][,4]
fem_size.labels <- data.frame(multcompLetters(fem_size.levels)['Letters'])
fem_size.labels$variable <- "Size at first brood (mm)"
for(i in 1:nrow(days.labels)){
  fem_size.labels$y_pos[i] <- max(pheno_dd[pheno_dd$family == rownames(fem_size.labels)[i],"size_at_first_brood"])+0.1*(max(pheno_dd[pheno_dd$family == rownames(fem_size.labels)[i],"size_at_first_brood"]))
}
fem_size.labels$family <- rownames(fem_size.labels)

# Compare first brood size between families...
pheno_tmp <- na.omit(pheno_dd[,c("family","first_brood_size_resid","temp","DOB")])
hist(pheno_tmp$first_brood_size_resid)
first_brood_glm <- glm(first_brood_size_resid~family+temp+DOB,
                       data=pheno_tmp,
                       family="gaussian")
first_brood_glm_step <- step(first_brood_glm)

# Assess model
shapiro.test(first_brood_glm_step$residuals)
#plot(first_brood_glm_step)
check_model <- simulateResiduals(fittedModel = first_brood_glm_step, n = 500)
plot(check_model)

# Get res
summary(first_brood_glm_step)
anova(first_brood_glm_step,test="F")
drop1(first_brood_glm_step,test="F")
first_brood_aov <- aov(first_brood_glm)
first_brood_tuk <- TukeyHSD(x=first_brood_aov, 'family', conf.level=0.95)
plot(first_brood_tuk , las=1 , col="brown")
plot(allEffects(first_brood_glm_step))

# Get % variance explained only by family and other effects
mod_rsq <- rsq.partial(first_brood_glm_step,adj = T,type="v")
round(sum(mod_rsq$partial.rsq),3)
mod_rsq$variable
round(mod_rsq$partial.rsq,3)

# Make significance groups...
first_brood.levels <- first_brood_tuk[["family"]][,4]
first_brood.labels <- data.frame(multcompLetters(first_brood.levels)['Letters'])
first_brood.labels$variable <- "First brood size (residual)"
for(i in 1:nrow(days.labels)){
  first_brood.labels$y_pos[i] <- max(pheno_dd[pheno_dd$family == rownames(first_brood.labels)[i],"size_at_first_brood"])+0.1*(max(pheno_dd[pheno_dd$family == rownames(first_brood.labels)[i],"size_at_first_brood"]))
}
first_brood.labels$family <- rownames(first_brood.labels)

# Compare interbrood between families...
hist(log(pheno_dd$interbrood))
tmp <- na.omit(pheno_dd[,c("interbrood","family","temp","DOB")])
interbrood_glm <- glm(interbrood~family+temp+DOB,
                      data=tmp,
                      family=Gamma(link="log"))
interbrood_glm_step <- step(interbrood_glm)

# Assess
shapiro.test(interbrood_glm_step$residuals)
#plot(interbrood_glm_step)
check_model <- simulateResiduals(fittedModel = interbrood_glm_step, n = 500)
plot(check_model)
testQuantiles(interbrood_glm_step) 

# Get results
summary(interbrood_glm_step)
anova(interbrood_glm_step,test="F")
drop1(interbrood_glm_step,test="F")

# Get % variance explained only by family and other effects
mod_rsq <- rsq.partial(interbrood_glm_step,adj = T,type="v")
round(sum(mod_rsq$partial.rsq),3)
mod_rsq$variable
round(mod_rsq$partial.rsq,3)

# Compare dry offspring weight between families...
hist((pheno_dd$dry_offspring_weight_1))
hist(log(pheno_dd$dry_offspring_weight_1))

# Which most normal
shapiro.test(pheno_dd$dry_offspring_weight_1)
shapiro.test(log(pheno_dd$dry_offspring_weight_1))

tmp <- na.omit(pheno_dd[,c("dry_offspring_weight_1","DOB","temp","family")])

# Remove bum individuals
off_weight_glm <- glm(dry_offspring_weight_1~family+DOB+temp,
                      data=tmp,
                      family = Gamma(link="log"))
off_weight_glm_step <- step(off_weight_glm)

# Test model
hist(off_weight_glm_step$residuals)
shapiro.test(off_weight_glm_step$residuals)
#plot(off_weight_glm_step)
check_model <- simulateResiduals(fittedModel = off_weight_glm_step, n = 500)
plot(check_model)

summary(off_weight_glm_step)
anova(off_weight_glm_step,test="F")
drop1(off_weight_glm_step,test="F")

# Get % variance explained only by family and other effects
mod_rsq <- rsq.partial(off_weight_glm_step,adj = T,type="v")
round(sum(mod_rsq$partial.rsq),3)
mod_rsq$variable
round(mod_rsq$partial.rsq,3)

# Just plot a correlation
cor.test(tmp$dry_offspring_weight_1,tmp$DOB,method="spearman")
#ggplot(tmp,aes(family,dry_offspring_weight_1))+geom_violin(draw_quantiles = 0.5)

# Get tukey labels for these...
off_weight_aov <- aov(off_weight_glm_step)
off_weight_tuk <- TukeyHSD(x=off_weight_aov, 'family', conf.level=0.95)
plot(off_weight_tuk , las=1 , col="brown")
plot(allEffects(off_weight_glm_step))

# Make significance groups...
off_weight.levels <- off_weight_tuk[["family"]][,4]
off_weight.labels <- data.frame(multcompLetters(off_weight.levels)['Letters'])
off_weight.labels$variable <- "First brood offspring weight (g)"
for(i in 1:nrow(days.labels)){
  off_weight.labels$y_pos[i] <- max(na.omit(pheno_dd[pheno_dd$family == rownames(off_weight.labels)[i],"dry_offspring_weight_1"]))+0.1*(max(na.omit(pheno_dd[pheno_dd$family == rownames(off_weight.labels)[i],"dry_offspring_weight_1"])))
}
off_weight.labels$family <- rownames(off_weight.labels)

# Rbind
tukey_labs <- rbind(fem_size.labels,first_brood.labels,off_weight.labels)

################
# Plot each...
plot_dd <- melt(pheno_dd,
                measure.vars=female_phenotypes_of_interest,
                id.vars=c("family"))
plot_dd <- plot_dd[plot_dd$variable != "dry_offspring_weight_2",]
plot_dd$variable <- gsub("age_first_brood","Age at first brood (Days)",plot_dd$variable)
plot_dd$variable <- gsub("size_at_first_brood","Size at first brood (mm)",plot_dd$variable)
plot_dd$variable <- gsub("first_brood_size_resid","First brood size (residual)",plot_dd$variable)
plot_dd$variable <- gsub("interbrood","Interbrood (Days)",plot_dd$variable)
plot_dd$variable <- gsub("dry_offspring_weight_1","First brood offspring weight (g)",plot_dd$variable)


female_pheno <- ggplot(plot_dd,aes(family,value,fill=family))+
  theme_bw()+
  geom_violin(alpha=0.5,draw_quantiles = c(0.25,0.5,0.75),show.legend = F)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        axis.text = element_text(size=16),
        strip.background = element_blank(),
        strip.text = element_text(size=11),
        axis.text.x = element_text(size=16,angle=45,hjust=1))+
  labs(y="Phenotype")+
  facet_wrap(~variable,ncol=2,nrow=3,scales = "free_y",strip.position = "right")+
  geom_text(data=tukey_labs,aes(x=family,y=y_pos,label=Letters),size=6)+
  scale_fill_manual(values=wes_palette("Rushmore1")[c(1,3:5)])

# Are female traits correlated?
female_cor <- rcorr(as.matrix(pheno_dd[,female_phenotypes_of_interest]),type = "spearman")
female_cor$r
cor_dd <- melt(female_cor$r)
cor_dd$p <- melt(female_cor$P)$value
cor_dd <- na.omit(cor_dd)
cor_dd$fdr <- p.adjust(cor_dd$p,method = "fdr")
cor_dd$signif <- NA
cor_dd[cor_dd$fdr < 0.05,"signif"] <- "*"
cor_dd[cor_dd$fdr < 0.01,"signif"] <- "**"
cor_dd[cor_dd$fdr < 0.001,"signif"] <- "***"

# Visualise
library(corrplot)
female_cor
corrplot(female_cor$r, method = "ellipse",type="upper",order="hclust",p.mat = female_cor$P,insig = "p-value", sig.level = -1)

# Write table
colnames(cor_dd) <- c("Phenotype 1","Phenotype 2","Spearmans rho","p","fdr","significance")
write.csv(cor_dd,
          "tables/TableSX_female_phenotype_correlations.csv",
          row.names = F,quote = F)

# Merge plots and expprt
all_phenos <- plot_grid(male_pheno,female_pheno,
                        ncol=2,nrow=1,labels = "AUTO",label_size = 20,rel_widths = c(1,2))
pdf("figs/FigureX_Phenotype_family_distributions.pdf",width=8,height=9)
all_phenos
dev.off()
