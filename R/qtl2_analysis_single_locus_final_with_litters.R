#########################################################
# Run QTL mapping with r/qtl2
lib<-c("data.table","ggplot2","qtl2","parallel","viridis","cowplot","see","RColorBrewer")
lapply(lib,library,character.only=T)

# Source
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

# These are the yamls to run...
female_yaml <- 'data/qtl_cross_females_informative.yaml'
male_yaml <-  'data/qtl_cross_males_informative.yaml'
sex_yaml <- 'data/qtl_cross_SEX_informative.yaml'

#########################################################
# Scan plotting function
qtl_plot <- function(scan_results=NULL,
                     map=NULL,
                     colour_vec=c("slateblue", "violetred","forestgreen","gold2"),
                     cutoff=NULL){
  
  # Make plotting map
  map2 <- data.frame(rbindlist(lapply(1:length(map),function(x){
    return(data.frame(chr=names(map)[x],
                      pos=map[[x]]))
  })))
  
  # Fetch results
  to_plot <- cbind(map2,scan_results)
  to_plot$chrN <- gsub("chr","",to_plot$chr)
  to_plot$chr_F <- factor(to_plot$chrN,levels=1:23)
  
  # Reverse positions on 'backwards' linkage groups to reflect the genome
  to_reverse <- paste0(c(1,2,5,6,7,8,12,13,14,16,17,18,19,20,22,23))
  for(i in to_reverse){
    max_tmp = max(to_plot[as.integer(to_plot$chrN) == i,"pos"])
    to_plot[as.integer(to_plot$chrN) == i,"pos"] <- max_tmp - to_plot[as.integer(to_plot$chrN) == i,"pos"]
  }
  
  
  # Melt
  to_plot_melt <- reshape2::melt(to_plot,
                                 value.name="LOD",
                                 id.vars=c("chr_F","pos"),
                                 measure.vars=colnames(scan_results))
  
  
  # Make basic plot
  p1 <- ggplot(to_plot_melt,aes(x=pos,y=LOD,colour=variable))+
    geom_line()+
    facet_wrap(~chr_F,nrow=1,strip.position = "bottom",scales = "free_x")+
    theme_minimal()+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size=20),
          strip.text = element_text(size=18),
          axis.text.y = element_text(size=18),
          panel.spacing.x=unit(0, "lines"),
          legend.title = element_blank(),
          legend.text = element_text(size=14))+
    labs(x="Genome Position (Chr)",y="LOD")+
    scale_colour_manual(breaks = colnames(scan_results),
                        values = colour_vec[1:ncol(scan_results)])+
    geom_hline(yintercept = cutoff,linetype="dashed")
  
}

# Prep female cross phenotypes
prep_female_cross <- function(cross_tmp){
  
  ### Edit phenotypes
  phenos <- data.frame(cross_tmp$pheno)
  
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
  #hist(phenos$first_brood_size_resid)
  
  # Log-transform all values, test, check hists
  to_transform <- c("age_first_brood","size_at_first_brood","interbrood","dry_offspring_weight_1")
  for(i in to_transform){
    phenos[,i] <- log(phenos[,i])
    #print(shapiro.test(phenos[,i]))
  }
  
  
  # Save
  cross_tmp$pheno<-phenos[,c("age_first_brood","size_at_first_brood","first_brood_size_resid","interbrood","dry_offspring_weight_1")]
  colnames(cross_tmp$pheno) <- c("Age at first brood","Size at first brood","First brood size","Interbrood","First brood offspring weight")
  return(cross_tmp)
}

# Get QTL effects
qtl_effects_plot <- function(cross_in,phenotype,qtl_chr,qtl_pos,qtl_covar,qtl_low,qtl_high){
  
  # Get kinship
  map <- insert_pseudomarkers(cross_in$gmap, step=1)
  pr <- calc_genoprob(cross_in, map, error_prob=0.002,cores=4)
  apr <- genoprob_to_alleleprob(pr)
  
  # Calculate a kinship matrix
  grid<-calc_grid(map = cross_in$gmap, step=1)
  pr_grid <- probs_to_grid(probs = pr, grid = grid)
  
  # Also calculate a loco kinship matrix
  kinship_loco <- calc_kinship(pr_grid, "loco", cores=4)
  
  pheno_tmp <- cross_in$pheno[,phenotype]
  names(pheno_tmp) <- rownames(cross_in$pheno)
  c2eff <- scan1coef(genoprobs = pr[,qtl_chr],
                     pheno = pheno_tmp,
                     kinship = kinship_loco[qtl_chr],
                     addcovar = qtl_covar)
  
  
  # Plot better
  plot_dd <- data.frame(map=unlist(map[qtl_chr]),
                        effect=c(c2eff[,1],c2eff[,2],c2eff[,3]),
                        geno=rep(colnames(c2eff)[1:3],each=length(unlist(map[qtl_chr]))))
  
  # Reverse positions on 'backwards' linkage groups to reflect the genome
  to_reverse <- paste0("chr",c(1,2,5,6,7,8,12,13,14,16,17,18,19,20,22,23))
  
  # Also plot the QTL effects at the specific locus
  g <- maxmarg(pr, map, chr=qtl_chr, pos=qtl_pos , return_char=TRUE)
  qtl_plot_dd <- data.frame(geno = g,
                            pheno = pheno_tmp)
  qtl_effect_rain <- ggplot(data = na.omit(qtl_plot_dd), aes(x = geno, y = pheno, fill = geno)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8,) +
    geom_point(aes(y = pheno,colour=geno), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.5) +
    expand_limits(x = 4) +
    guides(fill = FALSE) +
    guides(color = FALSE) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    coord_flip() +
    theme_bw() +
    theme(axis.title=element_text(size=16),
          axis.text=element_text(size=14),
          legend.text=element_text(size=14),
          legend.position = "top")+
    labs(x="",y=phenotype)
  
  # And the whole chr...
  if(qtl_chr %in% to_reverse){
    max_tmp <- max(plot_dd$map) 
    plot_dd$map <- max_tmp - plot_dd$map
    qtl_pos <- max_tmp - qtl_pos
    qtl_low <-  max_tmp - qtl_low
    qtl_high <-  max_tmp - qtl_high
  }
  qtl_effects <- ggplot(plot_dd,aes(x=map,y=effect,colour=geno))+
    geom_line()+
    theme_bw()+
    labs(y="Effect",x=paste0(qtl_chr," map pos"),colour="")+
    theme(axis.title=element_text(size=16),
          axis.text=element_text(size=14),
          legend.text=element_text(size=14),
          legend.position = "top")+
    geom_vline(xintercept=qtl_pos)+
    scale_color_brewer(palette = "Dark2")+
    annotate("rect",xmin=qtl_low,xmax=qtl_high,ymin=-Inf,ymax=Inf,fill="black",alpha=0.15)
  
  # Plot together
  # out <- plot_grid(qtl_effects,qtl_effect_rain,
  #                  ncol=1,nrow=2,rel_heights = c(2,1.5),labels = "AUTO",
  #                  axis = "tblr",align="v")
  return(list(qtl_effects,qtl_effect_rain))
  
}

calculate_map_distances <- function(input_map){
  clean_markers <- lapply(input_map,function(map){
    clean_tmp <- data.frame(chr=sapply(strsplit(names(map),"-"),'[[',1),
                            pos=as.integer(sapply(strsplit(names(map),"-"),'[[',2)),
                            cM=map)
    
    physical_dists <- sort(clean_tmp[clean_tmp$chr == names(sort(table(clean_tmp$chr)))[1],"pos"])
    physical_dist_next <- sapply(1:(length(physical_dists)-1),function(x) physical_dists[x+1]-physical_dists[x])
    cM_next <- sapply(1:(nrow(clean_tmp)-1),function(x) clean_tmp$cM[x+1]-clean_tmp$cM[x])
    return(list(physical_dist_next,cM_next))
  })
  
  all_physical_dists <- unlist(lapply(clean_markers,'[[',1))
  # hist(all_physical_dists)
  
  
  all_cM_dists <- unlist(lapply(clean_markers,'[[',2))
  # hist(all_cM_dists)
  
  return(list(median_physical=median(all_physical_dists,na.rm=T),
              median_cM=median(all_cM_dists,na.rm=T)))
}



#########################################################
# Fetch rearing covariates and amend IDs
female_rearing <- read.csv("data/female_rearing_covariates.csv",header=T)
male_rearing <- read.csv("data/male_rearing_covariates.csv",header=T)
metadata <- read.table("data/qtl_cross_ID_metadata_FINAL.tsv",header=T)

# Add in IDs
female_rearing$vcf_id <- NA
for(i in 1:nrow(female_rearing)){
  print(i)
  if(length(metadata[metadata$TRUE_ID == female_rearing$ID[i],"VCF_ID"]) != 0){
    female_rearing$vcf_id[i] <- metadata[metadata$TRUE_ID == female_rearing$ID[i],"VCF_ID"]
  }
}
male_rearing$vcf_id <- NA
for(i in 1:nrow(male_rearing)){
  print(i)
  if(length(metadata[metadata$TRUE_ID == male_rearing$ID[i],"VCF_ID"]) != 0){
    male_rearing$vcf_id[i] <- metadata[metadata$TRUE_ID == male_rearing$ID[i],"VCF_ID"]
  }
}

#########################################################
# Sex test
# Read in the cross data
cross <- read_cross2(sex_yaml)

# Summarise
summary(cross)

# Set up family as a covariate...
# pull out the cage covariate and turn it into a factor
cross$covar$cross <- NA
for(i in 1:nrow(cross$covar)){
  cross$covar$cross[i] <- metadata[metadata$VCF_ID == rownames(cross$covar)[i],"CROSS"]
}
table(cross$covar$cross)
family <- factor(cross$covar$cross)

# use model.matrix to create a numeric matrix of 0/1 covariates
family_matrix <- model.matrix( ~ family)

# omit the first column; drop=FALSE ensures it stays as a matrix
family_matrix <- family_matrix[,-1, drop=FALSE]
rownames(family_matrix) <- rownames(cross$covar)

# Process genotype probabilities etc.
# Insert pseudomarkers
map <- insert_pseudomarkers(cross$gmap, step=1)
pr <- calc_genoprob(cross, map, error_prob=0.002,cores=4)
apr <- genoprob_to_alleleprob(pr)

# Calculate a kinship matrix
grid<-calc_grid(map = cross$gmap, step=1)
pr_grid <- probs_to_grid(probs = pr, grid = grid)
kinship <- calc_kinship(pr_grid)

# Also calculate a loco kinship matrix
kinship_loco <- calc_kinship(pr_grid, "loco", cores=4)

######
# Just plot a basic scan first with family covariate
out <- scan1(pr, cross$pheno,cores=4,addcovar = family_matrix,model="binary")

# Find the peak
sex_peak <- find_peaks(out, map, threshold=55,drop = 1.5)

# Reverse these to reflect new map position
max(map$chr12) - sex_peak$pos
max(map$chr12) - sex_peak$ci_lo
max(map$chr12) - sex_peak$ci_hi

# Plot just the sex chr
sex_chr_effect <- qtl_effects_plot(cross_in = cross,
                                   phenotype = "sex",
                                   qtl_chr = "chr12",
                                   qtl_pos = sex_peak$pos,
                                   qtl_low = sex_peak$ci_lo,
                                   qtl_high=sex_peak$ci_hi,
                                   qtl_covar = family_matrix)

# Plot
sex_qtl <- qtl_plot(out,map)+
  theme(legend.position = "none")
pdf("figs/FigureSX_Sex_QTL.pdf",width = 10,height=6)
plot_grid(sex_qtl,
          plot_grid(plotlist = sex_chr_effect,ncol=2,axis = "tblr",align="h",labels=c("B","C")),
          ncol = 1,axis = "tblr",align="v",labels=c("A"))
dev.off()

#########################################################
##### Run females #####
# Read in the cross data
cross <- read_cross2(female_yaml)

# Summarise
summary(cross)
table(cross$covar$cross)

#### Summarise marker set #####
# Distance between markers...
fully_informative_marker_dists <- calculate_map_distances(cross$gmap)

# Markers per chrom
chr_markers <- matrix(ncol=2,nrow=length(cross$gmap))
chr_markers[,1] <- paste0("chr",1:23)
for(i in 1:nrow(chr_markers)){
  chr_markers[i,2] <- length(cross$gmap[[i]])
}
chr_markers <- data.frame(chr_markers)
chr_markers$X2 <- as.integer(chr_markers$X2)
chr_markers[order(-chr_markers$X2),]
median(chr_markers$X2)

###############################

# Insert pseudomarkers
map <- insert_pseudomarkers(cross$gmap, step=1)
pr <- calc_genoprob(cross, map, error_prob=0.002,cores=4)
apr <- genoprob_to_alleleprob(pr)

# Calculate a kinship matrix
grid <- calc_grid(map = cross$gmap, step=1)
pr_grid <- probs_to_grid(probs = pr, grid = grid)
kinship <- calc_kinship(pr_grid)

# Save female kinship for power analyses...
write.table(kinship,"outputs/qtl2_kinship_matrix_females.txt",
            sep="\t",quote=F)

# Also calculate a loco kinship matrix
kinship_loco <- calc_kinship(pr_grid, "loco", cores=4)

# Adjust phenos
cross <- prep_female_cross(cross)

# Set up family and rearing conditions as  covariate...
family <- factor(cross$covar$cross)

# use model.matrix to create a numeric matrix of 0/1 covariates
family_matrix <- model.matrix( ~ family)

# omit the first column; drop=FALSE ensures it stays as a matrix
family_matrix <- family_matrix[,-1, drop=FALSE]
rownames(family_matrix) <- rownames(cross$covar)

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

######
# First scans...

# Just plot a basic scan first with family covariate
out_add <- scan1(apr,kinship = kinship_loco, cross$pheno,cores=4,addcovar = family_rear_covar)

#######
# Permutation Test
# Additive
if(!(file.exists("outputs/female_qtl2_permutations_additive.rds"))){
  operm <- scan1perm(apr, cross$pheno,kinship_loco,addcovar = family_rear_covar, n_perm=1000,
                     chr_lengths=chr_lengths(map),cores=6)
  saveRDS(operm,"outputs/female_qtl2_permutations_additive.rds")
}

# Read back in
operm_additive <- readRDS("outputs/female_qtl2_permutations_additive.rds")

# Return 5% cut-offs
summary(operm_additive)
for(i in 1:ncol(cross$pheno)){
  print(find_peaks(out_add, map, threshold=summary(operm_additive)[,i],drop = 1.5))
}
# Or others
summary(operm_additive, alpha=c(0.1,0.05,0.01))

# Plots
female_add <- qtl_plot(out_add,map,
                       colour_vec=c("slateblue", "violetred","forestgreen","gold2","red2","blue2"))+
  theme(legend.position = "right")+
  scale_color_brewer(palette="Dark2")

# Add all cutoffs
for(i in 1:length(summary(operm_additive)[1,])){
  female_add <- female_add + geom_hline(yintercept=summary(operm_additive)[1,i],linetype="dashed",colour=brewer.pal(length(summary(operm_additive)[1,]),"Dark2")[i])
}

##### Whole-data QTL effects ####
# Get phenotype effects
female_peaks <- find_peaks(out_add, map, threshold=max(summary(operm_additive)),drop = 1.5)
# Size QTL
size_qtl_effects <- qtl_effects_plot(cross,female_peaks$lodcolumn[1],female_peaks$chr[1],female_peaks$pos[1],family_rear_covar,qtl_low =female_peaks$ci_lo[1],qtl_high = female_peaks$ci_hi[1])

# Offspring weight QTL
offspring_weight_qtl_effects <-  qtl_effects_plot(cross,female_peaks$lodcolumn[2],female_peaks$chr[2],female_peaks$pos[2],family_rear_covar,qtl_low =female_peaks$ci_lo[2],qtl_high = female_peaks$ci_hi[2])

# Merge these with original QTL plot to get full single scan results...
pdf("figs/Figure2_female_qtl_scans.pdf",width=10,height=10)
plot_grid(female_add+theme(legend.position ="top"),
          plot_grid(plot_grid(plotlist = offspring_weight_qtl_effects,ncol=1,nrow=2,axis = "tblr",align="v",labels=c("B","C"),label_size=20),
                    plot_grid(plotlist = size_qtl_effects,ncol=1,nrow=2,axis = "tblr",align="v",labels=c("D","E"),label_size=20)),
          ncol=1,nrow=2,rel_heights=c(1,1.5),labels = "A",label_size=20)
dev.off()


#### PVE of QTL... ####
summary(operm_additive)
calc_pve <- function(LOD, n.ind) {1 - 10 ^ ( - (2 / n.ind) * LOD) }
# Offspring weight
calc_pve(4.493234,nrow(cross$pheno))
# Reverse position to match map
max(map$chr19 - female_peaks$pos[2])
max(map$chr19 - female_peaks$ci_lo[2])
max(map$chr19 - female_peaks$ci_hi[2])

# Size
calc_pve(3.552253,nrow(cross$pheno))
# Reverse position to match map
max(map$chr22 - female_peaks$pos[1])
max(map$chr22 - female_peaks$ci_lo[1])
max(map$chr22 - female_peaks$ci_hi[1])

#########################################################
##### Repeat mappings but within each family ####
families <- unique(cross$covar$cross)

# Run over all families and export final figures...
family_figs <- lapply(families,function(family_id){
  print(family_id)
  
  # Read in cross
  cross_tmp <- read_cross2(paste0("data/qtl_cross_",family_id,"_females_informative.yaml"))
  
  # Prep female cross
  cross_tmp <- prep_female_cross(cross_tmp)
  
  # Get marker distances for these...
  family_informative_marker_dists <- calculate_map_distances(cross_tmp$gmap)
  
  # Insert pseudomarkers
  map <- insert_pseudomarkers(cross_tmp$gmap, step=1)
  pr <- calc_genoprob(cross_tmp, map, error_prob=0.002,cores=4)
  apr <- genoprob_to_alleleprob(pr)
  
  # Calculate a kinship matrix
  grid<-calc_grid(map = cross_tmp$gmap, step=1)
  pr_grid <- probs_to_grid(probs = pr, grid = grid)
  kinship <- calc_kinship(pr_grid)
  
  # Also calculate a loco kinship matrix
  kinship_loco <- calc_kinship(pr_grid, "loco", cores=4)
  
  # and rearing...
  phenos <- cross_tmp$pheno
  rear_mat <- matrix(nrow = nrow(phenos),ncol=2)
  rownames(rear_mat) <- rownames(phenos)
  colnames(rear_mat) <- c("DOB","mean_temp")
  for(i in 1:nrow(rear_mat)){
    if(length(na.omit(female_rearing[female_rearing$vcf_id == rownames(rear_mat)[i],"DOB"])) != 0){
      rear_mat[i,1] <- na.omit(female_rearing[female_rearing$vcf_id == rownames(rear_mat)[i],"DOB"])
    }
    if(length(na.omit(female_rearing[female_rearing$vcf_id == rownames(rear_mat)[i],"mean_temp"])) != 0){
      rear_mat[i,2] <- na.omit(female_rearing[female_rearing$vcf_id == rownames(rear_mat)[i],"mean_temp"])
    }}
  
  # Merge
  family_rear_covar <- rear_mat
  
  ######
  # First scans...
  
  # Just plot a basic scan first with family covariate
  out_add <- scan1(apr,kinship = kinship_loco, cross_tmp$pheno,cores=4,addcovar = family_rear_covar)
  
  #######
  # Permutation Test
  # Additive
  if(!(file.exists(paste0("outputs/female_qtl2_permutations_additive_",family_id,".rds")))){
    operm <- scan1perm(apr, cross_tmp$pheno,kinship_loco,addcovar = family_rear_covar, n_perm=1000,
                       chr_lengths=chr_lengths(map),cores=detectCores()-1)
    saveRDS(operm,paste0("outputs/female_qtl2_permutations_additive_",family_id,".rds"))
  }
  
  # Read back in
  operm_additive <- readRDS(paste0("outputs/female_qtl2_permutations_additive_",family_id,".rds"))
  
  # Get cutoffs
  summary(operm_additive, alpha=c(0.1,0.05,0.01))
  print(find_peaks(out_add, map, threshold=summary(operm_additive)[,1], drop = 1.5))
  
  # Retain and plot only variables for which we have significant peaks
  to_retain <- unique(find_peaks(out_add, map, threshold=summary(operm_additive)[,1], drop=1.5)[,2])
  
  # Re-Run model with subset variables
  if(length(to_retain > 0)){
    
    # Plots
    female_add_tmp <- qtl_plot(out_add,map,c("slateblue", "violetred","forestgreen","gold2","red2"))+
      ggtitle(label=family_id)+
      theme(legend.position = "top",
            title=element_text(size=16))+
      scale_color_brewer(palette = "Dark2")
    
    # Add all cutoffs
    for(i in 1:length(summary(operm_additive)[1,])){
      female_add_tmp <- female_add_tmp + geom_hline(yintercept=summary(operm_additive)[1,i],linetype="dashed",colour=brewer.pal(length(summary(operm_additive)[1,]),"Dark2")[i])
    }
    
  } else {
    female_add_tmp <- qtl_plot(out_add,map,c("slateblue", "violetred","forestgreen","gold2","red2"))+
      ggtitle(label=family_id)+
      theme(legend.position = "top",
            title=element_text(size=16))+
      scale_color_brewer(palette = "Dark2")
    
    # Add all cutoffs
    for(i in 1:length(summary(operm_additive)[1,])){
      female_add_tmp <- female_add_tmp + geom_hline(yintercept=summary(operm_additive)[1,i],linetype="dashed",colour=brewer.pal(length(summary(operm_additive)[1,]),"Dark2")[i])
    }
  }
  
  return(list(out_add,operm_additive,female_add_tmp,family_rear_covar,family_informative_marker_dists))
})

# Get marker Ns 
family_marker_N <- sapply(families,function(family_id){
  cross_tmp <- read_cross2(paste0("data/qtl_cross_",family_id,"_females_informative.yaml"))
  return(length(unlist(cross_tmp$gmap)))
})

# Assemble the marker distances and save to a table
family_dists <- lapply(family_figs,'[[',5)
all_marker_dists <- matrix(nrow = 5,ncol=3)
all_marker_dists[1,2:3] <- unlist(fully_informative_marker_dists)
all_marker_dists[,1] <- c(1220,family_marker_N)
for(i in 1:length(family_dists)){
  all_marker_dists[i+1,2:3] <- unlist(family_dists[[i]])
}
rownames(all_marker_dists) <- c("Fully informative",paste0(families," informative"))
colnames(all_marker_dists) <- c("Marker N","Median physical distance (kb)","Median genetic distance (cM)")
all_marker_dists[,2] <- round(all_marker_dists[,2]/1000,4)
all_marker_dists[,3] <- round(all_marker_dists[,3],4)

write.table(all_marker_dists,
            "tables/TableSX_marker_physical_cM_distance.txt",
            sep="\t",quote=F)

# Bring together all the family plots
family_qtl <- lapply(family_figs,function(x){return(x[[3]])})
family_legend <- get_legend(family_qtl[[1]])
qtl_family_plots <- plot_grid(family_legend,
                              family_qtl[[1]]+theme(axis.text.x = element_blank(),
                                                    axis.title.x = element_blank(),
                                                    legend.position="none"),
                              family_qtl[[4]]+theme(axis.text.x = element_blank(),
                                                    axis.title.x = element_blank(),
                                                    legend.position="none"),
                              family_qtl[[3]]+theme(axis.text.x = element_blank(),
                                                    axis.title.x = element_blank(),
                                                    legend.position="none"),
                              family_qtl[[2]]+theme(
                                legend.position="none"),
                              ncol=1,align = "v",axis = "tblr",rel_heights = c(1,4,4,4,4))
pdf("figs/FigureSX_single_family_qtl_scans.pdf",width=12,height=10)
qtl_family_plots
dev.off()

#########################################################
##### Get within family QTL figs and PVE #####
# YHF5QM6
# Get cross
YHF5QM6_cross <- prep_female_cross(read_cross2("data/qtl_cross_YHF5QM6_females_informative.yaml"))

# Get models 
YHF5QM6_mod <- family_figs[[which(families=="YHF5QM6")]][[1]]
YHF5QM6_perm <- family_figs[[which(families=="YHF5QM6")]][[2]]
YHF5QM6_covar <- family_figs[[which(families=="YHF5QM6")]][[4]]
YHF5QM6_map <- insert_pseudomarkers(YHF5QM6_cross$gmap, step=1)

# Get cutoffs
summary(YHF5QM6_perm, alpha=c(0.1,0.05,0.01))
YHF5QM6_peaks <- find_peaks(YHF5QM6_mod, YHF5QM6_map, threshold=summary(YHF5QM6_perm)[,i],drop = 1.5)

# Get plots
YHF5QM6_brood_size <- qtl_effects_plot(cross_in=YHF5QM6_cross,
                                       phenotype = "First brood size",
                                       qtl_chr="chr23",
                                       qtl_pos=25.422,
                                       qtl_low=9.164,
                                       qtl_high=57.011,
                                       qtl_covar = YHF5QM6_covar)
YHF5QM6_interbrood <- qtl_effects_plot(cross_in=YHF5QM6_cross,
                                       phenotype = "Interbrood",
                                       qtl_chr="chr12",
                                       qtl_pos=11.774,
                                       qtl_low=0.000,
                                       qtl_high=16.648,
                                       qtl_covar = YHF5QM6_covar)

# Get PVE
# Brood Size
calc_pve(LOD = 3.386414,nrow(YHF5QM6_covar))
# Reverse to match genome
max(YHF5QM6_map$chr23) - YHF5QM6_peaks$pos[1]
max(YHF5QM6_map$chr23) - YHF5QM6_peaks$ci_lo[1]
max(YHF5QM6_map$chr23) - YHF5QM6_peaks$ci_hi[1]

# Interbrood
calc_pve(LOD = 3.405899,nrow(YHF5QM6_covar))
# Reverse to match genome
max(YHF5QM6_map$chr12) - YHF5QM6_peaks$pos[2]
max(YHF5QM6_map$chr12) - YHF5QM6_peaks$ci_lo[2]
max(YHF5QM6_map$chr12) - YHF5QM6_peaks$ci_hi[2]

############### QF4YHM8
# Get cross
QF4YHM8_cross <- prep_female_cross(read_cross2("data/qtl_cross_QF4YHM8_females_informative.yaml"))

# Get models 
QF4YHM8_mod <- family_figs[[which(families=="QF4YHM8")]][[1]]
QF4YHM8_perm <- family_figs[[which(families=="QF4YHM8")]][[2]]
QF4YHM8_covar <- family_figs[[which(families=="QF4YHM8")]][[4]]
QF4YHM8_map <- insert_pseudomarkers(QF4YHM8_cross$gmap, step=1)

# Get cutoffs
summary(QF4YHM8_perm, alpha=c(0.1,0.05,0.01))
QF4YHM8_peaks <- find_peaks(QF4YHM8_mod, QF4YHM8_map, threshold=min(summary(QF4YHM8_perm)), drop = 1.5)


# Get plots
QF4YHM8_interbrood <- qtl_effects_plot(cross_in=QF4YHM8_cross,
                                       phenotype = "Interbrood",
                                       qtl_chr="chr14",
                                       qtl_pos=40.943,
                                       qtl_low=1.579,
                                       qtl_high=42.64,
                                       qtl_covar = QF4YHM8_covar)

plot_grid(plotlist = QF4YHM8_interbrood,
          ncol=1,axis = "tblr",align="v")

# Get PVE
# Interbrood
calc_pve(LOD = 3.558526,nrow(QF4YHM8_covar))
# Reverse to match genome
max(QF4YHM8_map$chr14) - QF4YHM8_peaks$pos[1]
max(QF4YHM8_map$chr14) - QF4YHM8_peaks$ci_lo[1]
max(QF4YHM8_map$chr14) - QF4YHM8_peaks$ci_hi[1]

# Plot all of them together
pdf("figs/FigureSX_per_family_female_qtl_effects.pdf",width=12,height=8)
plot_grid(
  plot_grid(plotlist = YHF5QM6_brood_size,ncol=1,axis = "tblr",align="v"),
  plot_grid(plotlist = YHF5QM6_interbrood,ncol=1,axis = "tblr",align="v"),
  plot_grid(plotlist = QF4YHM8_interbrood,ncol=1,axis = "tblr",align="v"),
  ncol=3,labels = "AUTO",label_size = 20)
dev.off()

#########################################################
##### Run males ####
# Read in the cross data
cross <- read_cross2(male_yaml)

# Summarise
summary(cross)
table(cross$covar$cross)

# log-transform male age
cross$pheno[,"Days"] <- log(cross$pheno[,"Days"])

# Insert pseudomarkers
map <- insert_pseudomarkers(cross$gmap, step=1)
pr <- calc_genoprob(cross, map, error_prob=0.002,cores=4)
apr <- genoprob_to_alleleprob(pr)

# Calculate a kinship matrix
grid<-calc_grid(map = cross$gmap, step=1)
pr_grid <- probs_to_grid(probs = pr, grid = grid)
kinship <- calc_kinship(pr_grid)

# Save male kinship for power analyses...
write.table(kinship,"outputs/qtl2_kinship_matrix_males.txt",
            sep="\t",quote=F)

# Also calculate a loco kinship matrix
kinship_loco <- calc_kinship(pr_grid, "loco", cores=4)

# Set up family and rearing conditions as  covariate...
family <- factor(cross$covar$cross)

# use model.matrix to create a numeric matrix of 0/1 covariates
family_matrix <- model.matrix( ~ family)

# omit the first column; drop=FALSE ensures it stays as a matrix
family_matrix <- family_matrix[,-1, drop=FALSE]
rownames(family_matrix) <- rownames(cross$covar)

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

######
# First scans...

# Just plot a basic scan first with family covariate
cross$pheno <- cross$pheno[,1:2]
colnames(cross$pheno) <- c("Days to maturity","Mature Length")
out_add <- scan1(apr,kinship = kinship_loco, cross$pheno[,1:2],cores=4,addcovar = family_rear_covar)

#####
# Permutation Test
# Additive
if(!(file.exists("outputs/male_qtl2_permutations_additive.rds"))){
  operm <- scan1perm(apr, cross$pheno,kinship_loco,addcovar = family_rear_covar, n_perm=1000,
                     chr_lengths=chr_lengths(map),cores=detectCores()-1)
  saveRDS(operm,"outputs/male_qtl2_permutations_additive.rds")
}

# Read back in
operm_additive <- readRDS("outputs/male_qtl2_permutations_additive.rds")

# Return 5% cut-offs
summary(operm_additive)
for(i in 1:2){
  print(find_peaks(out_add, map, threshold=summary(operm_additive)[,i], drop=1.5))
}
# Or others
summary(operm_additive, alpha=c(0.1,0.05,0.01))

# Plots
male_add <- qtl_plot(out_add,map)+
  theme(legend.position = "top")+
  scale_color_brewer(palette="Dark2")

# Add all the cutoffs
for(i in 1:length(summary(operm_additive)[1,])){
  male_add <- male_add + geom_hline(yintercept=summary(operm_additive)[1,i],linetype="dashed",colour=brewer.pal(length(summary(operm_additive)[1,]),"Dark2")[i])
}

# Save figure
pdf("figs/Figure3_male_traits_qtl_results.pdf",width=12,height=2.5)
male_add
dev.off()

##### Run within male families ######
# Run over all families and export final figures...
families <- unique(cross$covar$cross)
family_figs <- lapply(families,function(family_id){
  print(family_id)
  
  # Read in cross
  cross_tmp <- read_cross2(paste0("data/qtl_cross_",family_id,"_males_informative.yaml"))
  
  # Insert pseudomarkers
  map <- insert_pseudomarkers(cross_tmp$gmap, step=1)
  pr <- calc_genoprob(cross_tmp, map, error_prob=0.002,cores=4)
  apr <- genoprob_to_alleleprob(pr)
  
  # Calculate a kinship matrix
  grid<-calc_grid(map = cross_tmp$gmap, step=1)
  pr_grid <- probs_to_grid(probs = pr, grid = grid)
  kinship <- calc_kinship(pr_grid)
  
  # Also calculate a loco kinship matrix
  kinship_loco <- calc_kinship(pr_grid, "loco", cores=4)
  
  # and rearing...
  phenos <- cross_tmp$pheno
  rear_mat <- matrix(nrow = nrow(phenos),ncol=2)
  rownames(rear_mat) <- rownames(phenos)
  colnames(rear_mat) <- c("DOB","mean_temp")
  for(i in 1:nrow(rear_mat)){
    if(length(na.omit(male_rearing[male_rearing$vcf_id == rownames(rear_mat)[i],"DOB"])) != 0){
      rear_mat[i,1] <- na.omit(male_rearing[male_rearing$vcf_id == rownames(rear_mat)[i],"DOB"])
    }
    if(length(na.omit(male_rearing[male_rearing$vcf_id == rownames(rear_mat)[i],"mean_temp"])) != 0){
      rear_mat[i,2] <- na.omit(male_rearing[male_rearing$vcf_id == rownames(rear_mat)[i],"mean_temp"])
    }}
  
  # Merge
  family_rear_covar <- rear_mat
  
  ######
  # First scans...
  
  # Just plot a basic scan first with family covariate
  cross_tmp$pheno <- cross_tmp$pheno[,1:2]
  cross_tmp$pheno[,"Days"] <- log(cross_tmp$pheno[,"Days"])
  colnames(cross_tmp$pheno) <- c("Days to maturity","Mature Length")
  out_add <- scan1(apr,kinship = kinship_loco, cross_tmp$pheno,cores=4,addcovar = family_rear_covar)
  
  #######
  # Permutation Test
  # Additive
  if(!(file.exists(paste0("outputs/male_qtl2_permutations_additive_",family_id,".rds")))){
    operm <- scan1perm(apr, cross_tmp$pheno,kinship_loco,addcovar = family_rear_covar, n_perm=1000,
                       chr_lengths=chr_lengths(map),cores=detectCores()-1)
    saveRDS(operm,paste0("outputs/male_qtl2_permutations_additive_",family_id,".rds"))
  }
  
  # Read back in
  operm_additive <- readRDS(paste0("outputs/male_qtl2_permutations_additive_",family_id,".rds"))
  
  # Get cutoffs
  summary(operm_additive, alpha=c(0.1,0.05,0.01))
  print(find_peaks(out_add, map, threshold=summary(operm_additive)[,1], drop = 1.5))
  
  # Plots
  male_add_tmp <- qtl_plot(out_add,map,c("slateblue", "violetred","forestgreen","gold2","red2"))+
    ggtitle(label=family_id)+
    theme(legend.position = "top",
          title=element_text(size=16))+
    scale_color_brewer(palette = "Dark2")
  
  for(i in 1:length(summary(operm_additive)[1,])){
    male_add_tmp <- male_add_tmp + geom_hline(yintercept=summary(operm_additive)[1,i],linetype="dashed",colour=brewer.pal(length(summary(operm_additive)[1,]),"Dark2")[i])
  }
  
  
  return(list(out_add,operm_additive,male_add_tmp,family_rear_covar))
})

# Bring together all the family plots
family_qtl <- lapply(family_figs,function(x){return(x[[3]])})
family_legend <- get_legend(family_qtl[[1]])
qtl_family_plots <- plot_grid(family_legend,
                              family_qtl[[1]]+theme(axis.text.x = element_blank(),
                                                    axis.title.x = element_blank(),
                                                    legend.position="none"),
                              family_qtl[[4]]+theme(axis.text.x = element_blank(),
                                                    axis.title.x = element_blank(),
                                                    legend.position="none"),
                              family_qtl[[3]]+theme(axis.text.x = element_blank(),
                                                    axis.title.x = element_blank(),
                                                    legend.position="none"),
                              family_qtl[[2]]+theme(
                                legend.position="none"),
                              ncol=1,align = "v",axis = "tblr",rel_heights = c(1,4,4,4,4))
pdf("figs/FigureSX_single_family_qtl_scans_male.pdf",width=12,height=10)
qtl_family_plots
dev.off()

##### Male cross-specific
# Get cross
YHF6QM7_cross <- (read_cross2("data/qtl_cross_YHF6QM7_males_informative.yaml"))
colnames(YHF6QM7_cross$pheno)[1:2] <- c("Days to maturity","Mature length")

# Get models 
YHF6QM7_mod <- family_figs[[which(families=="YHF6QM7")]][[1]]
YHF6QM7_perm <- family_figs[[which(families=="YHF6QM7")]][[2]]
YHF6QM7_covar <- family_figs[[which(families=="YHF6QM7")]][[4]]
YHF6QM7_map <- insert_pseudomarkers(YHF6QM7_cross$gmap, step=1)

# Get cutoffs
summary(YHF6QM7_perm, alpha=c(0.1,0.05,0.01))
YHF6QM7_peaks <- find_peaks(YHF6QM7_mod, YHF6QM7_map, threshold=min(summary(YHF6QM7_perm)), drop = 1.5)


# Get plots
YHF6QM7_size <- qtl_effects_plot(cross_in=YHF6QM7_cross,
                                 phenotype = "Mature length",
                                 qtl_chr=YHF6QM7_peaks$chr,
                                 qtl_pos=YHF6QM7_peaks$pos,
                                 qtl_low=YHF6QM7_peaks$ci_lo,
                                 qtl_high=YHF6QM7_peaks$ci_hi,
                                 qtl_covar = YHF6QM7_covar)
pdf("figs/FigureSX_male_size_qtl_effect_YHF6QM7.pdf",width=6,height=6)
plot_grid(plotlist = YHF6QM7_size,
          ncol=1,axis = "tblr",align="v")
dev.off()

# Get PVE
# Interbrood
calc_pve(LOD = 3.462264,nrow(YHF6QM7_covar))
# Get posiitons given transformations
max(YHF6QM7_map$chr23) - YHF6QM7_peaks$pos
max(YHF6QM7_map$chr23) - YHF6QM7_peaks$ci_lo
max(YHF6QM7_map$chr23) - YHF6QM7_peaks$ci_hi


