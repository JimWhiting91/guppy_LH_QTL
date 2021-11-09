#########################################################################
# Script to visualise LepMAP linkage maps vs genetic distance
# Here we're visualising our v10 map

lib<-c("ggplot2","data.table")
lapply(lib,library,character.only=T)

# Read in the data for each LG and plot
LG_plots<-lapply(1:23,function(x){

  # Get the linkage group group
  cM<-read.table(paste0("data/v10/qtl_cross_full_pedigree_v10_inform123_order_LG",x,"_BEST.SA.txt"))[,1:2]
  mapped<-read.table(paste0("data/v10/qtl_cross_full_pedigree_v10_inform123_order_LG",x,"_BEST.SA.mapped.clean.txt"))[,1:2]
  colnames(cM)<-c("marker","cM")
  colnames(mapped)<-c("chr","bp")
  
  # Correct mapped
  mapped$bp<-gsub(pattern = "_",replacement = "",mapped$bp)
  mapped$bp<-gsub("[*].*$","",mapped$bp)
  
  # merge
  cM$mapped<-mapped$bp
  cM$original_chr<-mapped$chr
  
  # Correct
  cM$mapped<-as.numeric(cM$mapped)
  
  # Visualise
  ggplot(cM,aes(x=mapped,y=cM))+
    geom_point()+
    facet_wrap(~original_chr)
  
  # And filter for the main chr
  table(cM$original_chr)
  main_chr<-names(table(cM$original_chr))[c(table(cM$original_chr))==max(table(cM$original_chr))]
  out<-ggplot(cM,aes(x=mapped,y=cM,colour=original_chr))+
    geom_point()+
    labs(x=paste0(main_chr," pos"),
         y=paste0(main_chr," cM"))
  
  cM$LG<-paste0("LG",x)
  return(list(out,cM))
}
)


# Plot all to pdf
pdf("figs/preliminary_lepmap_STAR_comparisons_v10.pdf",width=10,height=10)
for(i in 1:23){
  print(LG_plots[[i]][[1]])
}
dev.off()

# Get the rest
cM_dd<-data.frame(rbindlist(lapply(LG_plots,function(chr){return(chr[[2]])})))

# Make a supp table describing position of scaffolds on "original_chr"
chr_scaf_supp <- data.frame(rbindlist(lapply(1:23,function(chr){
  
  # Get LG
  LG <- cM_dd[cM_dd$original_chr == paste0("chr",chr),"LG"]
  LG <- names(which.max(table(LG)))
  
  # Get all scafs on that LG
  scafs_tmp <- cM_dd[cM_dd$LG == LG,]
  scafs_tmp <- scafs_tmp[grep("chr",scafs_tmp$original_chr,invert = T),]
  
  scafs_out <- unique(scafs_tmp$original_chr)
  out_mat <- matrix(ncol=5,nrow=length(scafs_out))
  
  # Fill mat
  if(length(scafs_out) > 0){
  for(i in 1:length(scafs_out)){
    out_mat[i,1] <- paste0("chr",chr)
    out_mat[i,2] <- scafs_out[i]
    out_mat[i,3] <- min(scafs_tmp[scafs_tmp$original_chr==scafs_out[i],"cM"])
    out_mat[i,4] <- max(scafs_tmp[scafs_tmp$original_chr==scafs_out[i],"cM"])
    out_mat[i,5] <- ifelse(as.numeric(out_mat[i,3]) < abs(as.numeric(out_mat[i,3])-max(cM_dd[cM_dd$LG == LG,"cM"])),"Start","End")
  }
  }
  
  # Save as adata.frame
  out <- as.data.frame(out_mat)
  
  # Invert start end if needs be
  is_inverted <- cor.test(cM_dd[cM_dd$LG == LG,"mapped"],cM_dd[cM_dd$LG == LG,"cM"])$estimate
  if(is_inverted < 0){
    out[out$V5 == "End","V5"] <- "tmp"
    out[out$V5 == "Start","V5"] <- "End"
    out[out$V5 == "tmp","V5"] <- "Start"
  }
  
  # Return
  return(out)
})))

# Tidy and output
colnames(chr_scaf_supp) <- c("Chromosome","Scaffold","Min cM","Max cM","Relative position in assembly")
chr_scaf_supp$Scaffold <- gsub("F","F_",chr_scaf_supp$Scaffold)
write.table(chr_scaf_supp,
            "tables/TableSX_linkage_map_chromosome_position_of_unplaced_scaffolds.txt",
            row.names = F,quote = F,sep = "\t")
##########################################################
# Trimming...
# Plot all to pick trims

plots<-lapply(1:23,function(x){
  ggplot(cM_dd[cM_dd$LG == paste0("LG",x),],aes(x=mapped,y=cM,colour=original_chr))+
    geom_point()
})

# Go through and define LG Trims
LG=c(1:23)
lower=c(0,15,15,14,12,17,17,17,0,17,32,11,0,1,22,31,1,27,10,1,7,5,21)
upper=c(101,89,89,89,82,95,85,87,70,90,102,88,75,77,95,110,90,91,83,70,84,75,94)

# and check em
plots2<-lapply(1:23,function(x){
  ggplot(cM_dd[cM_dd$LG == paste0("LG",x),],aes(x=mapped,y=cM,colour=original_chr))+
    geom_point()+
    geom_hline(yintercept = c(lower[x],upper[x]))
})

# Now earmark markers to trim
trimmed_cM<-data.frame(rbindlist(lapply(1:23,function(x){
  tmp<-cM_dd[cM_dd$LG == paste0("LG",x),]
  
  # Assign trimming values
  tmp$trim<-"No"
  tmp[tmp$cM < lower[x] |
        tmp$cM > upper[x],"trim"]<-"Yes"
  return(tmp)
})))

# Write this output
write.table(trimmed_cM,
            "data/v10/trimming_guide.txt",
            row.names = F,quote = F,sep = "\t")

# Also do the subtracting...
subtracted<-data.frame(rbindlist(lapply(1:23,function(x){
  tmp<-trimmed_cM[trimmed_cM$LG == paste0("LG",x),]
  tmp2<-tmp[tmp$trim == "No",]
  
  # Subtract
  tmp2$cM<-tmp2$cM-min(tmp2$cM)
  
  return(tmp2)
})))

write.table(subtracted,
            "data/v10/trimmed_map.txt",
            row.names = F,quote = F,sep = "\t")

# And finally plot final maps
final_plots<-lapply(1:23,function(x){
  # Subset
  tmp<-subtracted[subtracted$LG == paste0("LG",x),]
  
  # Get main_chr
  main_chr<-names(table(tmp$original_chr))[c(table(tmp$original_chr))==max(table(tmp$original_chr))]
  
  # Plot
  ggplot(tmp,aes(x=mapped,y=cM,colour=original_chr))+
    geom_point()+
    xlab("Chr/Scaf Pos")+
    ggtitle(main_chr)
})

# Plot
pdf("~/Exeter/qtl_crosses/figs/lepmap_STAR_comparisons_v10_trimmed.pdf",width=10,height=10)
for(i in 1:23){
  print(final_plots[[i]])
}
dev.off()

#################################################
# Final map figures...
plot_dd <- data.frame(rbindlist(lapply(1:23,function(x){
  
  # Subset
  tmp<-subtracted[subtracted$LG == paste0("LG",x),]
  
  # Get main_chr
  main_chr<-names(table(tmp$original_chr))[c(table(tmp$original_chr))==max(table(tmp$original_chr))]
  tmp$main_chr <- main_chr
  tmp$chr_scaf <- tmp$original_chr
  tmp[grep("chr",tmp$chr_scaf,invert = T),"chr_scaf"] <- "scaf"
  tmp[grep("chr",tmp$chr_scaf),"chr_scaf"] <- "chr"
  
  #Return
  return(tmp)
})))

# Tidy up
plot_dd$chr_F <- factor(plot_dd$main_chr,levels = paste0("chr",1:23))

# Plot
genetic_map_fig <- ggplot(plot_dd,aes(mapped,cM,colour=chr_scaf))+
  geom_point(alpha=0.5,size=1)+
  facet_wrap(~chr_F,scales="free")+
  scale_colour_manual(breaks=c("chr","scaf"),
                      values = c("black","red2"))+
  theme_minimal()+
  theme(axis.text = element_blank(),
        axis.title = element_text(size=16),
        strip.text = element_text(size=16),
        legend.position = "right",
        legend.text = element_text(size=16))+
  labs(x="Male Guppy Genome Position (bp)",y="Genetic Map Position (cM)",colour="")

pdf("figs/FigureSX_genetic_map.pdf",width=10,height=10)
genetic_map_fig
dev.off()



