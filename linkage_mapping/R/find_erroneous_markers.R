#####################################################################
# Read in a file of LODscores and remove erroneous markers
lib<-c("pracma")
lapply(lib,library,character.only=T)

args<-commandArgs(TRUE)
file=as.character(args[1])
output=as.character(args[2])

# Get file
#file<-"data/qtl_cross_full_pedigree_order_LG1_LODscores.txt"
dd<-(read.table(file,header=T))
snps<-data.frame((dd[,2:ncol(dd)]))
rownames(snps)<-1:nrow(snps)

# Normalise each row
for(i in 1:nrow(snps)){
  snps[i,]<-(snps[i,]-min(snps[i,]))/(max(snps[i,]-min(snps[i,])))
}

# Logical vector for whether max LOD is less than 1 SD from the mean
to_keep<-unlist(sapply(1:nrow(snps),function(x){
  #print(x)
  vec<-as.numeric(snps[x,])
  
  # First find peaks
  peaks<-findpeaks(vec,minpeakdistance = length(vec)/4, minpeakheight = 0.80)
  
  if(is.null(peaks)){
    return(FALSE)
  } else {
  
  # Return logical
   if(max(vec) < mean(vec) + sd(vec) | nrow(peaks) > 1){
     return(FALSE)
   } else {
     return(TRUE)
   }
  }
}))

# Output the markers
out_dd<-data.frame(X1=(1:nrow(snps))[!(to_keep)])

# Write
write.table(out_dd,
            output,
            col.names = F,row.names = F,quote = F,sep = "\t")




