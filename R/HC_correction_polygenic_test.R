# HC Correction for regression of h2c and chrom size
lib <- c("data.table","qtl2","ggplot2")
sapply(lib,library,character.only=T)

###### FUNCTIONS #######
########### HC_correction ##########
# size; numeric vector of chromosome sizes 
# h2; numeric vector of chromosome heritabilities
# SE; numeric vector of standard errors for h2 estimates
# c; precision for adaptive permutation, used to calculate b
# alpha; false positive rate
# m; number of tests, used to calculate b
# r; minimum number of replicates 
# b; maximum number of replicatse default = round((1-alpha/m)/((c^2)*alpha/m))
# censored_to; value to which h2 estimates are censored to
# n; number of replciates. If NA adaptive permutation will be used, other wise n number of replicates will be used (appropriate when only a few tests need to be conducted).

HC_correction <- function(size, h2, SE, c=0.2, alpha=0.05, m=1, r=100,b=round((1-alpha/m)/((c^2)*alpha/m)), censored_to=0.000001, n=NA){
  data <- data.table(cbind(size, h2, SE))
  
  data <- data[apply(data, 1, function(x) !any(is.na(x))),] # remove rows with missing data
  
  obs <- cor.test(data$size, data$h2, alternative = 'greater')$p.value ## observed uncorrected p-value, one tailed OLS regression
  
  
  if(is.na(n)){
    ## run at lest r number of replicates
    exp <- sapply(1:r, function(x){
      data_re <- data.table(size=data[,size],h2=as.vector(sapply(data[,SE], function(x) rnorm(1, mean=0,sd=x))), SE=data[,SE])
      data_re[h2<censored_to, h2:=censored_to]
      cor.test(data_re$size, data_re$h2, alternative = 'greater')$p.value
    })
    
    counts <- table(exp<obs)
    
    if(length(counts)==1){
      counts <- c(counts, 0)
    }
    
    ## run until more than r more extreme values or b total number of replicates
    while(counts[2]<r & sum(counts)<b){
      data_re <- data.table(size=data[,size],h2=as.vector(sapply(data[,SE], function(x) rnorm(1, mean=0,sd=x))), SE=data[,SE])
      data_re[h2<censored_to, h2:=censored_to]
      exp <- cor.test(data_re$size, data_re$h2, alternative = 'greater')$p.value
      
      ifelse(exp<obs, counts[2] <- counts[2]+1, counts[1] <- counts[1]+1)
    }
    
    
    if(counts[2]==r){
      p <- r/sum(counts)
    }
    
    if(sum(counts)==b){
      p <- (counts[2]+1)/(b+1)
    }
    
    if(counts[2]!=r & sum(counts)!=b){
      p <- 'error'
    }
  }
  else{
    
    ## n number of replicates
    counts <- c(0,0)
    for(x in 1:n){
      data_re <- data.table(size=data[,size],h2=as.vector(sapply(data[,SE], function(x) rnorm(1, mean=0,sd=x))), SE=data[,SE])
      data_re[h2<censored_to,h2:=censored_to]
      exp <- cor.test(data_re$size, data_re$h2, alternative = 'greater')$p.value
      
      ifelse(exp<obs, counts[2] <- counts[2]+1, counts[1] <- counts[1]+1)
      
    }
    
    p <- (counts[2]+1)/(n+1)
    
  }
  
  return(c(p=obs,p_cor=p))
}

###### Correct GCTA results ########
# For each phenotype we need to read in per chromosome results, correct, and regress

# Get the new chromosome sizes from the merged scaffold vcf
scaf_merged_vcf <- read.vcfR("data/qtl_crosses_FINAL_inds.sorted.filtered.scafs_placed_with_chroms.vcf.gz")
chr_sizes <- sapply(1:23,function(x){
  max(as.integer(scaf_merged_vcf@fix[scaf_merged_vcf@fix[,1] == paste0("chr",x),2]))
})

# Get male/female qtl info
males <- read_cross2("data/qtl_cross_males_informative.yaml")
females <- read_cross2("data/qtl_cross_females_informative.yaml")

# # Fetch all the phenotypes
pheno_n <- 1:7

# Now HC-Correct each phenotype results
hc_corrected_res <- data.frame(rbindlist(lapply(pheno_n,function(x){
  set.seed(1234)
  # Set sex
  sex <- ifelse(x < 6,"female","male")
  x <- ifelse(sex == "male",x-5,x)
  
  # Read in each chromosome result and harvest info
  chr_res <- lapply(1:23,function(chr){
    
    res <- read.table(paste0("outputs/GCTA_scafs_merged_to_chrom_with_litters/",sex,"_qtl_crosses_heritability_pheno_",x,"_chr",chr,"_RES.hsq"),fill=T,header=T)
    
    se <- res$SE[4]
    var <- res$Variance[4]
    return(list(var,se))
  })
  
  # Vectorise
  all_var <- sapply(chr_res,function(x){return(x[[1]])})
  all_se <- sapply(chr_res,function(x){return(x[[2]])})
  
  # Fetch regression coefficients
  cor_res <- cor.test(all_var,chr_sizes,method = "pearson")
  
  # Perform HC corrected regression
  HC_res <- HC_correction(size = chr_sizes,
                          h2 = all_var,
                          SE = all_se)
  out <- matrix(ncol=3,nrow=1)
  out[1,1] <- cor_res$estimate
  out[1,2:3] <- HC_res
  colnames(out) <- c("r","p","p.corr")
  return(data.frame(out))
})))

# Add phenotype labels
pheno_vec <- c("Age at first brood (F)",
               "Size at first brood (F)",
               "First brood size (F)",
               "Interbrood period (F)",
               "First brood offspring weight (F)",
               "Age at maturity (M)",
               "Size at maturity (M)")
rownames(hc_corrected_res) <- pheno_vec
colnames(hc_corrected_res) <- c("Pearson's r","p","p (HC-corrected)")

# Tidy up
for(i in 1:3){
  hc_corrected_res[,i] <- round( hc_corrected_res[,i],3)
}

# Save table
write.csv(hc_corrected_res,
          "tables/TableX_polygenic_chrom_size_correlations.csv",
          quote = F)

# Read if we don't have
hc_corrected_res <- read.csv("tables/TableX_polygenic_chrom_size_correlations.csv")

#### Make correlations figure #####
plot_dd <- data.frame(rbindlist(lapply(pheno_n,function(x){
  # Set sex
  sex <- ifelse(x < 6,"female","male")
  y <- ifelse(sex == "male",x-5,x)
  
  # Read in each chromosome result and harvest info
  chr_res <- lapply(1:23,function(chr){
    
    res <- read.table(paste0("outputs/GCTA_scafs_merged_to_chrom_with_litters/",sex,"_qtl_crosses_heritability_pheno_",y,"_chr",chr,"_RES.hsq"),fill=T,header=T)
    
    se <- res$SE[4]
    var <- res$Variance[4]
    return(list(var,se))
  })
  
  # Vectorise
  all_var <- sapply(chr_res,function(x){return(x[[1]])})
  all_se <- sapply(chr_res,function(x){return(x[[2]])})
  
  # Do correlation
  cor_res <- cor.test(all_var,chr_sizes,method = "pearson")
  
  # Output 
  out <- data.frame(variance=all_var,
                    chr_size=chr_sizes,
                    pheno=pheno_vec[x])
  return(out)
})))

# Plot
cor_labs <- hc_corrected_res
cor_labs$label_txt <- paste0("r = ",hc_corrected_res$Pearson.s.r," p (HC) = ",hc_corrected_res$p..HC.corrected.)
cor_labs$Phenotype_F <- factor(cor_labs$X,levels=pheno_vec)
plot_dd$Phenotype_F <- factor(plot_dd$pheno,levels=pheno_vec)

hc_corrected_fig <- ggplot(plot_dd,aes(x=chr_size,y=variance))+
  facet_wrap(~Phenotype_F,ncol=3)+
  geom_point()+
  geom_smooth(method="lm",se=F,colour="red2")+
  labs(x="Chromosome Size (Mb)",y=expression(V[G]/V[P]))+
  theme_minimal()+
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16))+
  scale_x_continuous(breaks=c(1e7,2e7,3e7,4e7,5e7),
                     labels=c(10,20,30,40,50))+
  geom_label(data=cor_labs,aes(x=2e7,y=0.2,label=label_txt),hjust=0,vjust=1)

# Read back in chr partitioning fig
per_chr_heritability <- readRDS("figs/FigureX_perchr_heritability_GCTA_scafs_merged_to_chroms.rds")

# Plot on its own
pdf("figs/FigureSX_HC_corr_polygenic.pdf",width=10,height=10)
hc_corrected_fig
dev.off()

