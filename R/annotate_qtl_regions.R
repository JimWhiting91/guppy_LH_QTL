# Get gene information for QTL regions...
lib<-c("qtl2","ggplot2","biomaRt","data.table")
lapply(lib,library,character.only=T)

######################################################################################################
# Function takes a set of genome regions, maps back to genome and returns genes with NCBI, Uniprot and KEGG IDs etc
# Fetch the Uniprot IDs...
ensembl_digger<-function(outliers=NULL,
                         genome1=NULL,
                         genome2=NULL,
                         aln_block_length=1000,
                         aln_qual=40,
                         biomart_attributes=c("ensembl_gene_id","entrezgene_id","uniprot_gn_id","kegg_enzyme","chromosome_name","start_position","end_position"),
                         n.cores=1,
                         mart=NULL){
  
  # Firstly, take all the outlier regions and make a multi-fasta
  lapply(outliers,function(x){
    system(paste0("samtools faidx ",genome1," ",x," >> outputs/tmp_multi.fa"))
  })
  
  # Align to old genome
  system2('/Users/jimwhiting/bin/minimap2',
          args=c(paste0("-t ",n.cores),genome2,"outputs/tmp_multi.fa"),
          stdout ="outputs/tmp_multi.paf",wait=T)
  
  # Fetch regions
  aln<-read.table("outputs/tmp_multi.paf",fill=T)
  
  # Keep tidy!
  system(paste0("rm -f outputs/tmp_multi*"))
  
  # Only keep "high support alignment regions"
  aln <- aln[aln$V11 > aln_block_length & aln$V12 >= aln_qual,]
  aln <- aln[order(aln$V3),]
  regions<-paste0(aln$V6,":",as.integer(aln$V8),":",as.integer(aln$V9))
  
  # Pull uniprot genes from biomaRt for each region
  tmp_biomart<-getBM(attributes = biomart_attributes,
                     filters= "chromosomal_region",
                     values=regions,
                     mart=mart)
  
  
  # Return
  return(tmp_biomart)
}
######################################################################################################
# Set up Ensembl
ensembl <- useMart("ensembl")
guppy <- useDataset("preticulata_gene_ensembl",mart=ensembl)
zebrafish <- useDataset("drerio_gene_ensembl",mart=ensembl)
danio_universe<-getBM(attributes = c("ensembl_gene_id","drerio_homolog_ensembl_gene","drerio_homolog_orthology_type"),
                      mart=guppy)
danio_universe_one2one <- danio_universe[danio_universe$drerio_homolog_orthology_type == "ortholog_one2one",]

filters <- listFilters(guppy)
guppy_attributes <- listAttributes(guppy)

# Set biomart attributes to pull
biomart_atts<-c("external_gene_name","ensembl_gene_id","chromosome_name","start_position","end_position","go_id","name_1006","kegg_enzyme")

# QTL 'outliers'
qtl_regions <- list("chr19:18602889-19602889",
                    c("000111F_0:528180-1199617","chr22:23415429-24223839"),
                    c("000149F_0:47124-197258","chr12:24705290-24525856"))
names(qtl_regions) <- c("chr19_offspring_weight","chr22_female_size","chr12_female_interbrood")

# Fetch Ensembl info: Gene names, GO, Zebrafish orthologue for pathway info...
qtl_candidates <- data.frame(rbindlist(lapply(names(qtl_regions),function(qtl){
  
  # Run digger
  digger_res <- ensembl_digger(outliers=qtl_regions[[qtl]],
                               genome1="/Users/jimwhiting/Exeter/Genomes/STAR.chromosomes.release.fasta",
                               genome2="/Users/jimwhiting/Exeter/Genomes/Poecilia_reticulata.Guppy_female_1.0_MT.dna.toplevel.fa",
                               aln_block_length = 5000,
                               aln_qual=60,
                               biomart_attributes = biomart_atts,
                               n.cores=6,
                               mart = guppy)
  
  # Fetch danio orthologies
  digger_res$danio_orthologues <- NA
  for(i in 1:nrow(digger_res)){
   # print(i)
     if(nrow(danio_universe_one2one[danio_universe_one2one$ensembl_gene_id == digger_res$ensembl_gene_id[i],]) > 0){
      digger_res$danio_orthologues[i] <- danio_universe_one2one[danio_universe_one2one$ensembl_gene_id == digger_res$ensembl_gene_id[i],"drerio_homolog_ensembl_gene"]
   }
  }
  
  digger_res$qtl <- qtl
  
  # Get gene identifiers
  genes <- paste0(digger_res$chromosome_name,"_",digger_res$start_position)
  print(paste0("QTL ",qtl," has ",length(unique(genes))," genes"))
    
return(digger_res)
})))

# Write to a file
write.csv(qtl_candidates,
            "tables/TableSX_qtl_candidates.csv",
            quote = F,row.names = F)

# Write separate danio orthologies...
danio_orth <- data.frame(orthologue = unique(na.omit(qtl_candidates$danio_orthologues)))
write.table(danio_orth,
            "outputs/qtl_danio_orthologues.txt",
            row.names = F,quote = F,sep="\t",col.names = F)
