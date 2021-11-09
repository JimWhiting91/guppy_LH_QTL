# Compare QTL candidate genes to Sharma transcriptome
lib <- c("ggplot2","data.table","biomaRt","seqinr")
lapply(lib,library,character.only=T)

GENOME="~/Exeter/Genomes/Poecilia_reticulata.Guppy_female_1.0_MT.dna.toplevel.fa"
TRANSCRIPTOME="~/Exeter/Genomes/trin_cuff_v14_cdhit90.fa"

# Set up Ensembl
ensembl <- useMart("ensembl")
guppy <- useDataset("preticulata_gene_ensembl",mart=ensembl)
listAttributes(guppy)
biomart_attributes <- c("ensembl_gene_id","chromosome_name","start_position","end_position","cdna")

# Fetch candidates
candidates <- read.csv("tables/TableSX_qtl_candidates_annotated_v2.csv")

# Get candidates for each qtl
qtls <- unique(candidates$QTL)
qtl_candidates <- lapply(qtls,function(x){
  return(unique(candidates[candidates$QTL == x,]))
})

# Get the Sharma sex-biased transcripts
sharma_res <- read.csv("data/sharma_sex_biased_transcripts.csv")

# Run over QTLs
sex_biased_qtls <- lapply(c(1,3),function(i){
  
# Map back to the transcriptome one-by-one
chr19_qtl <- qtl_candidates[[i]]

transcriptome_hits <- data.frame(rbindlist(lapply(1:nrow(chr19_qtl),function(x){
  
  # Fetch cdna
  cdna <- getBM(attributes = biomart_attributes,
                filters = "ensembl_gene_id",
                values = chr19_qtl$Ensembl.Gene[x],
                mart=guppy)
  
  # Write to a fasta
  fasta_temp <- tempfile(pattern = "qtl_cdna", fileext = '.fa')
  paf_temp <- tempfile(pattern = "qtl_cdna", fileext = '.paf')
  
  on.exit({ unlink(fasta_temp) })
  on.exit({ unlink(paf_temp) })

  # Write to temp fasta
  seqinr::write.fasta(as.list(cdna$cdna),
                      names=paste0(cdna$ensembl_gene_id,"_",1:nrow(cdna)),
                      file.out=fasta_temp)
  
  # Do the alignments
  system(paste0("~/bin/minimap2 -t6 ",TRANSCRIPTOME," ",fasta_temp," > ",paf_temp),wait=T)
  
  # Read in and filter the PAF
  paf <- read.table(paf_temp,fill=T)
  
  # Only keep high quality alignments
  paf <- paf[paf$V12 == 60,]
  paf <- paf[order(paf$V1,paf$V3),]
  
  # Only keep aligments that are >20% length of transcript
  paf <- paf[paf$V11 > 0.2*paf$V7 | paf$V4 - paf$V3 > 0.2*paf$V2,]
  
  if(nrow(paf) > 0){
    paf$gene <- chr19_qtl$Ensembl.Gene[x]
    
    # Add the percentage of coverage
    paf$transcript_cover <- round((paf$V11/paf$V7)*100,2)
    paf$cdna_cover <- round(((paf$V4-paf$V3)/paf$V2)*100,2)
    
    # Add the percentage of coverage
    paf[paf$V11 > 0.2*paf$V7 | paf$V4 - paf$V3 > 0.2*paf$V2,]
    
    return(paf)
  } else {
    return(NULL)
  }
  
})))

# Filter for small transcripts
transcriptome_hits <- transcriptome_hits[transcriptome_hits$V7 > 200,]
transcriptome_hits <- transcriptome_hits[,c(1:12,19,20,21)]

# Compare with the sharma res...
female_biased_genes <- unique(transcriptome_hits[transcriptome_hits$V6 %in% sharma_res[sharma_res$logFC < 0,"Gene"],"gene"])
male_biased_genes <- unique(transcriptome_hits[transcriptome_hits$V6 %in% sharma_res[sharma_res$logFC > 0,"Gene"],"gene"])

# Add these markers back into supp table
chr19_qtl$sex_biased_expression_gonads <- ""
chr19_qtl[chr19_qtl$Ensembl.Gene %in% female_biased_genes,"sex_biased_expression_gonads"] <- "Female-biased"
chr19_qtl[chr19_qtl$Ensembl.Gene %in% male_biased_genes,"sex_biased_expression_gonads"] <- "Male-biased"


return(chr19_qtl)
})

# Fill the second one as well just because
tmp_qtl <- qtl_candidates[[2]]
tmp_qtl$sex_biased_expression_gonads <- ""

# Bind them
out_table <- rbind(sex_biased_qtls[[1]],tmp_qtl,sex_biased_qtls[[2]])

# Write
write.table(out_table,
          "tables/TableSX_qtl_candidates_annotated_v2_with_transcriptome_results.tsv",
          row.names = F,quote = F,sep="\t")

# Read in and get results
out_table <- read.table("tables/TableSX_qtl_candidates_annotated_v2_with_transcriptome_results.tsv",fill=T,sep="\t",header=T)

# Get counts of female-biased / male-biased for each qtl
unique(out_table[out_table$QTL == "chr19_offspring_weight" & out_table$sex_biased_expression_gonads == "Female-biased","Ensembl.Gene"  ])
unique(out_table[out_table$QTL == "chr19_offspring_weight" & out_table$sex_biased_expression_gonads == "Male-biased","Ensembl.Gene"  ])
unique(out_table[out_table$QTL == "chr12_female_interbrood" & out_table$sex_biased_expression_gonads == "Female-biased","Ensembl.Gene"  ])
unique(out_table[out_table$QTL == "chr12_female_interbrood" & out_table$sex_biased_expression_gonads == "Male-biased","Ensembl.Gene"  ])

# # Compare...
# transcriptome_hits[transcriptome_hits$gene %in% genes_with_sex_bias,]
# chr19_qtl[chr19_qtl$Ensembl.Gene %in% genes_with_sex_bias,]
# 
# # Examine closest gene
# transcriptome_hits[transcriptome_hits$gene == "ENSPREG00000018824",]
# sharma_res[sharma_res$Gene %in% transcriptome_hits[transcriptome_hits$gene == "ENSPREG00000018824","V6"],]
#            