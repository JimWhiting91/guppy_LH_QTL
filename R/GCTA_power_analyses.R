# Get power analysis inputs
library(data.table)

# Get GRM
males_grm <- read.table("~/Exeter/tmp/male_qtl_crosses_heritability_GRM.grm")
colnames(males_grm) <- c("id1","id2","shared_snps","gr")
females_grm <- read.table("~/Exeter/tmp/female_qtl_crosses_heritability_GRM.grm")
colnames(females_grm) <- c("id1","id2","shared_snps","gr")

# Get inds IDS and add cross family...
male_ids <- read.table("~/Exeter/tmp/male_qtl_crosses_heritability_GRM.grm.id")
colnames(male_ids) <- c("family1","vcf_id")
male_ids$id1 <- 1:nrow(male_ids)
male_ids$id2 <- 1:nrow(male_ids)
male_ids$family2 <- male_ids$family1

female_ids <- read.table("~/Exeter/tmp/female_qtl_crosses_heritability_GRM.grm.id")
colnames(female_ids) <- c("family1","vcf_id")
female_ids$id1 <- 1:nrow(female_ids)
female_ids$id2 <- 1:nrow(female_ids)
female_ids$family2 <- female_ids$family1

# Attach family to each grm
males_grm_merge <- merge(males_grm,male_ids[,c("id1","family1")],by="id1")
males_grm_merge <- merge(males_grm_merge,male_ids[,c("id2","family2")],by="id2")
males_grm_merge$unrelated <- "Yes"
males_grm_merge[males_grm_merge$family1 == males_grm_merge$family2,"unrelated"] <- "No"

females_grm_merge <- merge(females_grm,female_ids[,c("id1","family1")],by="id1")
females_grm_merge <- merge(females_grm_merge,female_ids[,c("id2","family2")],by="id2")
females_grm_merge$unrelated <- "Yes"
females_grm_merge[females_grm_merge$family1 == females_grm_merge$family2,"unrelated"] <- "No"

# Calculate off-diagonal variance
library(ggplot2)
ggplot(males_grm_merge,aes(x=gr,fill=unrelated))+geom_histogram()

male_unrelated_var <- var(males_grm_merge[males_grm_merge$family1 != males_grm_merge$family2,"gr"])
hist(males_grm_merge[males_grm_merge$family1 != males_grm_merge$family2,"gr"])
male_unrelated_var
female_unrelated_var <- var(females_grm_merge[females_grm_merge$family1 != females_grm_merge$family2,"gr"])
hist(females_grm_merge[females_grm_merge$family1 != females_grm_merge$family2,"gr"])
female_unrelated_var

# Build data.frames of male and female power -  values from GCTA power analysis shiny app, assuming alpha=0.05, variance of genetic relationships=male_var/female_var
male_power <- data.frame(h2=c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15,0.2,0.25,0.3,0.4,0.5),
                         power=c(0.0534,0.0638,0.0814,0.1064,0.1391,0.1795,0.2276,0.2826,0.3435,0.4089,0.7370,0.9330,0.9909,0.9994,1,1),
                         sex="Male (n=370)")

male_downsample_power <- data.frame(h2=c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15,0.2,0.25,0.3,0.4,0.5),
                                    power=c(0.0518,0.0572,0.0662,0.0790,0.0957,0.1163,0.1409,0.1697,0.2024,0.2389,0.4650,0.7040,0.8769,0.9628,0.9988,1),
                                    sex="Male (n=267)")

male_downsample200_power <- data.frame(h2=c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15,0.2,0.25,0.3,0.4,0.5),
                                       power=c(0.051,0.054,0.0591,0.0662,0.0754,0.0867,0.1003,0.1161,0.1342,0.1545,0.2889,0.464,0.6469,0.8008,0.9624,0.9967),
                                       sex="Male (n=200)")

female_power <- data.frame(h2=c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15,0.2,0.25,0.3,0.4,0.5),
                           power=c(0.0509,0.0537,0.0584,0.0651,0.0737,0.0843,0.0969,0.1116,0.1284,0.1474,0.2729,0.4389,0.6170,0.7731,0.9507,0.9947),
                           sex="Female (n=267)")

# Merge them together and plot
all_power <- rbind(male_power,male_downsample_power,male_downsample200_power,female_power)
gcta_power_plot <- ggplot(all_power,aes(x=h2,y=power,colour=sex))+
  geom_point()+
  geom_line()+
  theme_minimal()+
  labs(x=expression(h^2),y="Statistical Power",colour="Sex")+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=15),
        legend.text = element_text(size=14),
        legend.title = element_text(size=15))+
  scale_color_brewer(palette = "Dark2")


# Fetch power for downsampled males ---------------------------------------
male_downsample <- 100
nperm=100
male_downsample_unrelated_var <- sapply(1:nperm,function(iter){
  set.seed(iter)
  
  # Fetch males to keep
  males_to_keep <- male_ids[sample(1:nrow(male_ids),male_downsample),]
  
  # Downsample the GRM
  males_grm_merge_downsample <- males_grm_merge[males_grm_merge$id1 %in% males_to_keep$id1 &
                                                  males_grm_merge$id2 %in% males_to_keep$id1,]
  
  # Calculate unrelated variance and return...
  var(males_grm_merge_downsample[males_grm_merge_downsample$family1 != males_grm_merge_downsample$family2,"gr"])
})
mean(male_downsample_unrelated_var) # It's the same...


# Bring together with QTL Power -------------------------------------------

# Fetch male and female power matrices
male_power_matrix <- read.csv("data/male_qtl_power_matrix.csv")
colnames(male_power_matrix)[2:ncol(male_power_matrix)] <- seq(0.01,0.1,0.01)
rownames(male_power_matrix) <- male_power_matrix[,1]
male_power_matrix <- as.matrix(male_power_matrix[,2:ncol(male_power_matrix)])
male_power_qtl <- reshape2::melt(male_power_matrix)
colnames(male_power_qtl) <- c("Power","h2","sample_size")
male_power_qtl$sex <- "Male"

female_power_matrix <- read.csv("data/female_qtl_power_matrix.csv")
colnames(female_power_matrix)[2:ncol(female_power_matrix)] <- seq(0.01,0.1,0.01)
rownames(female_power_matrix) <- female_power_matrix[,1]
female_power_matrix <- as.matrix(female_power_matrix[,2:ncol(female_power_matrix)])
female_power_qtl <- reshape2::melt(female_power_matrix)
colnames(female_power_qtl) <- c("Power","h2","sample_size")
female_power_qtl$sex <- "Female"

# Bind together and plot
both_power_qtl <- rbind(male_power_qtl,female_power_qtl)
both_power_qtl$Power_F <- factor(both_power_qtl$Power,levels=sort(unique(both_power_qtl$Power)))
both_power_qtl$h2_F <- factor(both_power_qtl$h2,levels=sort(unique(both_power_qtl$h2)))

# Mark the cells to annotate
both_power_qtl$cell_annotate <- ""
both_power_qtl[both_power_qtl$sex=="Male" & both_power_qtl$sample_size <= 370,"cell_annotate"] <- "*"
both_power_qtl[both_power_qtl$sex=="Female" & both_power_qtl$sample_size <= 267,"cell_annotate"] <- "*"

power_matrix_plot <- ggplot(both_power_qtl,aes(x=Power_F,y=h2_F,fill=sample_size))+
  geom_tile()+
  facet_wrap(~sex,ncol=1)+
  scale_fill_viridis_c(trans="log10")+
  labs(x="Statistical Power",y=expression(h^2),fill="Required\nSample\nSize")+
  geom_text(aes(label=cell_annotate),colour="white",size=5)+
  theme_bw()+
  theme(axis.title=element_text(size=16),
        axis.text = element_text(size=14),
        axis.text.x = element_text(size=14,angle=45,hjust=1),
        strip.text = element_text(size=16))


# Plot QTL powers
qtl_powers <-read.csv("data/female_male_qtl_power.csv")
qtl_powers_plot <- ggplot(qtl_powers,aes(x=h2,y=Power,colour=Sex))+
  geom_point()+
  geom_line()+
  theme_minimal()+
  labs(x=expression(h^2),y="Statistical Power",colour="Sex")+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=15),
        legend.text = element_text(size=14),
        legend.title = element_text(size=15))+
  scale_color_brewer(palette = "Dark2")

# Plot them all together
pdf("figs/FigureSX_Power_analyses_summary.pdf",width=12,height=8)
plot_grid(
  plot_grid(gcta_power_plot,qtl_powers_plot,ncol=1,axis = "tblr",align="v",labels="AUTO",label_size=32,hjust = -0.1),
  power_matrix_plot,labels=c("","C"),label_size=32,rel_widths=c(1.5,1),
  ncol=2)
dev.off()

