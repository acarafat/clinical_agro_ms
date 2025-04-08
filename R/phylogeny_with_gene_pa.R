library(ggtree)
library(ggplot2)
library(phangorn)
library(data.table)
library(RColorBrewer)
library(ggthemes)
library(viridis)
library(dplyr)
library(treeio)
library(ggnewscale)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

tol21rainbow= c("#771155", "#AA4488", "#CC99BB", 
                "#114477", "#4477AA", "#77AADD", 
                "#117777", "#44AAAA", "#77CCCC", 
                "#117744", "#44AA77", "#88CCAA", 
                "#777711", "#AAAA44", "#DDDD77", 
                "#774411", "#AA7744", "#DDAA77", 
                "#771122", "#AA4455", "#DD7788")


# sampletypecolors <- tableau_color_pal(palette="Tableau 10", type="regular",direction=1)(10) # Tableau 10 color


sampletypecolors <- c(
  "hospital environment" = "#777711",
  "terrestrial" = "#114477",
  "animal" = "#DDDD70",
  "built environment" = "#77AADD",
  "plant" = "#44AA77",
  "clinical" = "#771122",
  "other" = "grey30",
  "soil" = "#CC99BB",
  "unknown" = "grey"
)

sampletypeshapes <- c(
  "hospital environment" = 21,
  "terrestrial" = 22,
  "animal" = 23,
  "built environment" = 12,
  "plant" = 25,
  "clinical" = 8,
  "other" = 9,
  "soil" = 24,
  "unknown" = 7
)

################
# Plot on tree #
################
snptree <- read.iqtree('results/tree/species1_calls_merged.vcf.nonpassnonhetfilt.vcf.vcffilter.vcf.tab.noref.fasta.treefile')
snptree@phylo <- midpoint(snptree@phylo)
meta <- read.table('results/tree/metadata.table.tsv', header=T, sep='\t')
meta$treelab <- paste(meta$strain, meta$short_location,sep="   ")

#original root from midpointing gets set as NA by phangorn midpoint
snptree@phylo$node.label[is.na(snptree@phylo$node.label)] <- "100/100"
#tip labels, only for plotting bootstrap as color
snptree@phylo$node.label[snptree@phylo$node.label == ""] <- "100/100"

# Add scale-bar
# Add Clade annotation

t <- ggtree(snptree) %<+% meta[, c(2:ncol(meta))] + 
  geom_tiplab(aes(label = treelab, x = x + 0.02), size=4) +
  #geom_tiplab2(aes(label=short_location), hjust = 0.5, angle = 0, size=2.5) +
  geom_tippoint(aes(x = x + 0.01, fill=sample_type, color=sample_type, shape=sample_type), size=2.5) +
  scale_fill_manual(values=sampletypecolors, name='Sample Type') +
  scale_color_manual(values=sampletypecolors, name='Sample Type') +
  scale_shape_manual(values=sampletypeshapes, name='Sample Type') +
  #labs(color = 'Sample source') +
  new_scale_fill() +
  geom_treescale(x= 0.02, y=50) +
  theme(legend.position = c(0.15, 0.80)) 

t$data$bootstrap <- '0'

t$data$SH_aLRT[is.na(t$data$SH_aLRT)] <- 100
t$data$UFboot[is.na(t$data$UFboot)] <- 100

t$data[which(t$data$SH_aLRT >= 80 & t$data$UFboot  >= 95),]$bootstrap <- '1'


t <- t + new_scale_color() +
  geom_tree(aes(color=bootstrap == '1')) +
  scale_color_manual(name='Bootstrap', values=setNames( c('black', 'grey'), c(T,F)), 
                     labels=c("Unsupported","UFBoot > 95% & SH-aLRT > 80%"))



t


####################
# Load PIRATE data #
####################

# Create presence-absence matrix for all contrast
gene_fam <- read.csv('results/PIRATE.gene_families.ordered.csv', header=T)
gene_fam.pw <- gene_fam[, c(1, 23:ncol(gene_fam))]

###################################################################################
# Parse PIRATE gene families list and subset the OLG to plot along with phylogeny #
###################################################################################


####################
# Clinical vs rest #
####################
clinical_vs_rest <- read.csv('results/scoary/clinical_vs_rest_10_06_2024_1743.results.csv', header=T)
clinical_vs_rest <- clinical_vs_rest[which(clinical_vs_rest$Bonferroni_p < 0.05),]
clinical_vs_rest$Gene <- gsub("__.*", "", clinical_vs_rest$Gene)


# Let's make a function for ploting scoary output with a tree
plot_scoary_with_tree <- function (tree, meta, pirate_gene_fam_pam, scoary_gene_list){
  # Subset gene_fam.pw for gene PAV matrix
  # clinical_vs_rest.pav > pav_mat
  # t.clinical_vs_rest.pav > pav_mat.t
  # gene_fam.pw > pirate_gene_fam_pam
  # clinical_vs_rest > scoary_gene_list
  # meta.clinical_vs_rest > meta.pam
  
  # tree <- t
  # meta <- meta
  # pirate_gene_fam_pam <- gene_fam.pw
  # scoary_gene_list <- combined
  
  pav_mat <- pirate_gene_fam_pam[which(pirate_gene_fam_pam$allele_name %in% scoary_gene_list$Gene),]
  rownames(pav_mat) <- pav_mat$allele_name
  
  pav_mat <- pav_mat[, c(2:ncol(pav_mat))]
  pav_mat[pav_mat != ''] <- "Present"
  pav_mat[pav_mat ==''] <- "Absent"
  
  # Transpose
  pav_mat.t <- as.data.frame(t(pav_mat))
  colnames(pav_mat.t) <- rownames(pav_mat)
  rownames(pav_mat.t) <- colnames(pav_mat)
  
  # Now need to match the PIRATE genome name with tree tip names
  pav_mat.t$genome <- rownames(pav_mat.t)
  
  meta.pam <- merge(meta, pav_mat.t, by.x='genome', by.y='genome')
  strainlist <- meta.pam$strain
  meta.pam <- meta.pam[, c(12:ncol(meta.pam)),]
  rownames(meta.pam) <- strainlist
  
  meta.pam <- meta.pam[scoary_gene_list$Gene]
  
  gheatmap(tree, meta.pam, offset = 0.2, width = 0.8,
           colnames_angle = 60, colnames_position = 'top', colnames_offset_y = 0.2, 
           custom_column_labels = scoary_gene_list$Non.unique.Gene.name,hjust = 0,
           font.size = 3) + vexpand(0.12) +
    scale_fill_manual(values=c('Present'='steelblue', 'Absent'='grey90'), 
                      labels=c('Present'='Present', 'Absent'='Absent'), name='Gene Presence/Absence') #+
    #labs(fill = "Gene PAV")
}



png("results/clinical_vs_rest.png", width = 1700, height = 1700)
plot_scoary_with_tree(t, meta, gene_fam.pw, clinical_vs_rest) 
dev.off()

###########################
# Curated contrast genes  #
###########################
combined <- read.csv("results/scoary/combined_contrast.csv", header=T)
combined$Gene <- gsub("__.*", "", combined$Gene)

pdf("Manuscript/Figures/Main/Figure_1._gene_PAV_combined_curated_contrast.pdf", width = 20, height = 16)
plot_scoary_with_tree(t, meta, gene_fam.pw, combined)
dev.off()


#####################
#####################
clinical_vs_plant <- read.csv('results/scoary/clinical_vs_plant_12_06_2024_1348.results.csv', header=T)
clinical_vs_plant <- clinical_vs_plant[which(clinical_vs_plant$Bonferroni_p < 0.05),]
clinical_vs_plant$Gene <- gsub("__.*", "", clinical_vs_plant$Gene)

dim(clinical_vs_plant)


png("results/clinical_vs_plant.png", width = 3000, height = 1700)
plot_scoary_with_tree(t, meta, gene_fam.pw, clinical_vs_plant) 
dev.off()

# Plant vs rest
plant_vs_rest <- read.csv('results/scoary/plant_vs_rest_10_06_2024_1743.results.csv', header=T)
plant_vs_rest <- plant_vs_rest[which(plant_vs_rest$Bonferroni_p < 0.05),]
plant_vs_rest$Gene <- gsub("__.*", "", plant_vs_rest$Gene)

dim(plant_vs_rest)

plot_scoary_with_tree(t, meta, gene_fam.pw, plant_vs_rest)
