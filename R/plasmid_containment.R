library(ggtree)
library(phangorn)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(viridis)
library(ggnewscale)
library(dplyr)

library(ggrepel)
library(ape)
library(picante)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


sampletypecolors <- tableau_color_pal(palette="Tableau 10", type="regular",direction=1)(10)
names(sampletypecolors) <- c("hospital environment","terrestrial","animal","built environment","plant","yellow","clinical","other","soil","unknown")


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

# Plasmid containment tree


################
# Plot on tree #
################
snptree <- midpoint(read.tree('results/tree/species1_calls_merged.vcf.nonpassnonhetfilt.vcf.vcffilter.vcf.tab.noref.fasta.treefile'))
meta <- read.table('results/tree/metadata.table.tsv', header=T, sep='\t')

plasmid <- read.table('results/plasmid/containment/merged_types_table.mcl.tab', header=T, sep='\t', check.names = FALSE)


meta <- merge(meta, plasmid, by.x= 'fullstrain', by.y='strain', all.x = T)

# Check everything in plasmid got merged
# If not found, manuallmy update the name from meta to plasmid
for (i in plasmid$strain) {
  if ( ! (i %in% meta$fullstrain)){
    print(i)
  }
}

t <- ggtree(snptree) %<+% meta[, c(2:ncol(meta))] + 
  geom_tiplab(aes(x = x + 0.015), size=4) +
  geom_tippoint(aes(x = x + 0.01, fill=sample_type, color=sample_type, shape=sample_type),size=4) +
  scale_fill_manual(values=sampletypecolors, name='Sample Type') +
  scale_color_manual(values=sampletypecolors, name='Sample Type') +
  scale_shape_manual(values=sampletypeshapes, name='Sample Type') +
  labs(color = 'Sample source') +
  theme(legend.position = c(0.15, 0.85)) +
  new_scale_fill() 

t

metaColLen <- length(colnames(meta))
meta.plasmid <- as.data.frame(meta[,12:metaColLen])
meta.plasmid[is.na(meta.plasmid)] <- "Not Found"
colnames(meta.plasmid) <- as.numeric(colnames(meta)[12:metaColLen])
rownames(meta.plasmid) <- meta$strain

meta.plasmid <- meta.plasmid[, c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                                 "11", "12", "13", "14", "15", "16", "17", "18", "19", '20',
                                 "21", "22", "23", "24", "25", "26", "27", "28", "29", "30",
                                 "31", "32", "33", "34", "35")]


t2 <- gheatmap(t, meta.plasmid, width=0.9, offset = 0.12, 
         colnames_angle = 90, colnames_position = 'top', hjust = 0, color='white', 
         colnames_offset_y = 0.5, font.size = 4, legend_title = "Plasmid")  + 
  vexpand(0.04) + #new_scale_fill() 
  scale_fill_manual(values=c('detected'='navy', 'inferred'='steelblue', 'Not Found'='grey90'), 
                    labels=c('detected'='Detected', 'inferred'='Inferred', 'Not Found'='NA')) +
  labs(fill = "Plasmid Category")


pdf("Manuscript/Figures/Main/Figure_3_plasmid_containment.pdf", width = 18, height = 18)
t2
dev.off()

# Plasmid comparison between clinical vs non-clinical strains
comparison <- meta.plasmid
comparison[comparison == "detected"] <- as.numeric(1)
comparison[comparison == "inferred"] <- as.numeric(1)
comparison[comparison == "Not Found"] <- as.numeric(0)

comparison$sum <- rowSums(sapply(comparison, function(x) as.numeric(as.character(x))))
comparison$strain <- rownames(comparison)
comparison <- merge(comparison, meta[c('strain','sample_type')], by='strain')
#comparison$sample_type[comparison$sample_type != 'clinical'] <- 'non_clinical'
#t.test(sum ~ sample_type, data = comparison)
head(comparison)

anova_results <- aov(sum ~ sample_type, data = comparison)
summary(anova_results)
# Run Tukey's HSD
post_hoc <- TukeyHSD(anova_results)
print(post_hoc)

################
# Need to find the association of PD between shared plasmids
plasmidDF <- meta.plasmid
plasmidDF[plasmidDF == 'Not Found'] = 'Absent'
plasmidDF[plasmidDF != 'Absent'] = 'Present'
plasmidDF$Strain <- rownames(plasmidDF)
head(plasmidDF)

plasmidDF <- reshape2::melt(plasmidDF, id='Strain')
colnames(plasmidDF) <- c('Strain', 'Plasmid', 'PAV')
plasmidDF <- plasmidDF[which(plasmidDF$PAV == 'Present'), ]

pComm <- t(table(plasmidDF$Strain, plasmidDF$Plasmid))

#
snp.clean.tree <- match.phylo.comm(snptree, pComm)$phy
snp.clean.comm <- match.phylo.comm(snptree, pComm)$comm

# Phylogenetic alpha diversity (Faith's PD)
site.pd <- pd(samp= snp.clean.comm, 
              tree= snp.clean.tree,
              include.root = F)

site.pd.df <- as.data.frame(site.pd)
site.pd.df$site <- rownames(site.pd)

cor.test(site.pd.df$PD, site.pd.df$SR)



ggplot(site.pd.df, aes(x=SR, y=PD, label=site)) + 
  geom_point(aes(size=SR), show.legend = F, alpha=0.3) +
  ylab("Phylogenetic Diversity") + xlab("Strain Richness") +
  geom_text_repel(hjust=0, vjust=0) +
  ggtitle("") 
  
