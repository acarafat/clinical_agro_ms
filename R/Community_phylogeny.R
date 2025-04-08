################################
# Community Phylogeny Analysis #
################################

library(ape)
library(picante)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(ggpubr)

`%nin%` <- Negate(`%in%`)

snp.tree <- read.tree('results/tree/species1_calls_merged.vcf.nonpassnonhetfilt.vcf.vcffilter.vcf.tab.noref.fasta.treefile')
mlsa.tree <- read.tree('results/combined_mlsa/G2.partition.rxml.treefile')

snp.meta <- read.table('results/tree/metadata.table.tsv', header=T, sep='\t')
mlsa.meta <- read.csv('results/combined_mlsa/combined_meta.csv', header=T)

snp.meta$sample_type2 <- snp.meta$sample_type
snp.meta[which(snp.meta$sample_type %nin% c('clinical', 'plant')),]$sample_type2 <- 'other_combined'


mlsa.meta$sample_type2 <- mlsa.meta$sample_type
mlsa.meta[which(mlsa.meta$sample_type %nin% c('clinical', 'plant')),]$sample_type2 <- 'other_combined'

################################### 
# Community analysis for MLSA tree #
###################################

#############################################
# Convert meta variable to community format #
#############################################

mlsa.comm <- t(as.matrix(table(mlsa.meta$fullstrain, mlsa.meta$sample_type)))

mlsa.clean.tree <- match.phylo.comm(mlsa.tree, mlsa.comm)$phy
mlsa.clean.comm <- match.phylo.comm(mlsa.tree, mlsa.comm)$comm

# Phylogenetic alpha diversity (Faith's PD)
site.pd <- pd(samp= mlsa.clean.comm, 
              tree= mlsa.clean.tree,
              include.root = F)

cor.test(site.pd$PD, site.pd$SR)

plot(site.pd$PD, site.pd$SR, xlab = "Phylogenetic Diversity", ylab = "Species Richness", pch = 16)


site.pd.df <- as.data.frame(site.pd)
site.pd.df$site <- rownames(site.pd)
head(site.pd.df)

# Test standardized effect size of PD of each community

site.ses.pd <- ses.pd(samp=mlsa.clean.comm,
                      tree=mlsa.clean.tree,
                      null.model="taxa.labels",
                      include.root = F,
                      runs=1000)



site.pd.df <- cbind(site.pd.df, site.ses.pd)

p1 <- ggplot(site.pd.df, aes(x=PD, y=SR, label=site)) + 
  geom_point(aes(color = ifelse(PD <= 0.05, 'red', 'blue')), show.legend = F) +
  xlab("Phylogenetic Diversity") + ylab("Strain Richness") +
  geom_text_repel(hjust=0, vjust=0) +
  ggtitle("A. Alpha diversity in all sample categories\n   Correlation between PD & SR = 0.77*") +
  theme_classic2()

##########
mlsa.comm <- t(as.matrix(table(mlsa.meta$fullstrain, mlsa.meta$sample_type2)))

mlsa.clean.tree <- match.phylo.comm(mlsa.tree, mlsa.comm)$phy
mlsa.clean.comm <- match.phylo.comm(mlsa.tree, mlsa.comm)$comm

# Phylogenetic alpha diversity (Faith's PD)
site.pd <- pd(samp= mlsa.clean.comm, 
              tree= mlsa.clean.tree,
              include.root = F)

cor.test(site.pd$PD, site.pd$SR)


site.pd.df <- as.data.frame(site.pd)
site.pd.df$site <- rownames(site.pd)


# Test standardized effect size of PD of each community

site.ses.pd <- ses.pd(samp=mlsa.clean.comm,
                      tree=mlsa.clean.tree,
                      null.model="taxa.labels",
                      include.root = F,
                      runs=1000)



site.pd.df <- cbind(site.pd.df, site.ses.pd)

p2 <- ggplot(site.pd.df, aes(x=PD, y=SR, label=site)) + 
  geom_point(aes(color = ifelse(PD <= 0.05, 'red', 'blue')), show.legend = F) +
  xlab("Phylogenetic Diversity") + ylab("Strain Richness") +
  geom_text_repel(hjust=0, vjust=0) +
  ggtitle("B. Alpha diversity in clinical, plant, and all other samples combined\n    Correlation between PD & SR = -0.31") +
  theme_classic2()



ggarrange(p1, p2)
######################################
# Phylogenetic community relatedness #
######################################
site.cophenDist <- cophenetic.phylo(mlsa.clean.tree)

# NRI Net Relatedness Index
site.ses.mpd <- ses.mpd(mlsa.clean.comm, site.cophenDist, null.model="taxa.labels", abundance.weighted = F, runs=100)
site.nri <- as.matrix(-1 * ((site.ses.mpd[,2] - site.ses.mpd[,3]) / site.ses.mpd[,4]))
rownames(site.nri) <- row.names(site.ses.mpd)
colnames(site.nri) <- "NRI"
site.nri #positive means more phylogenetic clustering

nri.nti.summary <- cbind(site.nri, site.ses.mpd[,c('mpd.obs.p')])

# NTI Net Taxon Index
site.ses.mntd <- ses.mntd(mlsa.clean.comm, site.cophenDist, null.model = 'taxa.labels', abundance.weighted = F, runs=100)
site.nti <- as.matrix(-1*((site.ses.mntd[,2] - site.ses.mntd[,3]) / site.ses.mntd[,4]))

rownames(site.nti) <- rownames(site.ses.mntd)
colnames(site.nti) <- 'NTI'
site.nti #positive means more phylogenetic clustering

nri.nti.summary <- cbind(nri.nti.summary, site.nti)
nri.nti.summary <- cbind(nri.nti.summary, site.ses.mntd[, 'mntd.obs.p'])
colnames(nri.nti.summary) = c('NRI', 'NRI p-value', 'NTI', 'NTI p-value')
View(nri.nti.summary)


################################### 
# Community analysis for SNP tree #
###################################

#############################################
# Convert meta variable to community format #
#############################################

snp.comm <- t(as.matrix(table(snp.meta$strain, snp.meta$sample_type)))

snp.clean.tree <- match.phylo.comm(snp.tree, snp.comm)$phy
snp.clean.comm <- match.phylo.comm(snp.tree, snp.comm)$comm

# Phylogenetic alpha diversity (Faith's PD)
site.pd <- pd(samp= snp.clean.comm, 
              tree= snp.clean.tree,
              include.root = F)

cor.test(site.pd$PD, site.pd$SR)

site.pd.df <- as.data.frame(site.pd)
site.pd.df$site <- rownames(site.pd)


# Test standardized effect size of PD of each community
site.ses.pd <- ses.pd(samp=snp.clean.comm,
                      tree=snp.clean.tree,
                      null.model="taxa.labels",
                      include.root = F,
                      runs=1000)



site.pd.df <- cbind(site.pd.df, site.ses.pd)


p1 <- ggplot(site.pd.df, aes(x=PD, y=SR, label=site)) + 
  geom_point(aes(color = ifelse(PD <= 0.05, 'red', 'blue')), show.legend = F) +
  xlab("Phylogenetic Diversity") + ylab("Strain Richness") +
  geom_text_repel(hjust=0, vjust=0) +
  ggtitle("A. Alpha diversity in all sample categories\n   Correlation between PD & SR = 0.85*") +
  theme_classic2()

snp.comm <- t(as.matrix(table(snp.meta$strain, snp.meta$sample_type2)))

snp.clean.tree <- match.phylo.comm(snp.tree, snp.comm)$phy
snp.clean.comm <- match.phylo.comm(snp.tree, snp.comm)$comm

# Phylogenetic alpha diversity (Faith's PD)
site.pd <- pd(samp= snp.clean.comm, 
              tree= snp.clean.tree,
              include.root = F)

cor.test(site.pd$PD, site.pd$SR)


site.pd.df <- as.data.frame(site.pd)
site.pd.df$site <- rownames(site.pd)


# Test standardized effect size of PD of each community
site.ses.pd <- ses.pd(samp=snp.clean.comm,
                      tree=snp.clean.tree,
                      null.model="taxa.labels",
                      include.root = F,
                      runs=1000)



site.pd.df <- cbind(site.pd.df, site.ses.pd)


p2 <- ggplot(site.pd.df, aes(x=PD, y=SR, label=site)) + 
  geom_point(aes(color = ifelse(PD <= 0.05, 'red', 'blue')), show.legend = F) +
  xlab("Phylogenetic Diversity") + ylab("Strain Richness") +
  geom_text_repel(hjust=0, vjust=0) +
  ggtitle("B. Alpha diversity in clinical, plant, and all other samples combined\n   Correlation between PD & SR = 0.67") +
  theme_classic2()

ggarrange(p1, p2)

######################################
# Phylogenetic community relatedness #
######################################
site.cophenDist <- cophenetic.phylo(snp.clean.tree)

# NRI Net Relatedness Index
site.ses.mpd <- ses.mpd(snp.clean.comm, site.cophenDist, null.model="taxa.labels", abundance.weighted = F, runs=100)
site.nri <- as.matrix(-1 * ((site.ses.mpd[,2] - site.ses.mpd[,3]) / site.ses.mpd[,4]))
rownames(site.nri) <- row.names(site.ses.mpd)
colnames(site.nri) <- "NRI"
site.nri #positive means more phylogenetic clustering

nri.nti.summary <- cbind(site.nri, site.ses.mpd[,c('mpd.obs.p')])

# NTI Net Taxon Index
site.ses.mntd <- ses.mntd(snp.clean.comm, site.cophenDist, null.model = 'taxa.labels', abundance.weighted = F, runs=100)
site.nti <- as.matrix(-1*((site.ses.mntd[,2] - site.ses.mntd[,3]) / site.ses.mntd[,4]))

rownames(site.nti) <- rownames(site.ses.mntd)
colnames(site.nti) <- 'NTI'
site.nti #positive means more phylogenetic clustering

nri.nti.summary <- cbind(nri.nti.summary, site.nti)
nri.nti.summary <- cbind(nri.nti.summary, site.ses.mntd[, 'mntd.obs.p'])
colnames(nri.nti.summary) <- c('NRI', 'NRI p-value', 'NTI', 'NTI p-value')
View(nri.nti.summary)

