library('phangorn')
library("ggtree")
library('ggnewscale')
library('viridis')
library('ggplot2')

##############################
# BV1 vs clinical G2 strains #
##############################

# Subset the BV1 strains
mlsa_tree <- read.tree('~/Documents/PostDoc_OSU/Projects/Avitis_genomics/data/fastTree_mlsa/avitis_mlsa.fasttree')
mlsa_tree <- midpoint(mlsa_tree)

# Tree
t0 <- ggtree(mlsa_tree, layout = "circular") + geom_tiplab()

t0

png(filename='results/big_fast_tree.png', width = 12000, height = 12000)
t0
dev.off()

# Update Biovar
getMRCA(mlsa_tree, c('CG971', 'W4'))
bv3_strain <- get_taxa_name(t0, node=4181)

getMRCA(mlsa_tree, c('Agrobacterium_tumefaciens_Yub001', 'Agrobacterium_tumefaciens_1D1609'))
bv1_strain <- get_taxa_name(t0, node=4314)

# BV1 strains
length(bv1_strain)

# Load scoary table
gene_pa <- read.table('/Users/arafat/Documents/PostDoc_OSU/Projects/Woods_clinical_agro/results/scoary/BV1_vs_pusense/PIRATE.gene_families.ordered.originalIDs.genePA.tsv', header=T)

bv1_pa <-  gene_pa[, which(bv1_strain %in% colnames(gene_pa))]

dim(bv1_pa)
colnames(bv1_pa)

bv1_strain # Pretty much contains all A pusense strains. If names does not match, need to update the dots "." to hiphen "-"


# Now get the list of A pusense clinical strains
apus_strains <- read.csv('results/scoary/BV1_vs_pusense/traits.clinical_vs_plants.csv')
apus_strains$X

dim(gene_pa[, which(apus_strains$X %in% colnames(gene_pa))])


# Check which A pus strains does not match with the gene_pa matrix
apus_strains$X %in% colnames(gene_pa)
# Swedish strains with "Inst2_Iso1_2014" is missing
# hmm ... making a new PIRATE run is a better decision ...

