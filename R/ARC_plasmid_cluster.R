############################
# Visualize AutoMLSA2 Tree #
############################

library('phangorn')
library("ggtree")
library('viridis')
library('ggnewscale')
library('ggplot2')
library('ggtree')
library('stringr')
library('tidyverse')
library('reshape2')

# Load AutoMLSA2 tree
mlsa <- midpoint(read.tree('results/tree/bv1_tree/BV1_mlsa.nex.treefile'))

droplist <- c("1D1410", "22-3674b3", "22-3674a2", "22-3674b1", "22-3574a1", "22-3674a1", "22-3674b4", 
              "22-3574b1", "22-3638c",  "22-3638B1", "22-3639A1", "22-3639C1", "22-3639B1", "22-3639C2", 
              "22-3639B2", "22-3639A2", "22-3638B2", "22-3638A2", "22-797-1", "22-798-1",  "22-692-1a", 
              "22-692-2",  "22-210-1",  "22-211-1",  "22-221-1",  "22-223-1",  "22-214-1",  "22-224-1", 
              "22-222-1",  "22-226-1",  "22-209-1", "1D1409", "1D1408", "1D1420", "1D1407", "1D1419", 
              "1D1342", "1D166",  "1D1418", "1D1410", "1D1341", "1D1424", "1D1589", "1D1340", "1D1434", 
              "1D1339", "1D1465", "1D1436", "1D1411", "1D1119", "C6.1A.A", "CG117", "CG974",  "CG653", "CG628",
              "CG1057", "CG49", "CG1056", "CG987", "CG160-95", "CG492", "CG490", "CG491", "P30_93", "H",
              "H53_94", "H3.2.1", "CJ94_95", "A3-95", "O80_94", "8924", "8996", "8628", "55", "229",
              "16-2014-1-2a", "16-172Ci", "AJW0106", "AJW0128", "16_6", "0d", "U1_95", "2022-3a", "N15-94",
              "N53-94", "315", "b.bak", "NSC", "8924", "371", "387", "IL21", "IL20", "6_12", "IVE-57_1",
              "IVE-48_1", "IVE-89_4", "IVE-190_1", "IVE-80_2", "47-2", "CFBP_2407", "CFBP_2642", 
              "CFBP_2732", "CFBP_2682", "47-3", "16-2014-1-2a", "16_6", "16-172Ci", "48-55", "48-63", "48-52", 
              "48-31", "585", "ICMP_19468", "ARF003_13.5.A", "I121_95", "21-3430-2c")

mlsa <- drop.tip(mlsa, droplist)

 
#plasmids <- read.table('results/mydataset_k21.csv.cluster.0.1.list', header=T)
plasmids <- read.table('results/tree/bv1_tree/mydataset_k21.csv.cluster.0.2.list', header=T)


# Update strain names
plasmids$strain <- plasmids$plasmid
plasmids$strain <- gsub("_filtered_plasmid.*", "", plasmids$strain)
plasmids$strain <- gsub("_pTi.*", "", plasmids$strain)
plasmids$strain <- gsub("_NZ_.*", "", plasmids$strain)
plasmids$strain <- gsub("_plasmid.*", "", plasmids$strain)
plasmids$strain <- gsub("_pAt.*", "", plasmids$strain)
plasmids$strain <- gsub("__.*", "", plasmids$strain)
plasmids$strain <- gsub("\\.plasmid.*", "", plasmids$strain)
plasmids$strain <- gsub("_LUKX0100000.*", "", plasmids$strain)
plasmids$strain <- gsub("_LT00973(.*)$", "", plasmids$strain)
plasmids$strain <- gsub("_HG518323(.*)$", "", plasmids$strain)
plasmids$strain <- gsub("_CAICSX020000002.1", "", plasmids$strain)



# Check if strain name from plasmid df does not match with tree tiplabs
t <- ggtree(mlsa)

plasmids$updated_strain <- gsub("\\.", "-", plasmids$strain)

# Update strain id
plasmids[which(plasmids$updated_strain == "Agrobacterium_leguminum_ST04-17-025"),]$updated_strain <- "Agrobacterium_leguminum_ST04.17.025"
plasmids[which(plasmids$updated_strain == "Agrobacterium_leguminum_ST07-17-026"),]$updated_strain <- "Agrobacterium_leguminum_ST07.17.026"
plasmids[which(plasmids$updated_strain == "Agrobacterium_pusense_ST15-13-056"),]$updated_strain <- "Agrobacterium_pusense_ST15.13.056"
plasmids[which(plasmids$updated_strain == "Agrobacterium_rhizogenes_B_4-1"),]$updated_strain <- "Agrobacterium_rhizogenes_B_4.1"
plasmids[which(plasmids$updated_strain == "Agrobacterium_rhizogenes_ST04-16-045"),]$updated_strain <- "Agrobacterium_rhizogenes_ST04.16.045"
plasmids[which(plasmids$updated_strain == "Agrobacterium_rhizogenes_ST07-17-004"),]$updated_strain <- "Agrobacterium_rhizogenes_ST07.17.004"
plasmids[which(plasmids$updated_strain == "Agrobacterium_rhizogenes_ST07-17-018"),]$updated_strain <- "Agrobacterium_rhizogenes_ST07.17.018"
plasmids[which(plasmids$updated_strain == "Agrobacterium_rhizogenes_ST07-17-029"),]$updated_strain <- "Agrobacterium_rhizogenes_ST07.17.029"
plasmids[which(plasmids$updated_strain == "Agrobacterium_rhizogenes_ST15-13-057"),]$updated_strain <- "Agrobacterium_rhizogenes_ST15.13.057"
plasmids[which(plasmids$updated_strain == "Agrobacterium_rhizogenes_ST15-16-020"),]$updated_strain <- "Agrobacterium_rhizogenes_ST15.16.020"
plasmids[which(plasmids$updated_strain == "Agrobacterium_rhizogenes_ST15-16-024"),]$updated_strain <- "Agrobacterium_rhizogenes_ST15.16.024"
plasmids[which(plasmids$updated_strain == "Agrobacterium_salinitolerans_ST04-16-212"),]$updated_strain <- "Agrobacterium_salinitolerans_ST04.16.212"
plasmids[which(plasmids$updated_strain == "Agrobacterium_salinitolerans_ST07-17-032"),]$updated_strain <- "Agrobacterium_salinitolerans_ST07.17.032"
plasmids[which(plasmids$updated_strain == "Agrobacterium_salinitolerans_ST15-13-006"),]$updated_strain <- "Agrobacterium_salinitolerans_ST15.13.006"
plasmids[which(plasmids$updated_strain == "Agrobacterium_salinitolerans_ST15-13-013"),]$updated_strain <- "Agrobacterium_salinitolerans_ST15-13-013"
plasmids[which(plasmids$updated_strain == "Agrobacterium_salinitolerans_ST15-13-091"),]$updated_strain <- "Agrobacterium_salinitolerans_ST15.13.091"
plasmids[which(plasmids$updated_strain == "Agrobacterium_salinitolerans_ST15-13-097"),]$updated_strain <- "Agrobacterium_salinitolerans_ST15.13.097"
plasmids[which(plasmids$updated_strain == "Agrobacterium_salinitolerans_ST15-16-021"),]$updated_strain <- "Agrobacterium_salinitolerans_ST15.16.021"
plasmids[which(plasmids$updated_strain == "Agrobacterium_salinitolerans_ST15-16-055"),]$updated_strain <- "Agrobacterium_salinitolerans_ST15.16.055"
plasmids[which(plasmids$updated_strain == "Agrobacterium_salinitolerans_ST15-13-013"),]$updated_strain <- "Agrobacterium_salinitolerans_ST15.13.013"
plasmids[which(plasmids$updated_strain == "Agrobacterium_sp_10MFCol1-1"),]$updated_strain <- "Agrobacterium_sp_10MFCol1.1"
plasmids[which(plasmids$updated_strain == "Agrobacterium_sp_224MFTsu3-1"),]$updated_strain <- "Agrobacterium_sp_224MFTsu3.1"
plasmids[which(plasmids$updated_strain == "Agrobacterium_sp_224MFTsu3-1"),]$updated_strain <- "Agrobacterium_sp_224MFTsu3.1"
plasmids[which(plasmids$updated_strain == "Agrobacterium_sp_33MFTa1-1_ASM432849v1"),]$updated_strain <- "Agrobacterium_sp_33MFTa1.1_ASM432849v1"
plasmids[which(plasmids$updated_strain == "Agrobacterium_sp_224MFTsu3-1"),]$updated_strain <- "Agrobacterium_sp_224MFTsu3.1"
plasmids[which(plasmids$updated_strain == "Agrobacterium_sp_SRR10754060_bin-19_metawrap_v1-3-0_MAG"),]$updated_strain <- "Agrobacterium_sp_SRR10754060_bin.19_metawrap_v1.3.0_MAG"
plasmids[which(plasmids$updated_strain == "Agrobacterium_sp_ST15-13-013"),]$updated_strain <- "Agrobacterium_sp_ST15.13.013"
plasmids[which(plasmids$updated_strain == "Agrobacterium_sp_ST15-13-015"),]$updated_strain <- "Agrobacterium_sp_ST15.13.015"
plasmids[which(plasmids$updated_strain == "Agrobacterium_sp_ST15-13-095"),]$updated_strain <- "Agrobacterium_sp_ST15.13.095"
plasmids[which(plasmids$updated_strain == "Agrobacterium_sp_ST15-16-024"),]$updated_strain <- "Agrobacterium_sp_ST15.16.024"
plasmids[which(plasmids$updated_strain == "Agrobacterium_sp_ST15-16-055"),]$updated_strain <- "Agrobacterium_sp_ST15.16.055"
plasmids[which(plasmids$updated_strain == "Agrobacterium_sp_SRR14536366_bin-3_metawrap_v1-3-0_MAG"),]$updated_strain <- "Agrobacterium_sp_SRR14536366_bin.3_metawrap_v1.3.0_MAG"
plasmids[which(plasmids$updated_strain == "Agrobacterium_tumefaciens_Ag_tumefaciens_EHA105_Agrobacterium_tumefaciens_Ag_tumefaciens_EHA105"),]$updated_strain <- "Agrobacterium_tumefaciens_Ag_tumefaciens_EHA105"
plasmids[which(plasmids$updated_strain == "Agrobacterium_tumefaciens_G_8-3"),]$updated_strain <- "Agrobacterium_tumefaciens_G_8.3"
plasmids[which(plasmids$updated_strain == "ARF003_13-5-A"),]$updated_strain <- "ARF003_13.5.A"
plasmids[which(plasmids$updated_strain == "C6-1A-A"),]$updated_strain <- "C6.1A.A"
plasmids[which(plasmids$updated_strain == "H3-2-1"),]$updated_strain <- "H3.2.1"
plasmids[which(plasmids$updated_strain == ""),]$updated_strain <- ""

# Update tree dataset
t$data[which(t$data$label == 'Agrobacterium_fabrum_GV3101__pMP90'),]$label <- 'Agrobacterium_fabrum_GV3101'
t$data[which(t$data$label == 'Agrobacterium_tumefaciens_M56__79'),]$label <- 'Agrobacterium_tumefaciens_M56'

# Check missing
`%ni%` <- Negate(`%in%`)

tiplabs <- t$data$label
missing <- c()

for (id in unique(plasmids$updated_strain)){
  if ( id %ni% tiplabs) {
    missing <- c(missing, id)
  }
}

missing


# Compare
plasmids[startsWith(plasmids$strain, 'ARF003_13'),]
as.data.frame(t$data[startsWith(t$data$label, 'Agrobacterium_tumefaciens_Ag_tumefaciens_EHA105'),])

t$data[which(t$data$label == 'Agrobacterium_fabrum_GV3101'),]$label

missing

# Update tips
plasmids[which(plasmids$strain == '2022.3a')]$strain <- "2022-3a"


## Characterize plasmids 
#plasmids <- plasmids %>% 
#  mutate(type = case_when(
#    grepl("pTi", plasmid) & grepl("pAt", plasmid, ignore.case = TRUE) ~ "pTi_pAt",
#    grepl("pTi", plasmid) ~ "pTi",
#    grepl("pAt", plasmid) ~ "pAt",
#    TRUE ~ NA_character_ # Assign NA if neither pTi nor pAt are present
#  ))

# Characterize plasmisd
pTi <- read.table("results/tree/bv1_tree/pTi_plasmid_file_match.tsv", sep='\t', header=F)
colnames(pTi) <- c('query', 'match')

pAt <- read.csv("results/tree/bv1_tree/pAt_containment.csv", header=T)

plasmids$type <- NA # 1 means present / uncharacterized plasmids

# Combine pTi
plasmids[which(plasmids$plasmid %in% pTi$match),]$type <- "pTi"

#for (i in pTi$match){
#  if (i %in% plasmids$plasmid){
#    plasmids[which(plasmids$plasmid == i),]$type <- "pTi"
#  } else {
#    print(i)
#  }
#}

pTi_clusters <- unique(plasmids[which(plasmids$type == "pTi"),]$cluster)
plasmids[which(plasmids$cluster %in% pTi_clusters),]$type <- "pTi"


# Combine pAt
head(pAt)

for (strain in pAt$genome){
  if (strain %in% plasmids$updated_strain) {
    cluster <- pAt[which(pAt$genome == strain),]$cluster[1]
    if (cluster %in% plasmids[which(plasmids$updated_strain == strain),]$cluster){
      plasmids[which(plasmids$updated_strain == strain & plasmids$cluster == cluster & is.na(plasmids$type)),]$type <- "pAt" 
    }
    else {
      plasmids <- rbind(plasmids, c(NA, cluster, strain, strain, "inferred_pAt"))
    }
  }
  else {
    print(strain)
  }
}


tail(plasmids)


## Plasmid_clusters
clusters <- plasmids[, c('updated_strain', 'cluster', 'type')]

head(clusters)

result <- clusters %>%
  mutate(value = case_when(
    is.na(type) ~ "1",   # If type is NA, assign "1"
    type == "pTi" ~ "pTi",
    type == "pAt" ~ "pAt",
    type == "inferred_pAt" ~ "inferred_pAt", 
    TRUE ~ "0"           # If no cluster, assign "0"
  )) %>%
  group_by(updated_strain, cluster) %>%  
  summarize(value = max(value), .groups = "drop") %>%
  dcast(updated_strain ~ cluster, value.var = "value", fill = "0") 



meta <- result[,2:length(colnames(result))]
meta <- meta[, order(as.numeric(names(meta)))]
rownames(meta) <- result$updated_strain

meta[] <- lapply(meta, as.character) # Ensure all values are character

# Print the final transformed meta data
#sum(as.numeric(meta$`171`))

meta <- meta[,1:170]

#result <- clusters %>%
#  mutate(value = 1) %>%
#  group_by(updated_strain, cluster) %>%  # Group by strain and cluster
#  summarize(value = max(value), .groups = "drop") %>%
#  reshape2::dcast(updated_strain ~ cluster, value.var = "value", fill = 0, fun.aggregate = length)
  
#meta <- result[,2:length(colnames(result))]
#rownames(meta) <- result$updated_strain
#meta[meta > 0] <- 1
#meta[] <- lapply(meta, as.character)


## plot plamid presence absence with tree
t2 <- gheatmap(t, meta, offset = 0.02,
               colnames_angle = 90, colnames_position = 'top', hjust = 0, color='white', 
               colnames_offset_y = 0.5, font.size = 3, legend_title = "Plasmid")  + 
  vexpand(0.04) + #new_scale_fill() 
  scale_fill_manual(values=c('1'='steelblue', '0'='grey90', 'pAt'='black', 'pTi'='red', 'inferred_pAt'='grey'), 
                    labels=c('1'='Present', '0'='Absent', 'pAt'='pAt', 'pTi'='pTi', 'inferred_pAt'='Inferred pAt'), na.value = 'white') +
  labs(fill = "Plasmid")


pdf("Manuscript/Figures/BV1_NCBI_plasmids.0.2.pdf", width = 50, height = 65)
t2 + geom_tiplab(size=4)
dev.off()






################
# Update meta  #
################
#meta <- read.csv('meta/Airtable_15Feb2024/Agrobacterium-All metadata.csv', header=T)
#meta <- meta[,c('strain_code', 'organism', 'oncogenic.plasmid', 'location', 'plant_host', 'plant_host_full', 'year')]
#meta <- meta[!duplicated(meta[ , c("strain_code")]),]

t <- ggtree(mlsa)

head(t$data$label)


####################

plasmids


match_string_ids <- function(df, col_name, id_vector, new_col_name) {
  new_col_name <- deparse(substitute(new_col_name))
  col_name <- deparse(substitute(col_name))
  
  df[[new_col_name]] <- NA_character_ # Initialize the new column with NA
  
  for (id in id_vector) {
    matches <- startsWith(df[[col_name]], id)
    print(df[matches,][[new_col_name]])
    df[matches,][[new_col_name]] <- id
  }
  return(df)
}

match_string_ids(plasmids, "plasmid", t$data$label, "strain")
