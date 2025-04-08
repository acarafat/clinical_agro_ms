library(tidyr)
library(dplyr)

###############################
# Old sourmash search results #
###############################
sourmash_comp_matrix <- read.csv("results/plasmid/containment/contain_plasmid.csv")
# Label the rows
rownames(sourmash_comp_matrix) <- colnames(sourmash_comp_matrix)
# Transform for plotting
sourmash_comp_matrix <- as.data.frame(sourmash_comp_matrix)


strainlist <- read.table("results/plasmid/containment/strainlist",header=F, stringsAsFactors=F) # all whole genome, complete and draft
plasmidlist <-read.table("results/plasmid/containment/plasmidlist",header=F, stringsAsFactors=F) # Known plasmid queries from complete genomes

# check if strainlist matches, otherwise update the mismatch
for (i in strainlist$V1){
  if (! (i %in% rownames(sourmash_comp_matrix)) ){
    print(i)
  }
}

# check if plasmidlist matches, otherwise update the mismatch
for (i in plasmidlist$V1){
  if (! (i %in% rownames(sourmash_comp_matrix)) ){
    print(i)
  }
}

filtered_mat <- sourmash_comp_matrix[strainlist$V1, plasmidlist$V1]
filtered_mat$strain <- row.names(filtered_mat)

tab_long <- gather(filtered_mat, plasmid, similarity,-strain)
tab_long_filt<-tab_long[tab_long$similarity >= 0.90,]

###############################
# New sourmash search results #
###############################
tab_long_filt <- read.csv('results/plasmid/containment/whole_genome_containment/contain_plasmid_genome.csv', header=T,stringsAsFactors=F)
colnames(tab_long_filt) <- c('plasmid', 'strain', 'similarity')

plasmid_types<-read.table("results/plasmid/containment/plasmid_types.mcl",header=T,stringsAsFactors=F)

tab_long_filt$plasmid_type <- plasmid_types[match(tab_long_filt$plasmid,plasmid_types$plasmid),]$plasmid_type
tab_long_filt$inference <- 2

head(tab_long_filt)

uniq_tab <- tab_long_filt[,c("strain","plasmid_type", "inference")] %>% distinct() %>% arrange(strain,plasmid_type)

write.table(uniq_tab,file="results/plasmid/containment/contain_plasmid.mcl.csv.types.long",row.names=F,quote=F,sep="\t")


###############
# PIVOT Shape #
###############
df <- read.table("results/plasmid/containment/mydataset_k21.csv.cluster.0.2.list", header=T)
df$strain <- gsub("__.*$", "", df$strain)



df_reformatted <- df %>%
  mutate(value = 1) %>%
  reshape2::dcast(strain~plasmid_type, value.var="value", fill=0, fun.aggregate=length) %>%
  mutate(across(where(is.numeric), ~ ifelse(. != 0, 1, 0))) # %>%
  #pivot_wider(names_from = plasmid_type, values_from = value, values_fill = list(value = 0)) 

dim(df_reformatted)

write.table(df_reformatted, "results/plasmid/containment/original_type_table.mcl.tab", sep = "\t", row.names = FALSE, quote = FALSE)


################## 
# Merge inferred #
##################
original_tab <- read.table("results/plasmid/containment/original_type_table.mcl.tab",header=T,stringsAsFactors=F, check.names = FALSE)

inferred_tab <- read.table("results/plasmid/containment/contain_plasmid.mcl.csv.types.long",header=T,stringsAsFactors=F)

# Process sourmash4 inferrence
#inferred_tab <- read.table("results/plasmid/containment/contain_plasmid_sourmash4.csv", header=T, sep=",")
#inferred_tab$inference <- 2
#inferred_tab <- inferred_tab[,c("match", "plasmid_type", "inference")]
#colnames(inferred_tab) <- c("strain", "plasmid_type", "inference")

orig_long <- gather(original_tab, plasmid_type, inference,-strain)


orig_long_filt <- orig_long[orig_long$inference > 0,]

orig_long_filt$inference <- "detected"
inferred_tab$inference <- "inferred"

merged_tab<- rbind(orig_long_filt, inferred_tab) %>% distinct(strain, plasmid_type, .keep_all = TRUE)

wide_tab <- spread(merged_tab, plasmid_type, inference, fill="Not Found")


dim(wide_tab)
head(wide_tab)

write.table(wide_tab,file="results/plasmid/containment/merged_types_table.mcl.tab",row.names=F,quote=F,sep="\t")
