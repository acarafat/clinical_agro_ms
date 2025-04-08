library(ggtree)
library(ggplot2)
library(phangorn)
library(data.table)
library(RColorBrewer)
library(viridis)
library(ggthemes)
library(treeio)
library(ggnewscale)
library(dplyr)
library(tangler)
library(phytools)
library(ggnewscale)


#########################
# Prepare combined meta #
#########################

#meta1 <- read.table('results/combined_mlsa/metadata.table.tsv', header=T, sep='\t')
#meta1 <- meta1[, c('fullstrain', 'strain', 'sample_type')]
#colnames(meta1) <- c('fullstrain', 'strain', 'sample_type')

#meta2 <- read.table('results/combined_mlsa/metadata.txt', header=T, sep='\t')
#meta2 <- meta2[, c('Strain', 'strain', 'sample_type')]
#colnames(meta2) <- c('fullstrain', 'strain', 'sample_type')

#meta <- rbind(meta1, meta2)

#write.csv(meta, 'results/combined_mlsa/combined_meta.csv', row.names = F, quote = F)


#sampletypecolors <- tableau_color_pal(palette="Tableau 10", type="regular",direction=1)(12)
#names(sampletypecolors) <- c("hospital environment","terrestrial","animal",
#                             "built environment","plant",'nematode',"clinical","other","soil","undocumented", "NA")

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



mlsatree <- read.iqtree('results/combined_mlsa/partition.rxml.treefile')

# Drop LMG17935, LMG90T1
mlsatree <- drop.tip(mlsatree, tip = c('LMG17935', 'LMG90T1'))


keepnode <- getMRCA(mlsatree@phylo, c("AGR11","RTH4"))

#G2.clade <- extract.clade(mlsatree@phylo, keepnode)
#G2.clade <- as.treedata(G2.clade)

mlsa_df <- as_tibble(mlsatree)

`%nin%` <- Negate(`%in%`)



# Subset the nodes of this clades
G2_nodes <- Descendants(mlsatree@phylo, keepnode, type = "all")
G2_nodes <- unlist(G2_nodes)

# Make a tiplist to drop from mlsatree
drop_list <- list()
j = 1

for (i in mlsa_df$node){
  if (i %nin% G2_nodes){
    if (! is.na(mlsa_df[i,'label'])){ # Subset so that the labels are not NA
      if (is.na(mlsa_df[i,'UFboot'])) { # Subset so that the UFboot is NA
        if ( as.character(mlsa_df[i,'label']) != "" ){ # The labels are not empty
          drop_list[[j]] <- as.character(mlsa_df[i,'label'])
          j = j + 1 
        }
      }
    }
  }
}

drop_list <- unlist(drop_list)

G2.clade  <- drop.tip(mlsatree, tip = drop_list)


#write.tree(as.phylo(G2.clade), 'results/combined_mlsa/G2.partition.rxml.treefile')



## Subset the MLSA associated data within the G2_nodes
#G2.clade@data <- mlsatree@data[which(mlsatree@data$node %in% G2_nodes), ]
#
#G2.clade.tree2 <- drop.tip(mlsatree, tip = )
#
#mlsatree@data[mlsatree@phylo$tip.label %in% 

#G2.clade@phylo <- midpoint(G2.clade@phylo)
#G2.clade@data




plot(mlsatree@phylo)


meta <- read.csv('results/combined_mlsa/combined_meta.csv', header=T)

t <- ggtree(G2.clade) %<+% meta + 
  geom_tiplab(aes(label = strain, x = x + 0.002), size=3) +
  geom_tippoint(aes( x = x + 0.001, color=sample_type, fill=sample_type, shape=sample_type),size=2)  +
  scale_fill_manual(values=sampletypecolors, name='Sample Type') +
  scale_color_manual(values=sampletypecolors, name='Sample Type') +
  scale_shape_manual(values=sampletypeshapes, name='Sample Type') +
  labs(color = 'Sample source') +
  theme(legend.position = c(0.2, 0.75)) #+
  #xlim(0, 1)

t
# Update bootstrap values

t$data$bootstrap <- 'Low'
t$data$SH_aLRT[is.na(t$data$SH_aLRT)] <- 100
t$data$UFboot[is.na(t$data$UFboot)] <- 100
t$data[which(t$data$SH_aLRT >= 80 & t$data$UFboot  >= 95),]$bootstrap <- 'High'

t2 <- t + new_scale_color() + 
  geom_tree(aes(color = bootstrap)) +
  scale_color_manual(name = 'Bootstrap support', 
    values = c('High' = 'black', 'Low' = 'grey'),
    labels = c('High' = 'UFBoot > 95% & \nSH-aLRT > 80%', 'Low' = 'other')
  ) +
  geom_treescale(x = 0.005, y = 70) + 
  xlim(0, 0.07)


pdf("Manuscript/Figures/Main/Figure_S1_MLSA_Tree_G2.pdf", width = 16, height = 18)
t2
dev.off()

G2.clade@phylo$Nnode


#######################################
# Cophylogeny of genome and MLSA tree #
#######################################

sampletypecolors <- tableau_color_pal(palette="Tableau 10", type="regular",direction=1)(10)
names(sampletypecolors) <- c("hospital environment","terrestrial","animal","built environment","plant","yellow","clinical","other","soil","unknown")


################
# Plot on tree #
################
snptree <- midpoint(read.tree('results/tree/species1_calls_merged.vcf.nonpassnonhetfilt.vcf.vcffilter.vcf.tab.noref.fasta.treefile'))

plot(snptree)
plot(G2.clade@phylo)

class(G2.clade@phylo)

treelist <- pre.rotate(midpoint.root(G2.clade@phylo), midpoint.root(snptree))

plot(G2.clade@phylo)
plot(treelist[[1]])
plot(snptree)
plot(treelist[[2]])

meta1 <- read.table('results/tree/metadata.table.tsv', header=T, sep='\t')
meta1 <- meta1[, c('strain', 'sample_type')]

meta2 <- read.csv('results/combined_mlsa/combined_meta.csv', header=T)

head(meta2)

treelist

t1 <- ggtree(treelist[[1]], ladderize=F) %<+% meta2[, c('fullstrain', 'strain', 'sample_type')]  + 
  geom_tiplab(aes(x=x+0.0025), size=1.5) +
  geom_tippoint(aes( x = x + 0.001, color=sample_type),size=1.5)  +
  scale_color_manual(values=sampletypecolors) +
  labs(color = 'Sample source') +
  theme(legend.position = c(0.90, 0.15)) 

t1$data$label = t1$data$strain
t1$data = t1$data[,-10]
t1

t2 <- ggtree(treelist[[2]], ladderize=F) %<+% meta2[, c('strain', 'sample_type')] + 
  geom_tiplab(aes( x = x + 0.002), size=1.5) +
  geom_tippoint(aes( x = x + 0.001, color=sample_type),size=1.5)  +
  scale_color_manual(values=sampletypecolors) +
  labs(color = 'Sample source') #+
  #geom_treescale()#+
  #theme(legend.position = c(0.1, 0.75))

t2

# Clinical strains
simple.tanglegram(t1, t2, sample_type, clinical, t2_pad = .5,
                  tiplab = T, lab_pad = 0.015, 
                  t2_y_scale = 1.5, 
                  t2_tiplab_size = 1.5, t2_tiplab_pad = 0.001)


# Custom color
tangler::simple.tanglegram(t1, t2, sample_type, clinical, t2_pad = 0.5,
                           tiplab = T, lab_pad = 0.015, x_hjust = 1, l_color = 'grey',
                           t2_y_pos = 0, t2_y_scale = 1.5, 
                           t2_tiplab_size = 1.5, t2_tiplab_pad = 0.001)

# No custom color, generate random color
tangler::simple.tanglegram(t1, t2, sample_type, clinical, t2_pad = 0.5,
                           tiplab = T, lab_pad = 0.015, x_hjust = 1, 
                           t2_y_pos = 0, t2_y_scale = 1.5, 
                           t2_tiplab_size = 1.5, t2_tiplab_pad = 0.001)

# All strains
tangler::common.tanglegram(t1, t2, column='sample_type', t2_pad = .3,
                  lab_pad = 0.015, tiplab=TRUE,
                  t2_tiplab_size = 1.5, t2_tiplab_pad = 0.001)


 
# Functions
pre.rotate <- function  (tree1, tree2) {
  cophylo <- phytools::cophylo(tree1, tree2)
  rotated_tree1 <- cophylo[[1]][[1]]
  rotated_tree2 <- cophylo[[1]][[2]]
  return(list(rotated_tree1, rotated_tree2))
}




simple.tanglegram <- function (tree1, tree2,  column, value, 
                               t2_pad=0.3, x_hjust=1, lab_pad = 2, 
                               l_color = NA, tiplab=F, t2_y_pos=0,  
                               t2_y_scale=1, t2_branch_scale=1, t2_tiplab_size=1, t2_tiplab_pad = 0) {
  # Update meta column variables for subsetting
  col_name <- deparse(substitute(column))
  parsed_value <- deparse(substitute(value))
  
  
  # Extract tree data
  d1 <- tree1$data
  d2 <- tree2$data
  
  
  # Update the associated variable
  d1$tree <-'t1'
  d2$tree <-'t2'
  
  # Define x coordinate for tree 2
  d2$x <- max(d2$x) - t2_pad*d2$x + max(d1$x)
  d2$y <- d2$y * t2_y_scale
  d2$y <- d2$y + t2_y_pos
  
  tree1$data$x <- tree1$data$x + t2_pad*max(d1$x)
  d1$x <- tree1$data$x
  
  
  # Draw cophylogeny
  pp <- tree1 + geom_tree(data=d2)
  
  # Combine tree associated data.frames
  dd1 <- rbind(d1, d2)
  dd1 <- as.data.frame(dd1[which(dd1$isTip == T),])
  
  # Conditionally join the tips from both tree
  conditional_subset <- dd1[which(dd1[,col_name] == parsed_value), ]
  conditional_subset$lab_x <- conditional_subset$x
  
  # Update label x position
  conditional_subset <- conditional_subset %>%
    dplyr::group_by(label) %>%
    dplyr::mutate(
      lab_x = case_when(
        lab_x == min(lab_x) ~ lab_x + lab_pad,
        lab_x == max(lab_x) ~ lab_x - lab_pad,
        TRUE ~ lab_x
      )
    ) %>%
    dplyr::ungroup()
  
  
  if (is.na(l_color)) {
    pp <- pp + new_scale_color() + ggplot2::geom_line(aes(x = lab_x, y = y, group = label, color = label), data = conditional_subset, show.legend = FALSE) +
    scale_color_viridis_d(option="turbo")  # Use a color scale for discrete colors
  } else {
    pp <- pp + ggplot2::geom_line(aes(lab_x, y, group=label), data=conditional_subset, color=l_color, show.legend = FALSE)
  }
  
  # Show tip-labels
  if (tiplab == T){
    pp + ggtree::geom_tiplab(aes(x = x - t2_tiplab_pad), size=t2_tiplab_size, data=d2, hjust=x_hjust)
  } else {
    pp
  }
  
}


# Clinical strains 
# without defined color
simple.tanglegram(t1, t2, sample_type, clinical, t2_pad = .5,
                  tiplab = T, lab_pad = 0.015, 
                  t2_y_scale = 1.5, 
                  t2_tiplab_size = 1.5, t2_tiplab_pad = 0.001)
# with defined color
simple.tanglegram(t1, t2, sample_type, clinical, t2_pad = .4,
                                     tiplab = T, lab_pad = 0.015, 
                                     t2_y_scale = 1.5, l_color = 'green2',
                                     t2_tiplab_size = 1.5, t2_tiplab_pad = 0.001)





####
common.tanglegram <- function(tree1, tree2, column, sampletypecolors=NA, 
                              t2_pad = 0.3, t2_y_scale = 1, t2_y_pos = 0, 
                              lab_pad = 2, tiplab = FALSE, t2_tiplab_size = 1, 
                              t2_tiplab_pad = 0) {
  
  # Extract tree data
  d1 <- tree1$data
  d2 <- tree2$data
  
  # Update the associated variable
  d1$tree <- 't1'
  d2$tree <- 't2'
  
  # Define x coordinate for tree 2
  d2$x <- max(d2$x) - t2_pad*d2$x + max(d1$x)
  d2$y <- d2$y * t2_y_scale
  d2$y <- d2$y + t2_y_pos
  
  tree1$data$x <- tree1$data$x + t2_pad*max(d2$x)
  d1$x <- tree1$data$x
  
  
  
  
  # Draw cophylogeny
  pp <- tree1 + geom_tree(data=d2, layout = "dendrogram") 
  
  # Merge tree data for tips only
  combined_data <- rbind(d1, d2) %>% filter(isTip == TRUE)
  
  # Create lines connecting the tips and assign colors by the specified column
  combined_data <- combined_data %>%
    group_by(label) %>%
    mutate(
      lab_x = case_when(
        tree == "t1" ~ x + lab_pad,
        tree == "t2" ~ x - lab_pad,
        TRUE ~ x
      )
    ) %>%
    ungroup()
  
  # Add connecting lines colored by the trait category using sampletypecolors
  pp <- pp + 
    geom_line(
      aes(
        x = lab_x, 
        y = y, 
        group = label, 
        color = .data[[column]]
      ),
      data = combined_data
    ) 
  
  if (missing(sampletypecolors) || is.null(sampletypecolors)) {
    pp <- pp + scale_color_viridis_d(option="turbo")   # Use random colors
  } else {
    pp <- pp + scale_color_manual(values = sampletypecolors)  # Use custom colors
  }
    

  
  # Optionally show tip labels
  if (tiplab) {
    pp <- pp + 
      ggtree::geom_tiplab(
        aes(x = x - t2_tiplab_pad),
        size = t2_tiplab_size,
        data = d2,
        hjust = 1
      )
  }
  
  return(pp)
}


common.tanglegram(t1, t2, column='sample_type', sampletypecolors, t2_pad = .2,
                  lab_pad = 0.015, tiplab=TRUE, t2_y_scale = 1.3,
                  t2_tiplab_size = 1.5, t2_tiplab_pad = 0.001)

common.tanglegram(t1, t2, column='sample_type',  t2_pad = .2,
                  lab_pad = 0.015, tiplab=TRUE, t2_y_scale = 1.3,
                  t2_tiplab_size = 1.5, t2_tiplab_pad = 0.001)



