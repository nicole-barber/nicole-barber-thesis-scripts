
## ---- Info ----

## script for
## pca
## phylogenetic signal (kmult)
## allometry

## requires
## resampled, slid and gpa adjusted landmarks
## phylogenetic tree in Newick format

## this script
## is for the astragalus


## ---- Setup ----

# load packages
library(tidyverse)
library(geomorph)
library(phytools)
library(rgl)
library(ape)
library(Rvcg)


# read in shape data
load("./data/import/shape-data/astragalus_species-gpa-coords.R")

# read in size data
load("./data/import/size-data/astragalus_species-csize.R")

# read in tree
phylo <- 
  read.tree(
    "./data/import/trees/species-52/consensus-tree-52"
  )

# read in species data
species_data <- 
  read_csv("./data/import/species-data/species-data.csv")


## ---- PCA ----

# perform PCA on GPA adjusted landmarks
pca_astrag <- 
  gm.prcomp(
    shape_data_means_astrag
  )

# write results of pca to txt file
# open sink
sink("./analysis/pca/astragalus-pca.txt")
# print summary of pca - it will write to that file
print(summary(pca_astrag))
# close sink
sink()


## plot

# extract PC scores for 95% of variance
pcs95_astrag <- 
  pca_astrag$x[,1:23]

# format PC scores for relevant variance as a tibble
pcscores_astrag <- 
  as_tibble(pcs95_astrag) %>%
  # add in an extra column with relevant label for each coordinate set
  # to facilitate merging pc scores with other data
  mutate(Species_ID = rownames(pcs95_astrag))

# pull out PCs 1 to 4 for plotting                   
pcscores_plot_astrag <- 
  pcscores_astrag %>% 
  select(Comp1, Comp2, Comp3, Comp4, Species_ID)

# join pca results with other species data
pcscores_plot_astrag <- 
  left_join(pcscores_plot_astrag, species_data)


# **************************************************************** #
# below is ggplot code for a basic morphospace                     #
# with species labels so you can see what you're looking at        #
# use script "plot-phylomorphospace.R                              #
# for plot including phylogeny                                     #
# **************************************************************** #


# source aesthetics for plotting
source("./scripts/plotting-aesthetics.R")


# plot PC1 and PC2
ggplot(
  pcscores_plot_astrag, 
  # specify PC1 on x axis, PC2 on y axis
  # point shape according to family
  # point colour according to genus
  aes(
    # PC1 on x axis
    x = Comp1, 
    # PC2 on y axis
    y = Comp2,
    # point shape according to family
    shape = Family,
    # point stroke weight according to family
    stroke = Family,
    # point colour according to genus
    colour = Genus,
    # fill colour according to genus 
    fill = Genus,
    label = Species_ID
  )
) +
  
  # add points and specify point size
  geom_point(
    size = 2.5
  ) +
  
  # manually select shapes for each family
  scale_shape_manual(
    values = point_shape_family13
  ) +
  
  # manually set point stroke weight for each family
  scale_discrete_manual(
    # call stroke aes specified in geom_points()
    aesthetics = "stroke",
    # set appropriate weight for each point shape
    values = point_stroke_family13
  ) +
  
  # manually select colour for each genus
  scale_color_manual(
    values = palette43
  ) +
  
  scale_fill_manual(
    values = palette43
  ) +
  
  # specify x and y labels and plot title
  labs(
    x = "Principal Component 1 (30.2%)",
    y = "Principal Component 2 (10.3%)",
    title = "Primate Astragalus Morphospace"
  ) +
  
  geom_text(
    size = 3,
    check_overlap = TRUE
  ) +
  
  # set theme for plot (white background with no gridlines)
  theme_bw()


# export  plot to file
ggsave(
  filename = 
    "./analysis/pca/astragalus-pc1-pc2.png",
  device = "png",
  width =  12,
  height = 8.5,
  units = "in",
  dpi = 300
)


# plot PC3 and PC4
ggplot(
  pcscores_plot_astrag, 
  # specify PC1 on x axis, PC2 on y axis
  # point shape according to family
  # point colour according to genus
  aes(
    # PC1 on x axis
    x = Comp3, 
    # PC2 on y axis
    y = Comp4,
    # point shape according to family
    shape = Family,
    # point stroke weight according to family
    stroke = Family,
    # point colour according to genus
    colour = Genus,
    # fill colour according to genus 
    fill = Genus,
    label = Species_ID
  )
) +
  
  # add points and specify point size
  geom_point(
    size = 2.5
  ) +
  
  # manually select shapes for each family
  scale_shape_manual(
    values = point_shape_family13
  ) +
  
  # manually set point stroke weight for each family
  scale_discrete_manual(
    # call stroke aes specified in geom_points()
    aesthetics = "stroke",
    # set appropriate weight for each point shape
    values = point_stroke_family13
  ) +
  
  # manually select colour for each genus
  scale_color_manual(
    values = palette43
  ) +
  
  scale_fill_manual(
    values = palette43
  ) +
  
  # specify x and y labels and plot title
  labs(
    x = "Principal Component 3 (10.2%)",
    y = "Principal Component 4 (9.2%)",
    title = "Primate Astragalus Morphospace"
  ) +
  
  geom_text(
    size = 2,
    check_overlap = TRUE
  ) +
  
  # set theme for plot (white background with no gridlines)
  theme_bw()


# export  plot to file
ggsave(
  filename = 
    "./analysis/pca/astragalus-pc3-pc4.png",
  device = "png",
  width =  12,
  height = 8.5,
  units = "in",
  dpi = 300
)


## ---- Phylogenetic signal ----

## calculate degree of phylogenetic signal in the data 
## using the Kmult statistic

physignal_astrag <- 
  physignal(
    # use Procrustes coords in 3D array
    shape_data_means_astrag, 
    # specify tree
    phylo, 
    print.progress = TRUE
  )

# write results to txt file
# open sink
sink("./analysis/kmult/phylogenetic-signal-astragalus.txt")
print(summary(physignal_astrag))
# close sink
sink()


## ---- Allometry ----

## test relationship between size (CSize) and procrustes shape data (coords)
## in a phylogenetic framework

# put relevant data into geomorph data frame
# for for easier function use
gdf_astrag <-
  geomorph.data.frame(
    Csize = csize_means_astrag, 
    shape = shape_data_means_astrag, 
    phy = phylo
  )

# perform regression of shape against size
# in a phylogenetic framework
pgls_astrag <-
  procD.pgls(
    # formula for linear model
    shape ~ log(Csize),
    # specify tree
    phy = phy,
    # specify geomoprh dataframe to use
    data = gdf_astrag,
    print.progress = TRUE
  )

# Evaluate the assumptions of the model
# par(mfrow = c(2,2))
# plot(pgls_astrag)
# par(mfrow = c(1,1))

# check residuals vs fitted
plot(
  pgls_astrag, 
  method="PredLine"
)

# check outliers
# plot(
#  pgls_calc,
#  outliers = TRUE
# )

# plot fit
plotAllometry(
  pgls_astrag, 
  size = csize_means_astrag, 
  # plot common allometric component (basically like PC1 for size)
  method = "CAC"
)

# write results to txt file
# open sink
sink("./analysis/pgls-allometry/pgls-astragalus.txt")
print(summary(pgls_astrag))
# close sink
sink()


