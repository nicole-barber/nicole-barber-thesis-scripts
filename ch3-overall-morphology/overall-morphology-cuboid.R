
## ---- Info ----

## script for
## pca
## phylogenetic signal (kmult)
## allometry

## requires
## resampled, slid and gpa adjusted landmarks
## phylogenetic tree in Newick format

## this script
## is for the cuboid


## ---- Setup ----

# load packages
library(tidyverse)
library(geomorph)
library(phytools)
library(rgl)
library(ape)
library(Rvcg)


# read in shape data
load("./data/import/shape-data/cuboid_species-gpa-coords.R")

# read in size data
load("./data/import/size-data/cuboid_species-csize.R")

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
pca_cuboid <- 
  gm.prcomp(
    shape_data_means_cuboid
  )

# write results of pca to txt file
# open sink
sink("./analysis/pca/cuboid-pca.txt")
# print summary of pca - it will write to that file
print(summary(pca_cuboid))
# close sink
sink()


## plot

# extract PC scores for 95% of variance
pcs95_cuboid <- 
  pca_cuboid$x[,1:11]

# format PC scores for relevant variance as a tibble
pcscores_cuboid <- 
  as_tibble(pcs95_cuboid) %>%
  # add in an extra column with relevant label for each coordinate set
  # to facilitate merging pc scores with other data
  mutate(Species_ID = rownames(pcs95_cuboid))

# pull out PCs 1 to 4 for plotting                   
pcscores_plot_cuboid <- 
  pcscores_cuboid %>% 
  select(Comp1, Comp2, Comp3, Comp4, Species_ID)

# join pca results with other species data
pcscores_plot_cuboid <- 
  left_join(pcscores_plot_cuboid, species_data)


# **************************************************************** #
# below is ggplot code for a basic morphospace                     #
# use script "plot-phylomorphospace.R                              #
# for plot including phylogeny                                     #
# **************************************************************** #


# source aesthetics for plotting
source("./scripts/plotting-aesthetics.R")


# plot PC1 and PC2
ggplot(
  pcscores_plot_cuboid, 
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
    fill = Genus
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
    x = "Principal Component 1 (37.7%)",
    y = "Principal Component 2 (17.8%)",
    title = "Primate Cuboid Morphospace"
  ) +
  
  # set theme for plot (white background with no gridlines)
  theme_bw()


# export  plot to file
ggsave(
  filename = 
    "./analysis/pca/cuboid-pc1-pc2.png",
  device = "png",
  width =  12,
  height = 8.5,
  units = "in",
  dpi = 300
)


# plot PC3 and PC4
ggplot(
  pcscores_plot_cuboid, 
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
    fill = Genus
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
    x = "Principal Component 3 (11.9%)",
    y = "Principal Component 4 (6.7%)",
    title = "Primate Cuboid Morphospace"
  ) +
  
  # set theme for plot (white background with no gridlines)
  theme_bw()


# export  plot to file
ggsave(
  filename = 
    "./analysis/pca/cuboid-pc3-pc4.png",
  device = "png",
  width =  12,
  height = 8.5,
  units = "in",
  dpi = 300
)


## ---- Phylogenetic signal ----

## calculate degree of phylogenetic signal in the data 
## using the Kmult statistic

physignal_cuboid <- 
  physignal(
    # use Procrustes coords in 3D array
    shape_data_means_cuboid, 
    # specify tree
    phylo, 
    print.progress = TRUE
  )

# write results to txt file
# open sink
sink("./analysis/kmult/phylogenetic-signal-cuboid.txt")
print(summary(physignal_cuboid))
# close sink
sink()


## ---- Allometry ----

## test relationship between size (CSize) and procrustes shape data (coords)
## in a phylogenetic framework

# put relevant data into geomorph data frame
# for for easier function use
gdf_cuboid <-
  geomorph.data.frame(
    Csize = csize_means_cuboid, 
    shape = shape_data_means_cuboid, 
    phy = phylo
  )

# perform regression of shape against size
# in a phylogenetic framework
pgls_cuboid <-
  procD.pgls(
    # formula for linear model
    shape ~ log(Csize),
    # specify tree
    phy = phy,
    # specify geomoprh dataframe to use
    data = gdf_cuboid,
    print.progress = TRUE
  )

# Evaluate the assumptions of the model
# par(mfrow = c(2,2))
# plot(pgls_cuboid)
# par(mfrow = c(1,1))

# check residuals vs fitted
plot(
  pgls_cuboid, 
  method="PredLine"
)

# check outliers
# plot(
#  pgls_calc,
#  outliers = TRUE
# )

# plot fit
plotAllometry(
  pgls_cuboid, 
  size = csize_means_cuboid, 
  # plot common allometric component (basically like PC1 for size)
  method = "CAC"
)

# write results to txt file
# open sink
sink("./analysis/pgls-allometry/pgls-cuboid.txt")
print(summary(pgls_cuboid))
# close sink
sink()

