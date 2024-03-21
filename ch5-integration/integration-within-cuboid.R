
## ---- Info ----

## script for
## calculating overall morphological integration
## and
## calculating integration between individual facets
## within the cuboid


## ---- Setup ----

# load packages
library(geomorph)
library(tidyverse)
library(ape)

# load gpa shape data
load("data/import/shape-data/cuboid_species-gpa-coords.R")

# load curve data
load("data/import/curve-and-landmark-data/cuboid_curve-info.R")

# read in curve table
curve_table_cuboid <- 
  read_csv("./data/import/curve-and-landmark-data/cuboid-curve-table.csv")

# read in tree
phylo <- 
  read.tree("./data/import/trees/species-52/consensus-tree-52")


## ---- Define landmarks for articular facets ----

# distal facet
lm_cuboid_distal <- 
  # extract fixed and curve landmarks 
  # which belong to correct facet from curve table 
  # i.e. match correct facet in facet column
  my_curves_cuboid$Curves[which(curve_table_cuboid$Facet%in%"distal")]%>%
  # simplify extracted list of curves 
  # to produce single vector 
  # containing all LMs belonging to dist. medial facet
  unlist(.)%>%
  # remove duplicate LM values
  unique(.)%>%
  # sort list of LMs into ascending order
  sort(.)

# calcaneal facet
lm_cuboid_calc <- 
  my_curves_cuboid$Curves[which(curve_table_cuboid$Facet%in%"calcaneal")]%>%
  unlist(.)%>%
  unique(.)%>%
  sort(.)

# medial facets
# combined medial and distal medial facets
# as distal medial facet is not always present
# and 
lm_cuboid_med <- 
  my_curves_cuboid$Curves[which(curve_table_cuboid$Facet%in%c("medial","dist_medial"))]%>%
  unlist(.)%>%
  unique(.)%>%
  sort(.)


## ---- Integration between distal and calcaneal facets ----

# two block pls
pls_cuboid_df_cf <- 
  integration.test(
    shape_data_means_cuboid[lm_cuboid_distal,,],
    A2 = shape_data_means_cuboid[lm_cuboid_calc,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/cuboid/pls/pls-cuboid-df-cf.txt")
print(summary(pls_cuboid_df_cf))
sink()


# phylogenetic partial least squares
ppls_cuboid_df_cf <- 
  phylo.integration(
    shape_data_means_cuboid[lm_cuboid_distal,,],
    A2 = shape_data_means_cuboid[lm_cuboid_calc,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/cuboid/ppls/ppls-cuboid-df-cf.txt")
print(summary(ppls_cuboid_df_cf))
sink()


## ---- Integration between distal and medial facets ----

# two block pls
pls_cuboid_df_mf <- 
  integration.test(
    shape_data_means_cuboid[lm_cuboid_distal,,],
    A2 = shape_data_means_cuboid[lm_cuboid_med,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/cuboid/pls/pls-cuboid-df-mf.txt")
print(summary(pls_cuboid_df_mf))
sink()


# phylogenetic partial least squares
ppls_cuboid_df_mf <- 
  phylo.integration(
    shape_data_means_cuboid[lm_cuboid_distal,,],
    A2 = shape_data_means_cuboid[lm_cuboid_med,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/cuboid/ppls/ppls-cuboid-df-mf.txt")
print(summary(ppls_cuboid_df_mf))
sink()


## ---- Integration between medial and calcaneal facets ----

# two block pls
pls_cuboid_mf_cf <- 
  integration.test(
    shape_data_means_cuboid[lm_cuboid_med,,],
    A2 = shape_data_means_cuboid[lm_cuboid_calc,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/cuboid/pls/pls-cuboid-mf-cf.txt")
print(summary(pls_cuboid_mf_cf))
sink()


# phylogenetic partial least squares
ppls_cuboid_mf_cf <- 
  phylo.integration(
    shape_data_means_cuboid[lm_cuboid_med,,],
    A2 = shape_data_means_cuboid[lm_cuboid_calc,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/cuboid/ppls/ppls-cuboid-mf-cf.txt")
print(summary(ppls_cuboid_mf_cf))
sink()


## ---- Global integration ----

glob_int_cuboid <- 
  globalIntegration(
    shape_data_means_cuboid,
    ShowPlot = TRUE
  )
