
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
load("data/import/shape-data/calcaneus_species-gpa-coords.R")

# load curve data
load("data/import/curve-and-landmark-data/calcaneus_curve-info.R")
 
# read in curve table
curve_table_calc <- 
  read_csv("./data/import/curve-and-landmark-data/calcaneus-curve-table.csv")

# read in tree
phylo <- 
  read.tree("./data/import/trees/species-52/consensus-tree-52")


## ---- Define landmarks for articular facets ----

# ectal facet
lm_calc_ectal <- 
  # extract fixed and curve landmarks 
  # which belong to correct facet from curve table 
  # i.e. match correct facet column
  my_curves_calc$Curves[which(curve_table_calc$facet%in%"ectal")]%>%
  # simplify extracted list of curves 
  # to produce single vector 
  # containing all LMs belonging to dist. medial facet
  unlist(.)%>%
  # remove duplicate LM values
  unique(.)%>%
  # sort list of LMs into ascending order
  sort(.)

# sustentaculum facet
lm_calc_sust <- 
  my_curves_calc$Curves[which(curve_table_calc$facet%in%"sustentaculum")]%>%
  unlist(.)%>%
  unique(.)%>%
  sort(.)

# cuboid facet
lm_calc_cuboid <- 
  my_curves_calc$Curves[which(curve_table_calc$facet%in%"cuboid")]%>%
  unlist(.)%>%
  unique(.)%>%
  sort(.)

# posterior surface
lm_calc_post <- 
  my_curves_calc$Curves[which(curve_table_calc$facet%in%"posterior_surface")]%>%
  unlist(.)%>%
  unique(.)%>%
  sort(.)


## ---- Integration between ectal and sustentacular facets ----

# two block pls
pls_calc_ef_sf <- 
  integration.test(
    shape_data_means_calc[lm_calc_ectal,,],
    A2 = shape_data_means_calc[lm_calc_sust,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/calcaneus/pls/pls-calcaneus-ef-sf.txt")
print(summary(pls_calc_ef_sf))
sink()


# phylogenetic partial least squares
ppls_calc_ef_sf <- 
  phylo.integration(
    shape_data_means_calc[lm_calc_ectal,,],
    A2 = shape_data_means_calc[lm_calc_sust,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/calcaneus/ppls/ppls-calcaneus-ef-sf.txt")
print(summary(ppls_calc_ef_sf))
sink()


## ---- Integration between ectal and cuboid facets ----

# two block pls
pls_calc_ef_cf <- 
  integration.test(
    shape_data_means_calc[lm_calc_ectal,,],
    A2 = shape_data_means_calc[lm_calc_cuboid,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/calcaneus/pls/pls-calcaneus-ef-cf.txt")
print(summary(pls_calc_ef_cf))
sink()


# phylogenetic partial least squares
ppls_calc_ef_cf <- 
  phylo.integration(
    shape_data_means_calc[lm_calc_ectal,,],
    A2 = shape_data_means_calc[lm_calc_cuboid,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/calcaneus/ppls/ppls-calcaneus-ef-cf.txt")
print(summary(ppls_calc_ef_cf))
sink()


## ---- Integration between ectal facet and posterior surface ----

# two block pls
pls_calc_ef_ps <- 
  integration.test(
    shape_data_means_calc[lm_calc_ectal,,],
    A2 = shape_data_means_calc[lm_calc_post,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/calcaneus/pls/pls-calcaneus-ef-ps.txt")
print(summary(pls_calc_ef_ps))
sink()


# phylogenetic partial least squares
ppls_calc_ef_ps <- 
  phylo.integration(
    shape_data_means_calc[lm_calc_ectal,,],
    A2 = shape_data_means_calc[lm_calc_post,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/calcaneus/ppls/ppls-calcaneus-ef-ps.txt")
print(summary(ppls_calc_ef_ps))
sink()


## ---- Integration between sustentacular facet and cuboid facet ----

# two block pls
pls_calc_sf_cf <- 
  integration.test(
    shape_data_means_calc[lm_calc_sust,,],
    A2 = shape_data_means_calc[lm_calc_cuboid,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/calcaneus/pls/pls-calcaneus-sf-cf.txt")
print(summary(pls_calc_sf_cf))
sink()


# phylogenetic partial least squares
ppls_calc_sf_cf <- 
  phylo.integration(
    shape_data_means_calc[lm_calc_sust,,],
    A2 = shape_data_means_calc[lm_calc_cuboid,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/calcaneus/ppls/ppls-calcaneus-sf-cf.txt")
print(summary(ppls_calc_sf_cf))
sink()


## ---- Integration between sustenticular facet and posterior surface ----

# two block pls
pls_calc_sf_ps <- 
  integration.test(
    shape_data_means_calc[lm_calc_sust,,],
    A2 = shape_data_means_calc[lm_calc_post,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/calcaneus/pls/pls-calcaneus-sf-ps.txt")
print(summary(pls_calc_sf_ps))
sink()


# phylogenetic partial least squares
ppls_calc_sf_ps <- 
  phylo.integration(
    shape_data_means_calc[lm_calc_sust,,],
    A2 = shape_data_means_calc[lm_calc_post,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/calcaneus/ppls/ppls-calcaneus-sf-ps.txt")
print(summary(ppls_calc_sf_ps))
sink()


## ---- Integration between cuboid facet and posterior surface ----

# two block pls
pls_calc_cf_ps <- 
  integration.test(
    shape_data_means_calc[lm_calc_cuboid,,],
    A2 = shape_data_means_calc[lm_calc_post,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/calcaneus/pls/pls-calcaneus-cf-ps.txt")
print(summary(pls_calc_cf_ps))
sink()


# phylogenetic partial least squares
ppls_calc_cf_ps <- 
  phylo.integration(
    shape_data_means_calc[lm_calc_cuboid,,],
    A2 = shape_data_means_calc[lm_calc_post,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/calcaneus/ppls/ppls-calcaneus-cf-ps.txt")
print(summary(ppls_calc_cf_ps))
sink()


## ---- Global integration ----

glob_int_calc <- 
  globalIntegration(
    shape_data_means_calc,
    ShowPlot = TRUE
  )
