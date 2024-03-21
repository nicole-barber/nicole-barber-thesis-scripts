
## ---- Info ----

## script for
## calculating morphological integration
## between astragalus, calcaneus, cuboid and navicular


## ---- Setup ----

# load packages
library(geomorph)
library(tidyverse)
library(ape)


# load gpa shape data

# astragalus
load("data/import/shape-data/astragalus_species-gpa-coords.R")

# calcaneus
load("data/import/shape-data/calcaneus_species-gpa-coords.R")

# cuboid
load("data/import/shape-data/cuboid_species-gpa-coords.R")

# navicular
load("data/import/shape-data/navicular_species-gpa-coords.R")


# read in tree
phylo <- 
  read.tree("./data/import/trees/species-52/consensus-tree-52")


## ---- Astragalus and calcaneus integration ----

# two block pls
pls_astrag_calc <-
  integration.test(
    shape_data_means_astrag,
    A2 = shape_data_means_calc,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
# open sink
sink("./analysis/integration-between-elements/pls/pls-astragalus-calcaneus.txt")
print(summary(pls_astrag_calc))
# close sink
sink()


# phylogenetic partial least squares
ppls_astrag_calc <-
  phylo.integration(
    shape_data_means_astrag,
    A2 = shape_data_means_calc,
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
# open sink
sink("./analysis/integration-between-elements/ppls/ppls-astragalus-calcaneus.txt")
print(summary(ppls_astrag_calc))
# close sink
sink()


## ---- Astraagalus and cuboid integration ----

# two block pls
pls_astrag_cuboid <-
  integration.test(
    shape_data_means_astrag,
    A2 = shape_data_means_cuboid,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
# open sink
sink("./analysis/integration-between-elements/pls/pls-astragalus-cuboid.txt")
print(summary(pls_astrag_cuboid))
# close sink
sink()


# phylogenetic partial least squares
ppls_astrag_cuboid <- 
  phylo.integration(
    shape_data_means_astrag,
    A2 = shape_data_means_cuboid,
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
# open sink
sink("./analysis/integration-between-elements/ppls/ppls-astragalus-cuboid")
print(summary(ppls_astrag_cuboid))
sink()


## ---- Astragalus and navicular integration ----

# two block pls
pls_astrag_navic <-
  integration.test(
    shape_data_means_astrag,
    A2 = shape_data_means_navic,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
# open sink
sink("./analysis/integration-between-elements/pls/pls-astragalus-navicular.txt")
print(summary(pls_astrag_navic))
# close sink
sink()


# phylogenetic partial least squares
ppls_astrag_navic <-
  phylo.integration(
    shape_data_means_astrag,
    A2 = shape_data_means_navic,
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
# open sink
sink("./analysis/integration-between-elements/ppls/ppls-astragalus-navicular.txt")
print(summary(ppls_astrag_navic))
# close sink
sink()


## ---- Calcaneus and cuboid integration ----

# two block pls
pls_calc_cuboid <-
  integration.test(
    shape_data_means_calc,
    A2 = shape_data_means_cuboid,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
# open sink
sink("./analysis/integration-between-elements/pls/pls-calcaneus-cuboid.txt")
print(summary(pls_calc_cuboid))
# close sink
sink()


# phylogenetic partial least squares
ppls_calc_cuboid <-
  phylo.integration(
    shape_data_means_calc,
    A2 = shape_data_means_cuboid,
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
# open sink
sink("./analysis/integration-between-elements/ppls/ppls-calcaneus-cuboid.txt")
print(summary(ppls_calc_cuboid))
# close sink
sink()


## ---- Calcaneus and navicular integration ----

# two block pls
pls_calc_navic <- 
  integration.test(
    shape_data_means_calc,
    A2 = shape_data_means_navic,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-between-elements/pls/pls-calcaneus-navicular.txt")
print(summary(pls_calc_navic))
sink()


# phylogenetic partial least squares
ppls_calc_navic <- 
  phylo.integration(
    shape_data_means_calc,
    A2 = shape_data_means_navic,
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-between-elements/ppls/ppls-calcaneus-navicular.txt")
print(summary(ppls_calc_navic))
sink()


## ---- Cuboid and navicular integration ----

# two block pls
pls_cuboid_navic <- 
  integration.test(
    shape_data_means_cuboid,
    A2 = shape_data_means_navic,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write result to .txt file
sink("./analysis/integration-between-elements/pls/pls-cuboid-navicular.txt")
print(summary(pls_cuboid_navic))
sink()


# phylogenetic partial least squares
ppls_cuboid_navic <- 
  phylo.integration(
    shape_data_means_cuboid,
    A2 = shape_data_means_navic,
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-between-elements/ppls/ppls-cuboid-navicular.txt")
print(summary(ppls_cuboid_navic))
sink()


## ---- Compare integration outputs ----

pls_comparison <-
  compare.pls(
    pls_astrag_calc,
    pls_astrag_cuboid,
    pls_astrag_navic,
    pls_calc_cuboid,
    pls_calc_navic,
    pls_cuboid_navic
  )
