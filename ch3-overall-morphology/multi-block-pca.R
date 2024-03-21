
## ---- Info ----

## script for
## performing regularised consensus PCA (RCPCA)
## using morphoBlocks

## requires
## reampled and slid landmarks
## OR here gpa adjusted landmarks as used here
## vector of centroid sizes from gpa 
## if using already gpa transformed coords


## ---- Setup ----

# if morphoBlocks not already installed
# devtools::install_github("aharmer/morphoBlocks")

# load packages
library(morphoBlocks)
library(tidyverse)


# source custom scripts
source("scripts/plotting-aesthetics.R")


# import shape data

# astragalus
load("data/import/shape-data/astragalus_species-gpa-coords.R")

# calcaneus
load("data/import/shape-data/calcaneus_species-gpa-coords.R")

# cuboid
load("data/import/shape-data/cuboid_species-gpa-coords.R")

# navicular
load("data/import/shape-data/navicular_species-gpa-coords.R")


# import size data
# astragalus
load("data/import/size-data/astragalus_species-csize.R")

# calcaneus
load("data/import/size-data/calcaneus_species-csize.R")

# cuboid
load("data/import/size-data/cuboid_species-csize.R")

# navicular
load("data/import/size-data/navicular_species-csize.R")


# import species data
species_data <- 
  read_csv("data/import/species-data/species-data.csv")


## ---- Generate data blocks ----

# astragalus
block_astragalus <-
  formatBlock(
    # 3D array containing astragalus landmarks
    shape_data_means_astrag,
    # include vector of centroid size as GPA already performed
    cs = csize_means_astrag,
    # number of landmark dimensions
    k = 3,
    # do not need to perform GPA as already done to exctract species means
    gpa = FALSE
  )

# calcaneus
block_calcaneus <-
  formatBlock(
    # 3D array containing astragalus landmarks
    shape_data_means_calc,
    # include vector of centroid size as GPA already performed
    cs = csize_means_calc,
    # number of landmark dimensions
    k = 3,
    # do not need to perform GPA as already done to exctract species means
    gpa = FALSE
  )

# cuboid
block_cuboid <-
  formatBlock(
    # 3D array containing astragalus landmarks
    shape_data_means_cuboid,
    # include vector of centroid size as GPA already performed
    cs = csize_means_cuboid,
    # number of landmark dimensions
    k = 3,
    # do not need to perform GPA as already done to exctract species means
    gpa = FALSE
  )

# navicular
block_navicular <-
  formatBlock(
    # 3D array containing astragalus landmarks
    shape_data_means_navic,
    # include vector of centroid size as GPA already performed
    cs = csize_means_navic,
    # number of landmark dimensions
    k = 3,
    # do not need to perform GPA as already done to exctract species means
    gpa = FALSE
  )


## ---- Combine data blocks ----

# combine individual element blocks into a list
blocklist_tarsals <-
  combineBlocks(
    blocks = c(
      block_astragalus,
      block_calcaneus,
      block_cuboid,
      block_navicular
    )
  )


## ---- Analyse data blocks ----

tarsals_result <-
  analyseBlocks(
    blocklist_tarsals,
    # specify number of principal components - same as no. of species
    ncomp = 52
  )

## ******************************************************************* ##
## You will get the following warning message:                         ##
## In rgcca(blockList, C = C, tau = tau, ncomp = rep(ncomp, J + 1),  : ##
##        Argument C is deprecated, use connection instead.            ##
## ******************************************************************* ##

# plot consensus space with GC1 and GC2
scoresPlot(
  tarsals_result,
  pcol = palette_morphoblocks
)


# plot loadings
loadingsPlot(
  tarsals_result,
  comp = 1,
  cex.3d = 10
)



