

## ---- Info ----

## script for
## estimating disparity in landmark data (morpho.disparity)
## estimating evolutionary rates from landmark data (compare.evol.rates)

## requires
## GPA adjusted landmarks
## phylogeny in Newick format


## ---- Setup ----

# load packages
library(ape)
library(geomorph)

# load gpa shape data
load("./data/import/shape-data/calcaneus_species-gpa-coords.R")

# read in tree
phylo <- 
  read.tree("./data/import/trees/species-52/consensus-tree-52")

# read in species data
species_data <- 
  read_csv("./data/import/species-data/species-data.csv")


## ---- Evolutionary rates ----

## rates per family

# assign family as factor
fam_factor <- 
  species_data$Family

# add names
names(fam_factor) <- 
  species_data$Species_ID

# calculate rates per family
rates_fam_calc <- 
  compare.evol.rates(
    # soecify shape data (Procrustes aligned coords)
    shape_data_means_calc,
    # specify tree
    phylo,
    # specify groups you want to compare between
    fam_factor,
    print.progress = TRUE
  )

# write results summary to txt file
sink("./analysis/geomorph-rates/rates-family-calcaneus.txt")
print(summary(rates_fam_calc))
sink()


## rates per infraorder
## (platyrrhine, catarrhine, lemurs, lorises, tarsiers)

# assign infraorder as factor
infra_factor <- 
  species_data$Infraorder

# add names
names(infra_factor) <- 
  species_data$Species_ID 

# calculate rates per infraorder
rates_infra_calc <- 
  compare.evol.rates(
    # specify shape data (Procrustes aligned coords)
    shape_data_means_calc,
    # specify tree
    phylo,
    # specify groups to compare between
    infra_factor,
    print.progress = TRUE
  )

# write results summary to txt file
sink("./analysis/geomorph-rates/rates-infraorder-calcaneus.txt")
print(summary(rates_infra_calc))
sink()


## rates per parvorder
## (platyrrhine, catarrhine, strepsirrhine, tarsiers)

# assign parvorder as factor
parv_factor <- 
  species_data$Parvorder

# add names
names(parv_factor) <- 
  species_data$Species_ID 

## calculate rates per parvorder
rates_parv_calc <- 
  compare.evol.rates(
    # specify shape data (Procrustes aligned coords)
    shape_data_means_calc,
    # specify tree
    phylo,
    # specify groups to compare between
    parv_factor,
    print.progress = TRUE
  )

# write results summary to txt file
sink("./analysis/geomorph-rates/rates-parvorder-calcaneus.txt")
print(summary(rates_parv_calc))
sink()


## ---- Morphological disparity ----

## disparity per family

disp_fam_calc <- 
  morphol.disparity(
    # specify data
    # ~ 1 tells it to use overall mean
    shape_data_means_calc ~ 1,
    # specify groups you want to compare
    groups = species_data$Family
  )

# write results summary to txt file
sink("./analysis/disparity/disparity-family-calcaneus.txt")
print(summary(disp_fam_calc))
sink()


## disparity per infraorder (ish)
## platyrrhine, catarrhine, Malagasy lemurs, lorises

disp_infra_calc <- 
  morphol.disparity(
    # specify data
    # ~ 1 tells it to use overall mean
    shape_data_means_calc ~ 1,
    # specify groups you want to compare
    groups = species_data$Infraorder
  )

# write results summary to txt file
sink("./analysis/disparity/disparity-infraorder-calcaneus.txt")
print(summary(disp_infra_calc))
sink()


## disparity per parvorder (ish)
## platyrrhine, catarrhine, strepsirrhine

disp_parv_calc <- 
  morphol.disparity(
    # specify data
    # ~ 1 tells it to use overall mean
    shape_data_means_calc ~ 1,
    # specify groups you want to compare                                           
    groups = species_data$Parvorder
  )

# write results summary to txt file
sink("./analysis/disparity/disparity-parvorder-calcaneus.txt")
print(summary(disp_parv_calc))
sink()

