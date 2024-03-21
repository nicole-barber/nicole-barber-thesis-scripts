
## ---- Info ----

## script for
## calculating overall morphological integration
## and
## calculating integration between individual facets
## within the astragalus


## ---- Setup ----

# load packages
library(geomorph)
library(tidyverse)
library(ape)


# load gpa shape data
load("data/import/shape-data/navicular_species-gpa-coords.R")

# load curve data
load("data/import/curve-and-landmark-data/navicular_curve-info.R")

# read in curve table
curve_table_navic <- 
  read_csv("./data/import/curve-and-landmark-data/navicular-curve-table.csv")

# read in tree
phylo <- 
  read.tree("./data/import/trees/species-52/consensus-tree-52")


## ---- Define landmarks for articular facets ----

# proximal facet
lm_navic_prox <- 
  # extract fixed and curve landmarks 
  # which belong to correct facet from curve table 
  # i.e. match correct facet in facet column
  my_curves_navic$Curves[which(curve_table_navic$facet%in%"proximal")]%>%
  # simplify extracted list of curves 
  # to produce single vector 
  # containing all LMs belonging to dist. medial facet
  unlist(.)%>%
  # remove duplicate LM values
  unique(.)%>%
  # sort list of LMs into ascending order
  sort(.)

# lateral cuneiform facet
lm_navic_lat_cun <- 
  my_curves_navic$Curves[which(curve_table_navic$facet%in%"lateral_cuneiform")]%>%
  unlist(.)%>%
  unique(.)%>%
  sort(.)

# medial-intermediate cuneiform facet
lm_navic_med_int <- 
  my_curves_navic$Curves[which(curve_table_navic$facet%in%"medial_intermediate_cuneiform")]%>%
  unlist(.)%>%
  unique(.)%>%
  sort(.)


## ---- Integration between proximal and lateral cuneiform facets ----

# two block pls
pls_navic_pf_lcf <- 
  integration.test(
    shape_data_means_navic[lm_navic_prox,,],
    A2 = shape_data_means_navic[lm_navic_lat_cun,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/navicular/pls/pls-navic-pf-lcf.txt")
print(summary(pls_navic_pf_lcf))
sink()


# phylogenetic partial least squares
ppls_navic_pf_lcf <- 
  phylo.integration(
    shape_data_means_navic[lm_navic_prox,,],
    A2 = shape_data_means_navic[lm_navic_lat_cun,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/navicular/ppls/ppls-navic-pf-lcf.txt")
print(summary(ppls_navic_pf_lcf))
sink()


## ---- Integration between proximal and medial-intermediate cuneiform facet ----

# two block pls
pls_navic_pf_micf <- 
  integration.test(
    shape_data_means_navic[lm_navic_prox,,],
    A2 = shape_data_means_navic[lm_navic_med_int,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/navicular/pls/pls-navic-pf-micf.txt")
print(summary(pls_navic_pf_micf))
sink()


# phylogenetic partial least squares
ppls_navic_pf_micf <- 
  phylo.integration(
    shape_data_means_navic[lm_navic_prox,,],
    A2 = shape_data_means_navic[lm_navic_med_int,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/navicular/ppls/ppls-navic-pf-micf.txt")
print(summary(ppls_navic_pf_micf))
sink()


## ---- Integration between lateral and medial-intermediate cuneiform facets ----

# two block pls
pls_navic_lcf_micf <- 
  integration.test(
    shape_data_means_navic[lm_navic_lat_cun,,],
    A2 = shape_data_means_navic[lm_navic_med_int,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/navicular/pls/pls-navic-lcf-micf.txt")
print(summary(pls_navic_lcf_micf))
sink()


# phylogenetic partial least squares
ppls_navic_lcf_micf <- 
  phylo.integration(
    shape_data_means_navic[lm_navic_lat_cun,,],
    A2 = shape_data_means_navic[lm_navic_med_int,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/navicular/ppls/ppls-navic-lcf-micf.txt")
print(summary(ppls_navic_lcf_micf))
sink()


## ---- Global integration ----

glob_int_navic <- 
  globalIntegration(
    shape_data_means_navic,
    ShowPlot = TRUE
  )
