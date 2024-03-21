
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
load("data/import/shape-data/astragalus_species-gpa-coords.R")

# load curve data
load("data/import/curve-and-landmark-data/astragalus_curve-info.R")

# read in curve table
curve_table_astrag <- 
  read_csv("./data/import/curve-and-landmark-data/astragalus-curve-table.csv")

# read in tree
phylo <- 
  read.tree("./data/import/trees/species-52/consensus-tree-52")


## ---- Define landmarks for articular facets ----

# medial malleolar facet
lm_astrag_med_mal <- 
  # extract fixed and curve landmarks 
  # which belong to correct facet from curve table 
  # i.e. match correct facet column
  my_curves_astrag$Curves[which(curve_table_astrag$facet%in%"medial_malleolar")]%>%
  # simplify extracted list of curves 
  # to produce single vector 
  # containing all LMs belonging to dist. medial facet
  unlist(.)%>%
  # remove duplicate LM values
  unique(.)%>%
  # sort list of LMs into ascending order
  sort(.)

# lateral malleolar facet
lm_astrag_lat_mal <- 
  my_curves_astrag$Curves[which(curve_table_astrag$facet%in%"lateral_malleolar")]%>%
  unlist(.)%>%
  unique(.)%>%
  sort(.)

# trochlea
lm_astrag_troc <- 
  my_curves_astrag$Curves[which(curve_table_astrag$facet%in%"trochlea")]%>%
  unlist(.)%>%
  unique(.)%>%
  sort(.)

# navicular facet
lm_astrag_navic <- 
  my_curves_astrag$Curves[which(curve_table_astrag$facet%in%"navicular")]%>%
  unlist(.)%>%
  unique(.)%>%
  sort(.)

# sustentacular facet
lm_astrag_sust <- 
  my_curves_astrag$Curves[which(curve_table_astrag$facet%in%"sustentaculum")]%>%
  unlist(.)%>%
  unique(.)%>%
  sort(.)

# ectal facet
lm_astrag_ectal <- 
  my_curves_astrag$Curves[which(curve_table_astrag$facet%in%"ectal")]%>%
  unlist(.)%>%
  unique(.)%>%
  sort(.)

# flexor hallucis longus groove
lm_astrag_fhlg <- 
  my_curves_astrag$Curves[which(curve_table_astrag$facet%in%"fhlg")]%>%
  unlist(.)%>%
  unique(.)%>%
  sort(.)


## ---- Integration between medial malleolar and lateral malleolar facets ----

# two block pls
pls_astrag_mmf_lmf <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_med_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_lat_mal,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-mmf-lmf.txt")
print(summary(pls_astrag_mmf_lmf))
sink()


# phylogenetic partial least squares
ppls_astrag_mmf_lmf <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_med_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_lat_mal,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-mmf-lmf.txt")
print(summary(ppls_astrag_mmf_lmf))
sink()


## ---- Integration between medial malleolar and ectal facets ----

# two block pls
pls_astrag_mmf_ef <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_med_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_ectal,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-mmf-ef.txt")
print(summary(pls_astrag_mmf_ef))
sink()


# phylogenetic partial least squares
ppls_astrag_mmf_ef <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_med_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_ectal,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-mmf-ef.txt")
print(summary(ppls_astrag_mmf_ef))
sink()


## ---- Integration between medial malleolar and navicular facets ----

# two block pls
pls_astrag_mmf_nf <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_med_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_navic,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-mmf-nf.txt")
print(summary(pls_astrag_mmf_nf))
sink()


# phylogenetic partial least squares
ppls_astrag_mmf_nf <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_med_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_navic,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-mmf-nf.txt")
print(summary(ppls_astrag_mmf_nf))
sink()


## ---- Integration between medial malleolar and sustentacular facets ----

# two block pls
pls_astrag_mmf_sf <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_med_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_sust,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-mmf-sf.txt")
print(summary(pls_astrag_mmf_sf))
sink()


# phylogenetic partial least squares
ppls_astrag_mmf_sf <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_med_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_sust,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-mmf-sf.txt")
print(summary(ppls_astrag_mmf_sf))
sink()


## ---- Integration between medial malleolar and trochlear facets ----

# two block pls
pls_astrag_mmf_tf <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_med_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_troc,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-mmf-tf.txt")
print(summary(pls_astrag_mmf_tf))
sink()


# phylogenetic partial least squares
ppls_astrag_mmf_tf <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_med_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_troc,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-mmf-tf.txt")
print(summary(ppls_astrag_mmf_tf))
sink()


## ---- Integration between medial malleolar facet and fhlg ----

# two block pls
pls_astrag_mmf_fhlg <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_med_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_fhlg,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-mmf-fhlg.txt")
print(summary(pls_astrag_mmf_fhlg))
sink()


# phylogenetic partial least squares
ppls_astrag_mmf_fhlg <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_med_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_fhlg,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-mmf-fhlg.txt")
print(summary(ppls_astrag_mmf_fhlg))
sink()


## ---- Integration between lateral malleolar and ectal facets ----

# two block pls
pls_astrag_lmf_ef <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_lat_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_ectal,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-lmf-ef.txt")
print(summary(pls_astrag_lmf_ef))
sink()


# phylogenetic partial least squares
ppls_astrag_lmf_ef <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_lat_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_ectal,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-lmf-ef.txt")
print(summary(ppls_astrag_lmf_ef))
sink()


## ---- Integration between lateral malleolar and navicular facets ----

# two block pls
pls_astrag_lmf_nf <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_lat_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_navic,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-lmf-nf.txt")
print(summary(pls_astrag_lmf_nf))
sink()


# phylogenetic partial least squares
ppls_astrag_lmf_nf <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_lat_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_navic,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-lmf-nf.txt")
print(summary(ppls_astrag_lmf_nf))
sink()


## ---- Integration between lateral malleolar and sustentacular facets ----

# two block pls
pls_astrag_lmf_sf <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_lat_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_sust,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-lmf-sf.txt")
print(summary(pls_astrag_lmf_sf))
sink()


# phylogenetic partial least squares
ppls_astrag_lmf_sf <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_lat_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_sust,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-lmf-sf.txt")
print(summary(ppls_astrag_lmf_sf))
sink()


## ---- Integration between lateral malleolar and trochlear facets ----

# two block pls
pls_astrag_lmf_tf <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_lat_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_troc,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-lmf-tf.txt")
print(summary(pls_astrag_lmf_tf))
sink()


# phylogenetic partial least squares
ppls_astrag_lmf_tf <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_lat_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_troc,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-lmf-tf.txt")
print(summary(ppls_astrag_lmf_tf))
sink()


## ---- Integration between lateral malleolar facet and fhlg ----

# two block pls
pls_astrag_lmf_fhlg <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_lat_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_fhlg,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-lmf-fhlg.txt")
print(summary(pls_astrag_lmf_fhlg))
sink()


# phylogenetic partial least squares
ppls_astrag_lmf_fhlg <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_lat_mal,,],
    A2 = shape_data_means_astrag[lm_astrag_fhlg,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-lmf-fhlg.txt")
print(summary(ppls_astrag_lmf_fhlg))
sink()


## ---- Integration between ectal and navicular facets ----

# two block pls
pls_astrag_ef_nf <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_ectal,,],
    A2 = shape_data_means_astrag[lm_astrag_navic,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-ef-nf.txt")
print(summary(pls_astrag_ef_nf))
sink()


# phylogenetic partial least squares
ppls_astrag_ef_nf <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_ectal,,],
    A2 = shape_data_means_astrag[lm_astrag_navic,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-ef-nf.txt")
print(summary(ppls_astrag_ef_nf))
sink()


## ---- Integration between ectal and sustentacular facets ----

# two block pls
pls_astrag_ef_sf <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_ectal,,],
    A2 = shape_data_means_astrag[lm_astrag_sust,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-ef-sf.txt")
print(summary(pls_astrag_ef_sf))
sink()


# phylogenetic partial least squares
ppls_astrag_ef_sf <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_ectal,,],
    A2 = shape_data_means_astrag[lm_astrag_sust,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-ef-sf.txt")
print(summary(ppls_astrag_ef_sf))
sink()


## ---- Integration between ectal and trochlar facets ----

# two block pls
pls_astrag_ef_tf <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_ectal,,],
    A2 = shape_data_means_astrag[lm_astrag_troc,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-ef-tf.txt")
print(summary(pls_astrag_ef_tf))
sink()


# phylogenetic partial least squares
ppls_astrag_ef_tf <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_ectal,,],
    A2 = shape_data_means_astrag[lm_astrag_troc,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-ef-tf.txt")
print(summary(ppls_astrag_ef_tf))
sink()


## ---- Integration between ectal facet and fhlg ----

# two block pls
pls_astrag_ef_fhlg <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_ectal,,],
    A2 = shape_data_means_astrag[lm_astrag_fhlg,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-ef-fhlg.txt")
print(summary(pls_astrag_ef_fhlg))
sink()


# phylogenetic partial least squares
ppls_astrag_ef_fhlg <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_ectal,,],
    A2 = shape_data_means_astrag[lm_astrag_fhlg,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-ef-fhlg.txt")
print(summary(ppls_astrag_ef_fhlg))
sink()


## ---- Integration between navicular and sustentacular facets ----

# two block pls
pls_astrag_nf_sf <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_navic,,],
    A2 = shape_data_means_astrag[lm_astrag_sust,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-nf-sf.txt")
print(summary(pls_astrag_nf_sf))
sink()


# phylogenetic partial least squares
ppls_astrag_nf_sf <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_navic,,],
    A2 = shape_data_means_astrag[lm_astrag_sust,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-nf-sf.txt")
print(summary(ppls_astrag_nf_sf))
sink()


## ---- Integration between navicular and trochlear facets ----

# two block pls
pls_astrag_nf_tf <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_navic,,],
    A2 = shape_data_means_astrag[lm_astrag_troc,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-nf-tf.txt")
print(summary(pls_astrag_nf_tf))
sink()


# phylogenetic partial least squares
ppls_astrag_nf_tf <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_navic,,],
    A2 = shape_data_means_astrag[lm_astrag_troc,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-nf-tf.txt")
print(summary(ppls_astrag_nf_tf))
sink()


## ---- Integration between navicular facet and fhlg ----

# two block pls
pls_astrag_nf_fhlg <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_navic,,],
    A2 = shape_data_means_astrag[lm_astrag_fhlg,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-nf-fhlg.txt")
print(summary(pls_astrag_nf_fhlg))
sink()


# phylogenetic partial least squares
ppls_astrag_nf_fhlg <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_navic,,],
    A2 = shape_data_means_astrag[lm_astrag_fhlg,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-nf-fhlg.txt")
print(summary(ppls_astrag_nf_fhlg))
sink()


## ---- Integration between sustentacular and trochlear facets ----

# two block pls
pls_astrag_sf_tf <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_sust,,],
    A2 = shape_data_means_astrag[lm_astrag_troc,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-sf-tf.txt")
print(summary(pls_astrag_sf_tf))
sink()


# phylogenetic partial least squares
ppls_astrag_sf_tf <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_sust,,],
    A2 = shape_data_means_astrag[lm_astrag_troc,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-sf-tf.txt")
print(summary(ppls_astrag_sf_tf))
sink()


## ---- Integration between sustentacular facet and fhlg ----

# two block pls
pls_astrag_sf_fhlg <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_sust,,],
    A2 = shape_data_means_astrag[lm_astrag_fhlg,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-sf-fhlg.txt")
print(summary(pls_astrag_sf_fhlg))
sink()


# phylogenetic partial least squares
ppls_astrag_sf_fhlg <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_sust,,],
    A2 = shape_data_means_astrag[lm_astrag_fhlg,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-sf-fhlg.txt")
print(summary(ppls_astrag_sf_fhlg))
sink()


## ---- Integration between trochlear facet and fhlg ----

# two block pls
pls_astrag_tf_fhlg <- 
  integration.test(
    shape_data_means_astrag[lm_astrag_troc,,],
    A2 = shape_data_means_astrag[lm_astrag_fhlg,,],
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/pls/pls-astragalus-tf-fhlg.txt")
print(summary(pls_astrag_tf_fhlg))
sink()


# phylogenetic partial least squares
ppls_astrag_tf_fhlg <- 
  phylo.integration(
    shape_data_means_astrag[lm_astrag_troc,,],
    A2 = shape_data_means_astrag[lm_astrag_fhlg,,],
    phy = phylo,
    partition.gp = NULL,
    seed = NULL,
    print.progress = TRUE
  )

# write results to .txt file
sink("./analysis/integration-within-elements/astragalus/ppls/ppls-astragalus-tf-fhlg.txt")
print(summary(ppls_astrag_tf_fhlg))
sink()


## ---- Global integration ----

glob_int_astrag <- 
  globalIntegration(
    shape_data_means_astrag,
    ShowPlot = TRUE
  )

