
## ---- Info ----

## script for
## reading in raw landmarks
## resampling curves
## sliding curves
## adjusting curves n.b. source separate script
## performing GPA
## extracting species mean shapes
## extracting species mean csize

## requires
## raw landmarks for each specimen in .pts format (IDAV Landmark Editor)
## mesh for each specimen in .ply format (ASCII)
## table specifying curve LMs in .csv format
## table with specimen data in .csv format
## table with species data in .csv format

## this script
## is for the astragalus


## ---- Setup ----

# if SURGE not installed already, run
# install.packages("devtools")
# library("devtools")
# devtools::install_github("rnfelice/SURGE")

# load packages
library(tidyverse)
library(SURGE)
library(geomorph)
library(Morpho)
library(Rvcg)
library(rgl)
library(abind)

# source custom functions
source("./scripts/custom-functions.R")

# read in specimen data
specimen_data <- 
  read_csv("./data/all-elements/import/specimen-data.csv")

# read in species data
species_data <- 
  read_csv("./data/all-elements/import/species-data.csv")


## ---- Locate and define landmarks and curves ----

## load curve table and specify pts files

# import table defining your curves
curve_table_astrag <- 
  read_csv("./data/astragalus/import/astragalus-curve-table.csv")

# identify the folder where your pts files are
ptsfolder_astrag <- 
  "./data/astragalus/import/pts"

# import the pts file names
ptslist_astrag <- 
  dir(
    ptsfolder_astrag, 
    pattern = ".pts", 
    recursive = F
  )

# use curve table to define curves in Morpho format
my_curves_astrag <- 
  create_curve_info(
    curve_table_astrag,
    # specify number of fixed landmarks
    n_fixed = 17
  )

# save and export my_curves output for later
save(
  my_curves_astrag, 
  file = "./data/astragalus/output/astragalus_curve-info.R"
)


## ---- Read in landmarks, resample curves ----

# set working directory to the folder where your pts files are
setwd(ptsfolder_astrag)

# import pts files and 
# resample curves to have the same number of points
subsampled_lm_astrag <- 
  import_chkpt_data(
    ptslist_astrag, 
    my_curves_astrag, 
    subsampl = TRUE
  )

# reset working directory to project folder 
# (by moving up four directories from your pts folder)
setwd("./../../../..")

# set path to ply files
plypath_astrag <- 
  "./data/astragalus/import/ply/"

# check to make sure subsampled curves look okay on each specimen
# using modified checkLM function that scales landmarks to centroid size
checkLM.mod(
  subsampled_lm_astrag,
  begin = 1,
  pt.size = NULL,
  pt.size.scale = 0.001,
  path = plypath_astrag,
  suffix = ".ply",
  render = "s"
)


## --- Slide curves ----

# slide landmarks
slid_astrag <- 
  slider3d(
    subsampled_lm_astrag, 
    SMvector = my_curves_astrag$Sliding.LMs,
    outlines = my_curves_astrag$Curves, 
    sur.path = plypath_astrag, 
    meshlist = paste(
      "./data/astragalus/import/ply/", 
      dimnames(subsampled_lm_astrag)[[3]], 
      ".ply",sep = ""
    ), 
    ignore = NULL,
    sur.type = "ply", 
    tol = 1e-10, 
    deselect = FALSE, 
    inc.check = FALSE,
    recursive = TRUE, 
    iterations = 3, 
    initproc = TRUE,
    pairedLM = 0, 
    bending = TRUE,
    fixRepro = FALSE,
    stepsize = 0.1
  )


# check what slid curves look like on each specimen
# and make a note of any curves
# you will need to adjust later
checkLM.mod(
  slid_astrag$dataslide,
  begin = 1,
  pt.size = NULL,
  pt.size.scale = 0.001, 
  path = plypath_astrag,
  suffix = ".ply",
  render = "s"
)


# assign slid landmarks to a new object
slidLMs_astrag <- 
  slid_astrag$dataslide


## ---- Adjust curves ----

# correct any curves that have slid out of place 
# or that need adjusting
source("./scripts/curve-corrections-astragalus.R")


# check corrected curve landmarks
checkLM.mod(
  slidLMs_astrag,
  begin = 1,
  pt.size = NULL,
  pt.size.scale = 0.001, 
  path = plypath_astrag,
  suffix = ".ply",
  render = "s"
)


## ---- Individual GPA ----

# perform GPA on all specimens
gpa_individual_astrag <- 
  gpagen(slidLMs_astrag)


# adjust curves again if needed
# source("./scripts/curve-correction-astragalus.R")


# save and export gpa output
save(
  gpa_individual_astrag, 
  file = "./data/astragalus/output/astragalus_individual-gpa-output.R"
)


# extract shape data from gpa
shape_data_individual_astrag <- 
  gpa_individual_astrag$coords


# save and export gpa adjusted landmarks
# as R object
save(
  shape_data_individual_astrag, 
  file = "./data/astragalus/output/astragalus_individual-gpa-coords.R"
)

# as csv
write.csv(
  x = two.d.array(shape_data_individual_astrag),
  file = "./data/astragalus/output/astragalus_individual-gpa-coords.csv"
)


## extract size data (cs) from gpa
csize_individual_astrag <-
  gpa_individual_astrag$Csize


# save and export individual cs
# as R object
save(
  csize_individual_astrag,
  file = "./data/astragalus/output/astragalus_individual-csize.R"
)

# as csv
write.csv(
  x = csize_individual_astrag,
  file = "./data/astragalus/output/astragalus_individual-csize.csv"
)


## ---- Extract species mean shapes ----

# extract species ID for each specimen
species_factor <-
  # convert to a factor
  as.factor(
    # pull out relevant column of specimen data table
    specimen_data$Species_ID
  )


# check levels for your species factor
# should be same as number of species you have
levels(species_factor)


# find species means
species_means_astrag <-
  # apply over
  lapply(
    # the list from 1:number of levels in species_factor object i.e. 52
    1:length(levels(species_factor)),
    # the following function
    function(x) 
      # find mean shape
      mshape(
        # of individual gpa coords
        shape_data_individual_astrag[,,which(
          # for everything in each level of species_factor object
          species_factor==levels(species_factor)[x]
        )]
      )
  )


# rename species means
names(species_means_astrag) <-
  # to match levels of species factor 
  levels(species_factor)


# extract and convert species means to a 3D array
shape_data_means_astrag <-
  # combine arrays
  abind(
    species_means_astrag,
    # along third dimension i.e. number of configurations/specimens
    along = 3
  )


# check the dimensions of the final object
# to make sure number of landmarks and species is correct 
dim(shape_data_means_astrag)


# save and export species mean shapes
# as R object
save(
  shape_data_means_astrag, 
  file = "./data/astragalus/output/astragalus_species-gpa-coords.R"
)

# as csv
write.csv(
  x = two.d.array(shape_data_means_astrag),
  file = "./data/astragalus/output/astragalus_species-gpa-coords.csv"
)


## ---- Extract species mean csize ----

# find species means
csize_means_astrag <-
  # apply over
  lapply(
    # the list from 1:number of levels in species_factor object i.e. 52
    1:length(levels(species_factor)),
    # the following function
    function(x) 
      # find mean centroid size
      mean(
        # of individual centroid sizes
        csize_individual_astrag[which(
          # for everything in each level of species_factor object
          species_factor==levels(species_factor)[x]
        )]
      )
  )


# reformat mean centroid sizes to numeric object
# since it is currently a list
csize_means_astrag <- 
  as.numeric(csize_means_astrag)

# check
is.numeric(csize_means_astrag)


# rename cs means
names(csize_means_astrag) <-
  # to match levels of species factor 
  levels(species_factor)


# save and export species mean centroid size
# as R object
save(
  csize_means_astrag, 
  file = "./data/astragalus/output/astragalus_species-csize.R"
)

# as csv
write.csv(
  csize_means_astrag,
  file = "./data/astragalus/output/astragalus_species-csize.csv"
)
