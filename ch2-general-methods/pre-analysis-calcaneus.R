
## ---- Info ----

## script for
## reading in raw landmarks
## resampling curves
## sliding curves
## adjusting curves n.b. source separate script
## performing GPA
## extracting species mean shapes

## requires
## raw landmarks for each specimen in .pts format (IDAV Landmark Editor)
## mesh for each specimen in .ply format (ASCII)
## table specifying curve LMs in .csv format
## table with specimen data in .csv format
## table with species data in .csv format

## this script
## is for the calcaneus


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
curve_table_calc <- 
  read_csv("./data/calcaneus/import/calcaneus-curve-table.csv")

# identify the folder where your pts files are
ptsfolder_calc <- 
  "./data/calcaneus/import/pts"

# import the pts file names
ptslist_calc <- 
  dir(
    ptsfolder_calc, 
    pattern = ".pts", 
    recursive = F
  )

# use curve table to define curves in Morpho format
my_curves_calc <- 
  create_curve_info(
    curve_table_calc,
    # specify number of fixed landmarks
    n_fixed = 17
  )

# save and export my_curves output for later
save(
  my_curves_calc, 
  file = "./data/calcaneus/output/calcaneus_curve-info.R"
)


## ---- Read in landmarks, resample curves ----

# set working directory to the folder where your pts files are
setwd(ptsfolder_calc)

# import pts files and 
# resample curves to have the same number of points
subsampled_lm_calc <- 
  import_chkpt_data(
    ptslist_calc, 
    my_curves_calc, 
    subsampl = TRUE
  )

# reset working directory 
# (by moving up four directories from your pts folder)
setwd("./../../../..")

# set path to ply files
plypath_calc <- 
  "./data/calcaneus/import/ply/"

# check to make sure subsampled curves look okay on each specimen
# using modified checkLM function that scales landmarks to centroid size
checkLM.mod(
  subsampled_lm_calc,
  begin = 1,
  pt.size = NULL,
  pt.size.scale = 0.001,
  path = plypath_calc,
  suffix = ".ply",
  render = "s"
)


## --- Slide curves ----

# slide landmarks
slid_calc <- 
  slider3d(
    subsampled_lm_calc, 
    SMvector = my_curves_calc$Sliding.LMs,
    outlines = my_curves_calc$Curves, 
    sur.path = plypath_calc, 
    meshlist = paste(
      "./data/calcaneus/import/ply/", 
      dimnames(subsampled_lm_calc)[[3]], 
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
  slid_calc$dataslide,
  begin = 1,
  pt.size = NULL,
  pt.size.scale = 0.001, 
  path = plypath_calc,
  suffix = ".ply",
  render = "s"
)


# assign slid landmarks to a new object
slidLMs_calc <- 
  slid_calc$dataslide


## ---- Adjust curves ----

# correct any curves that have slid out of place 
# or that need adjusting
source("./scripts/curve-corrections-calcaneus.R")


# check corrected curve landmarks
checkLM.mod(
  slidLMs_calc,
  begin = 1,
  pt.size = NULL,
  pt.size.scale = 0.001, 
  path = plypath_calc,
  suffix = ".ply",
  render = "s"
)


## ---- Individual GPA ----

# perform GPA on all specimens
gpa_individual_calc <- 
  gpagen(slidLMs_calc)


# save and export gpa output
save(
  gpa_individual_calc, 
  file = "./data/calcaneus/output/calcaneus_individual-gpa-output.R"
)


# extract shape data from gpa
shape_data_individual_calc <- 
  gpa_individual_calc$coords


# adjust curves again if needed
# source("./scripts/curve-correction-calcaneus.R")


# save and export gpa adjusted landmarks
# as R object
save(
  shape_data_individual_calc, 
  file = "./data/calcaneus/output/calcaneus_individual-gpa-coords.R"
)

# as csv
write.csv(
  x = two.d.array(shape_data_individual_calc),
  file = "./data/calcaneus/output/calcaneus_individual-gpa-coords.csv"
)


## extract size data (cs) from gpa
csize_individual_calc <-
  gpa_individual_calc$Csize


# save and export individual cs
# as R object
save(
  csize_individual_calc,
  file = "./data/calcaneus/output/calcaneus_individual-csize.R"
)

# as csv
write.csv(
  x = csize_individual_calc,
  file = "./data/calcaneus/output/calcaneus_individual-csize.csv"
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
species_means_calc <-
  # apply over
  lapply(
    # the list from 1:number of levels in species_factor object i.e. 52
    1:length(levels(species_factor)),
    # the following function
    function(x) 
      # find mean shape
      mshape(
        # of individual gpa coords
        shape_data_individual_calc[,,which(
          # for everything in each level of species_factor object
          species_factor==levels(species_factor)[x]
        )]
      )
  )


# rename species means
names(species_means_calc) <-
  # to match levels of species factor 
  levels(species_factor)


# extract and convert species means to a 3D array
shape_data_means_calc <-
  # combine arrays
  abind(
    species_means_calc,
    # along third dimension i.e. number of configurations
    along = 3
  )


# check the dimensions of the final object
# to make sure number of landmarks and species is correct 
dim(shape_data_means_calc)


# save and export species mean shapes
# as R object
save(
  shape_data_means_calc, 
  file = "./data/calcaneus/output/calcaneus_species-gpa-coords.R"
)

# as csv
write.csv(
  x = two.d.array(shape_data_means_calc),
  file = "./data/calcaneus/output/calcaneus_species-gpa-coords.csv"
)


## ---- Extract species mean csize ----

# find species means
csize_means_calc <-
  # apply over
  lapply(
    # the list from 1:number of levels in species_factor object i.e. 52
    1:length(levels(species_factor)),
    # the following function
    function(x) 
      # find mean centroid size
      mean(
        # of individual centroid sizes
        csize_individual_calc[which(
          # for everything in each level of species_factor object
          species_factor==levels(species_factor)[x]
        )]
      )
  )


# reformat mean centroid sizes to numeric object
# since it is currently a list
csize_means_calc <- 
  as.numeric(csize_means_calc)

# check
is.numeric(csize_means_calc)


# rename cs means
names(csize_means_calc) <-
  # to match levels of species factor 
  levels(species_factor)


# save and export species mean centroid size
# as R object
save(
  csize_means_calc, 
  file = "./data/calcaneus/output/calcaneus_species-csize.R"
)

# as csv
write.csv(
  csize_means_calc,
  file = "./data/calcaneus/output/calcaneus_species-csize.csv"
)

