
## ---- Info ----

## script for
## adjusting landmarks
## that belong to variably present regions

## requires 
## resampled and slid LMs as a 3D array
## i.e. dataslide output of Morpho::slider3d
## table specifying variably present regions in .csv format

## this script
## is for the cuboid


## ---- Setup ----

# read in table specifying absent regions
absent_cuboid <- 
  read.csv("./data/cuboid/import/cuboid-absent.csv")


## ---- Define variable regions ----

## define which landmarks belong to variably present facets
## n.b. cannot do in exactly the same way as with variable curves
## because you need the fixed as well as curve LMs

# extract fixed and curve landmarks 
# which belong to dist. medial facet from curve table 
# i.e. match dist_medial in Facet column
lm_dist_medial_cuboid <- 
  my_curves_cuboid$Curves[which(curve_table_cuboid$Facet%in%"dist_medial")]%>%
  # simplify extracted list of curves 
  # to produce single vector 
  # containing all LMs belonging to dist. medial facet
  unlist(.)%>%
  # remove duplicate LM values
  unique(.)%>%
  # sort list of LMs into ascending order
  sort(.)


## ---- Move LMs for variable regions ----

## use a for loop to search for specimens with absent distal medial facet
## and move associated landmarks for that facet to a single point

# for each row in 'absent' table
for (i in 1:nrow(absent_cuboid)) {
  # if is.na is NOT true for selected facet 
  # ('!' operator reverses logical function)
  if(!is.na(absent_cuboid$Dist_Medial[i]))
    # move all landmarks for distal medial facet in that specimen
    slidLMs_cuboid[lm_dist_medial_cuboid,,i] <- 
      # to position of landmark 13
      # using a matrix of same length as list of facet LMs
      # with correct x, y and z coords (i.e. that match L13 coords)
      matrix(
        # [landmark 13, all dimensions, specimen number]
        slidLMs_cuboid[13,,i],
        nrow = length(lm_dist_medial_cuboid),
        ncol = 3,
        byrow = TRUE
      )
}
