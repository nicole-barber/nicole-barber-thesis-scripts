
## ---- Info ----

## script for
## adjusting curves
## for variable facets

## requires 
## resampled and slid LMs as a 3D array
## i.e. dataslide output of Morpho::slider3d
## table specifying curves that need adjusting in .csv format

## this script
## is for the astragalus


## ---- Setup ----

# import table specifying curves to be adjusted
curveadjust_astrag <-
  read.csv("./data/astragalus/import/astragalus-curve-adjust.csv")


## ---- Define curves to be corrected ----

# define which landmarks belong to curves to be moved
# using curve table
curve10_astrag <- 
  my_curves_astrag$Curve.in[[10]]

curve12_astrag <- 
  my_curves_astrag$Curve.in[[12]]

curve13_astrag <- 
  my_curves_astrag$Curve.in[[13]]


## ---- Correct curves ----

# use a for loop
# to replace landmarks for each curve
# using a matrix of same length as list of curve LMs
# but with coordinates of relevant landmark


## curve 10

# for each row in 'curveadjust' table
for (i in 1:nrow(curveadjust_astrag)) {
  # if is.na is NOT true for selected curve 
  # ('!' operator reverses logical function)
  if(!is.na(curveadjust_astrag$Curve10[i]))
    # move all landmarks for curve 10 in that specimen
    slidLMs_astrag[curve10_astrag,,i] <- 
      # to position of landmark 10
      # using a matrix of same length as list of curve LMs
      # with correct x, y and z coords (i.e. that match L10 coords)
      matrix(
        # [landmark 10, all dimensions, specimen number]
        slidLMs_astrag[10,,i],
        nrow = length(curve10_astrag),
        ncol = 3,
        byrow = TRUE
      )
}


## curve 12

# for each row in 'curveadjust' table
for (i in 1:nrow(curveadjust_astrag)) {
  # if is.na is NOT true for selected curve 
  # ('!' operator reverses logical function)
  if(!is.na(curveadjust_astrag$Curve12[i]))
    # move all landmarks for curve 12 in that specimen
    slidLMs_astrag[curve12_astrag,,i] <- 
      # to position of landmark 5
      # using a matrix of same length as list of curve LMs
      # with correct x, y and z coords (i.e. that match L5 coords)
      matrix(
        # [landmark 12, all dimensions, specimen number]
        slidLMs_astrag[5,,i],
        nrow = length(curve12_astrag),
        ncol = 3,
        byrow = TRUE
      )
}


## curve 13

# for each row in 'curveadjust' table
for (i in 1:nrow(curveadjust_astrag)) {
  # if is.na is NOT true for selected curve 
  # ('!' operator reverses logical function)
  if(!is.na(curveadjust_astrag$Curve13[i]))
    # move all landmarks for curve 13 in that specimen
    slidLMs_astrag[curve13_astrag,,i] <- 
      # to position of landmark 7
      # using a matrix of same length as list of curve LMs
      # with correct x, y and z coords (i.e. that match L7 coords)
      matrix(
        # [landmark 7, all dimensions, specimen number]
        slidLMs_astrag[7,,i],
        nrow = length(curve13_astrag),
        ncol = 3,
        byrow = TRUE
      )
}

