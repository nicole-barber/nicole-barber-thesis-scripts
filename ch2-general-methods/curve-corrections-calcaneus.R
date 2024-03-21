
## ---- Info ----

## script for
## adjusting curves
## for variable facets

## requires 
## resampled and slid LMs as a 3D array
## i.e. dataslide output of Morpho::slider3d
## table specifying curves that need adjusting in .csv format

## this script
## is for the calcaneus


## ---- Setup ----

# import table specifying curves to be adjusted
curveadjust_calc <-
  read.csv("./data/calcaneus/import/calcaneus-curve-adjust.csv")


## ---- Define curves to be corrected ----

# define which landmarks belong to curves to be moved
# using curve table
curve6_calc <- 
  my_curves_calc$Curve.in[[6]]

curve8_calc <- 
  my_curves_calc$Curve.in[[8]]

curve12_calc <- 
  my_curves_calc$Curve.in[[12]]


## ---- Correct curves ----

# use a for loop
# to replace landmarks for each curve
# using a matrix of same length as list of curve LMs
# but with coordinates of relevant landmark


## curve 6

# for each row in 'curveadjust' table
for (i in 1:nrow(curveadjust_calc)) {
  # if is.na is NOT true for selected curve 
  # ('!' operator reverses logical function)
  if(!is.na(curveadjust_calc$Curve6[i]))
    # move all landmarks for curve 6 in that specimen
    slidLMs_calc[curve6_calc,,i] <- 
      # to position of landmark 7
      # using a matrix of same length as list of curve LMs
      # with correct x, y and z coords (i.e. that match L7 coords)
      matrix(
        # [landmark 7, all dimensions, specimen number]
        slidLMs_calc[7,,i],
        nrow = length(curve6_calc),
        ncol = 3,
        byrow = TRUE
      )
}


## curve 8

# for each row in 'curveadjust' table
for (i in 1:nrow(curveadjust_calc)) {
  # if is.na is NOT true for selected curve 
  # ('!' operator reverses logical function)
  if(!is.na(curveadjust_calc$Curve8[i]))
    # move all landmarks for curve 8 in that specimen
    slidLMs_calc[curve8_calc,,i] <- 
      # to position of landmark 6
      # using a matrix of same length as list of curve LMs
      # with correct x, y and z coords (i.e. that match L5 coords)
      matrix(
        # [landmark 12, all dimensions, specimen number]
        slidLMs_calc[6,,i],
        nrow = length(curve8_calc),
        ncol = 3,
        byrow = TRUE
      )
}


## curve 12

# for each row in 'curveadjust' table
for (i in 1:nrow(curveadjust_calc)) {
  # if is.na is NOT true for selected curve 
  # ('!' operator reverses logical function)
  if(!is.na(curveadjust_calc$Curve12[i]))
    # move all landmarks for curve 12 in that specimen
    slidLMs_calc[curve12_calc,,i] <- 
      # to position of landmark 10
      # using a matrix of same length as list of curve LMs
      # with correct x, y and z coords (i.e. that match L7 coords)
      matrix(
        # [landmark 7, all dimensions, specimen number]
        slidLMs_calc[10,,i],
        nrow = length(curve12_calc),
        ncol = 3,
        byrow = TRUE
      )
}

