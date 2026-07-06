# Introduction to `acre`

Just placeholder-ish for now.

## Installation

If you haven’t already, make sure to install the `acre` using
`devtools`.

``` r

library(devtools)
install_github("b-steve/acre")
```

## Reading in data

The first step is to combine all of these data sources together into an
R object using the read.acre() function. It is your job to create this
object. There are five arguments you’ll need to use:

1.  `captures`: the data frame with the detection data. traps: the
    object with the listening post locations.

2.  `control.mask`: you need to specify the maximum distance at which
    you can possibly detect a gibbon. In this case, the maximum feasible
    distance is 3000m. Use `control.mask = list(buffer = 3000)` for this
    argument.

3.  `loc.cov`: a data frame with columns x and y, specifying locations
    at which spatial covariates have been measured, and then a further
    column for each spatial covariates themselves. For this step you
    should provide the data frame containing the measured covariate
    values, rather than the interpolated values. The function will
    complete the interpolation for you.

4.  `dist.cov`: a data frame containing locations of objects of
    interest, from which you want to construct a spatial covariate for
    the distance to the nearest object. This needs to be a list, where
    each component name relates to the type of object, and the component
    itself is a data frame with columns named x and y specifying the
    locations of these objects. In this case, we just have to obtain the
    distance to the nearest village for each point in PPWS, so you can
    use `dist.cov = list(village = villages.df)`.

## Visualizing data

Multiple exploratory plotting tools available.

See…

1.  `plot(data.acre, type = "survey")`

2.  `plot(data.acre, type = "capt")`

3.  `plot(data.acre, type = "covariates")`

## Fitting models

Maybe fit a basic model here? I think that in this introduction we
should definitely use an IHD example, but maybe one that is as simple as
possible?

Just because being able to plot the detection surface, and seeing how it
varies is super neat, rather than just one uniform-coloured detection
surface blob.

1.  `fit <- fit.acre(data.acre)`

2.  `summary(fit)`

3.  `plot(data.acre, type = "detfn")`

4.  `plot(data.acre, type = "Dsurf")`

## More…(not sure what to title this section)

Here explain all the additional functionality provided, and then suggest
seeing the tutorials for further guidance?
