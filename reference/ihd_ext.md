# Data to showcase a "ihd_ext" demo

For the demonstration of the model with inhomogeneous density surface,
and half normal as detection function. One of the coefficients of the
detection function, "sigma", is modeled with additional covariates

## Usage

``` r
ihd_ext
```

## Format

a list with all necessary input for the model:

- detfn:

  detechtion function: "hn" - half normal

- model:

  a list with the name of coefficients to be modeled with additional
  covariates as the elements names, and each element contains a formula
  for that modeled coefficient.

- traps:

  a data frame with the coordinates of the acoustic detectors

- control.mask:

  a list with the basic argument "buffer", which is the max detectable
  distance, to create the detectable area from the coordinates of the
  detectors

- trap.cov:

  a data frame contains detectors related covariates for the modeled
  coefficients, it must contains the indices of detectors, if the
  detectors are different between any survey sessions, the indices of
  sessions should be provided as well

- loc.cov:

  a data frame contains location related covariates for the modeled
  coefficients, it must contains column 'x' and 'y' as coordinates.

- capt:

  a data frame of the capture history, contains columns as follows:  
    
  "session" - the indices of the survey sessions;  
  "ID" - the indices of the calls been detected;  
  "occasion" - not been used, could be ignored;  
  "trap" - the indices of the detectors.  

- fix:

  a list, contains the coefficient which been fixed instead of estimated
  by the model

## Source

created from the simulation
