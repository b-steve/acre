# Data to showcase a "ind_toa_hhn" demo

For the demonstration of the model with inhomogeneous density surface,
time of arrival as extra information, and hazard half normal as
detection function. The survey is carried with 2 sessions, and each
records in the capture history contains identity of the detected
individual.

## Usage

``` r
ind_toa_hhn
```

## Format

a list with all necessary input for the model:

- detfn:

  detechtion function: "hhn" - hazard half normal

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

- loc.cov:

  a data frame contains location related covariates for the modeled
  coefficients, it must contains column 'x' and 'y' as coordinates.

- control_create_cpat:

  a list with arguments for the function "create.capt()" which converts
  the data input to the data to be fed into the model.

- capt:

  a data frame of the capture history, contains columns as follows:  
    
  "session" - the indices of the survey sessions;  
  "ID" - the indices of the calls been detected;  
  "occasion" - not been used, could be ignored;  
  "trap" - the indices of the detectors;  
  "animal_ID" - the indices or identity of the detected individual;  
  "toa" - the time of arrival.  

## Source

created from the simulation
