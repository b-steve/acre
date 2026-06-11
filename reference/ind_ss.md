# Data to showcase a "ind_ss" demo

For the demonstration of the model with inhomogeneous density surface,
signal strength as extra information and detection function. One of the
coefficients of detection function, "b0.ss", is modeled by additional
covariates. The survey is carried with 3 sessions, and each records in
the capture history contains identity of the detected individual.

## Usage

``` r
ind_ss
```

## Format

a list with all necessary input for the model:

- detfn:

  detechtion function: "ss" - signal strength

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

- ss.opts:

  a list contains signal strength model related options. Here it
  contains "cutoff", which indicates the threshold of a distance of 100%
  detection.

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
  "ss" - the detected signal strength.  

## Source

created from the simulation
