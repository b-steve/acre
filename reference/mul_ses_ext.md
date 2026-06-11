# Data to showcase a "mul_ses_ext" demo

For the demonstration of the model with bearing and distance as extra
information, and half normal as detection function. Two of the
coefficients of the detection function, "g0" and "sigma", are modeled
with additional covariates. The survey is carried with 2 sessions.

## Usage

``` r
mul_ses_ext
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

  a list with a data frame with the coordinates of the acoustic
  detectors as each of its element

- control.mask:

  a list with the basic argument "buffer", which is the max detectable
  distance, to create the detectable area from the coordinates of the
  detectors

- session.cov:

  a data frame contains survey sessions related covariates for the
  modeled coefficients, it must contains the indices of sessions

- trap.cov:

  a data frame contains detectors related covariates for the modeled
  coefficients, it must contains the indices of detectors, if the
  detectors are different between any survey sessions, the indices of
  sessions should be provided as well

- capt:

  a data frame of the capture history, contains columns as follows:  
    
  "session" - the indices of the survey sessions;  
  "ID" - the indices of the calls been detected;  
  "occasion" - not been used, could be ignored;  
  "trap" - the indices of the detectors;  
  "bearing" - extra information, the direction of the location where a
  call is detected;  
  "dist" - extra information, the distance to and the location where a
  call is detected.  

- sv:

  a list, contains the coefficient which been assigned a start value for
  modeling

## Source

created from the simulation
