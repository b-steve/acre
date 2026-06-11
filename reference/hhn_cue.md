# Data to showcase a "hhn_cue" demo

For the demonstration of the model with hazard half normal detection
function, and cue rate as additional information to convert the call
density to individuals density.

## Usage

``` r
hhn_cue
```

## Format

a list with all necessary input for the model:

- detfn:

  detechtion function: "hhn" - hazard half normal

- traps:

  a data frame with the coordinates of the acoustic detectors

- control.mask:

  a list with the basic argument "buffer", which is the max detectable
  distance, to create the detectable area from the coordinates of the
  detectors

- survey.length:

  a numeric vector or a scalar contains the length of each session.

- cue.rates:

  a numeric vector. contains the recorded cue rates in a series of time
  periods with identical length. The length should be equal to the unit
  length in the "survey.length".

- capt:

  a data frame of the capture history, contains columns as follows:  
    
  "session" - the indices of the survey sessions;  
  "ID" - the indices of the calls been detected;  
  "occasion" - not been used, could be ignored;  
  "trap" - the indices of the detectors.  

## Source

created from the simulation
