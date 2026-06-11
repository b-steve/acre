# Data to showcase a "simple_hhn" demo

For the demonstration of the model with hazard half normal detection
function.

## Usage

``` r
simple_hhn
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

- capt:

  a data frame of the capture history, contains columns as follows:  
    
  "session" - the indices of the survey sessions;  
  "ID" - the indices of the calls been detected;  
  "occasion" - not been used, could be ignored;  
  "trap" - the indices of the detectors.  

## Source

created from the simulation
