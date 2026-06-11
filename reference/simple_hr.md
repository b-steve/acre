# Data to showcase a "simple_hr" demo

For the demonstration of the model with hazard rate detection function.

## Usage

``` r
simple_hr
```

## Format

a list with all necessary input for the model:

- detfn:

  detechtion function: "hr" - hazard rate

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

- sv:

  a list, contains the coefficient which been assigned a start value for
  modeling

## Source

created from the simulation
