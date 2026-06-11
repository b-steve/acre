# Data to showcase a "mul_ses" demo

For the demonstration of the model with bearing and distance as extra
information, and half normal as detection function. The survey is
carried with 2 sessions.

## Usage

``` r
mul_ses
```

## Format

a list with all necessary input for the model:

- detfn:

  detechtion function: "hn" - half normal

- traps:

  a list with a data frame with the coordinates of the acoustic
  detectors as each of its element

- control.mask:

  a list with the basic argument "buffer", which is the max detectable
  distance, to create the detectable area from the coordinates of the
  detectors

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
