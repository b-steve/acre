# Data to showcase a "ss" demo

For the demonstration of the model with signal strength detection
function, and signal strength as additional information.

## Usage

``` r
ss
```

## Format

a list with all necessary input for the model:

- detfn:

  detechtion function: "ss" - signal strength

- traps:

  a data frame with the coordinates of the acoustic detectors

- control.mask:

  a list with the basic argument "buffer", which is the max detectable
  distance, to create the detectable area from the coordinates of the
  detectors

- ss.opts:

  a list contains signal strength model related options. Here it
  contains "cutoff", which indicates the threshold of a distance of 100%
  detection.

- capt:

  a data frame of the capture history, contains columns as follows:  
    
  "session" - the indices of the survey sessions;  
  "ID" - the indices of the calls been detected;  
  "occasion" - not been used, could be ignored;  
  "trap" - the indices of the detectors;  
  "ss" - the detected signal strength.  

- sv:

  a list, contains the coefficient which been assigned a start value for
  modeling

## Source

created from the simulation
