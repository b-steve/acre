# Data to showcase a "ind_ss_log" demo

For the demonstration of the model with signal strength as extra
information and detection function with the link function of log. The
survey is carried with 3 sessions, and each records in the capture
history contains identity of the detected individual.

## Usage

``` r
ind_ss_log
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
  detection, and "ss.link", which indicates the link function for the
  signal strength detection function.

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
