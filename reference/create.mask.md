# Create mask object

Creates a mask object to use with the function read.acre().

## Usage

``` r
create.mask(traps, buffer, ...)
```

## Arguments

- traps:

  a matrix or a data frame with two columns or a list of such matrices
  or data frames for a multi-session model. Each row in a matrix/data
  frame provides Cartesian coordinates (in metres) for the location of a
  detector. In a list of matrices or data frames, each element of the
  list corresponds to the detector location of a different session. If
  the detector locations stayed the same across several sessions, only
  one matrix/data frame is required.

- buffer:

  a scalar, the furthest distance that a detector could detect (in
  metres)
