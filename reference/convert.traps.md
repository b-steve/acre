# Convert traps object

Converts an `acre` traps matrix to a `secr` traps object.

## Usage

``` r
convert.traps(traps, ss = FALSE)
```

## Arguments

- traps:

  a matrix or a data frame, contains one session's detectors'
  coordinates

- ss:

  Logical, set to `TRUE` if a signal strength detection function is to
  be used.

## Value

An object of class `traps` comprising a data frame of x- and
y-coordinates, the detector type ('single', 'multi', 'proximity',
'count', 'polygon' etc.), and possibly other attributes.

## Details

The returned object is suitable for use as the `traps` argument of the
function make.capthist.
