# A helper function to obtain the distances of the nearest points from a data frame

A helper function to obtain the distances of the nearest points from a
data frame

## Usage

``` r
dist_nearest(from, to, col_name = "dist_nearest")
```

## Arguments

- from:

  a matrix or a data frame with columns "x" and "y" contains the
  coordinates of the start points.

- to:

  a matrix or a data frame with columns "x" and "y" contains the
  coordinates of the end points.

- col_name:

  a character, contains the new column name for the distances in the
  output, default is dist_nearest.

## Value

a data frame as the same as "from", but with an extra column with
assigned column name. For each row, it contains the distance to the
nearest point in the "to" data set.
