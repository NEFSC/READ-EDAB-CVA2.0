# Bhattacharyya Coefficient of two distributions

Calculates the Bhattacharyya Coefficient (BC) between two distributions.
BC is the sum of the square root of the multiplication of the elements
in bin i for both distributions. It ranges from 0 (no overlap) to 1
(full overlap). v.0.1. Used within `eval_brt`. Included to ensure
functionality.

## Usage

``` r
bhatt_coeff_brt(x, y, bw = bw.nrd0, ...)
```

## Source

From Camrin Brawn (WHOI): https://zenodo.org/records/7971532.

## Arguments

- x:

  a numeric vector of length \>= 2.

- y:

  a numeric vector of length \>= 2.

- bw:

  can be either a fixed value of bins or a function to calculate the
  bandwidth of each bin (see ?bw.nrd). Default is bw.nrd0.

- ...:

  any optional arguments to be passed to the given bw function.
