# Complete a set of Kraken reports by filling in zero counts (Not tested!)

This function uses
[`tidyr::complete()`](https://tidyr.tidyverse.org/reference/complete.html)
to add rows with zeros for taxa that missing in some samples. Note: The
resulting data frame will no longer be compatible with
[`kraken_add_taxonomy()`](https://securebio.github.io/naoutils-rpkg/reference/kraken_add_taxonomy.md).

## Usage

``` r
kraken_complete(x, explicit = TRUE)
```

## Arguments

- x:

  A data frame of Kraken sample reports

- explicit:

  Boolean passed to
  [`tidyr::complete()`](https://tidyr.tidyverse.org/reference/complete.html)
