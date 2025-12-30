# Add full taxonomy path to a Kraken report data frame

This function requires that rows are in the order of the original Kraken
report, so cannot be used after
[`kraken_complete()`](https://securebio.github.io/naoutils-rpkg/reference/kraken_complete.md).

## Usage

``` r
kraken_add_taxonomy(x)
```

## Arguments

- x:

  A data frame containing the fields "rank_code", "scientific_name", and
  "rank_level" with rows in the original order of the Kraken sample
  report(s), as returned by
  [`read_kraken_reports()`](https://securebio.github.io/naoutils-rpkg/reference/read_kraken_reports.md).

## Details

You can apply this function to a data frame of Kraken reports for
multiple samples; however, doing so is typically inefficient due to
needing to recompute paths for taxa that appear in multiple samples.
Instead, consider generating taxonomy paths from the original database
and joining these to the kraken report data frame.
