# Read one or more Kraken2 sample reports

This function can read one or more Kraken(2) sample reports, as well as
output of `kraken2-inspect`. Kraken2 reports with and without minimizer
data are supported.

## Usage

``` r
read_kraken_reports(files, minimizers, id = "file", num_threads = 1)
```

## Arguments

- files:

  Paths to files

- minimizers:

  Boolean specifying whether the files include minimizer data

- id:

  A string or NULL. If a string, create column with that name for
  storing the files or the names of the `files` vector

- num_threads:

  Parameter passed onto
  [`readr::read_tsv()`](https://readr.tidyverse.org/reference/read_delim.html)

## Details

The column names in the resulting data frame are
"clade_fragments_percent", "clade_fragments", "node_fragments",
"minimizers_total", "minimizers_distinct", "rank_code", "taxid",
"scientific_name"; the minimizer columns are only included if
`minimizers = TRUE`.

This function should work on Kraken and Kraken2 reports as well as on
the output of `kraken2-inspect`, but we've not tested it on original
Kraken reports. Currently there is no check as to whether the minimizer
data is in the reports or not; if `minimizers` is set incorrectly, the
function will return without error but the reports will be parsed
incorrectly.

- [Kraken2
  docs](https://github.com/DerrickWood/kraken2/wiki/Manual#sample-report-output-format)

- [Kraken
  docs](https://ccb.jhu.edu/software/kraken/MANUAL.html#output-format)
