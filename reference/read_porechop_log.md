# Read useful information from a Porechop log file

Read useful information from a Porechop log file

## Usage

``` r
read_porechop_log(file)
```

## Arguments

- file:

  Path to Porechop (ABI) log file

## Value

A list containing:

- `file`: Character. Name of input file provided to Porechop

- `reads_total`: Integer. Total number of reads processed

- `reads_trimmed_start`: Integer. Number of reads with adapters trimmed
  from start

- `reads_trimmed_end`: Integer. Number of reads with adapters trimmed
  from end

- `reads_split`: Integer. Number of reads split based on middle adapters

- `adapter_matches`: Tibble with columns: adapter_set, best_start_perc,
  best_end_perc, selected

- `adapter_sequences`: Tibble with columns: sequence_id, sequence
