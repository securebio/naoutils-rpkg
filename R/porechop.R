# Functions for working with Porechop and Porechop ABI log files

#' Read useful information from a Porechop log file
#'
#' @param file Path to Porechop (ABI) log file
#'
#' @return A list containing:
#'   * `reads_total`: Integer. Total number of reads processed
#'   * `reads_trimmed_start`: Integer. Number of reads with adapters trimmed from start
#'   * `reads_trimmed_end`: Integer. Number of reads with adapters trimmed from end
#'   * `reads_split`: Integer. Number of reads split based on middle adapters
#'   * `adapter_matches`: Tibble with columns: adapter_set, best_start_perc, best_end_perc, selected
#'   * `adapter_sequences`: Tibble with columns: sequence_id, sequence
#' @export
read_porechop_log <- function(file) {
  # Input validation
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  # Helper function to extract count stats with str_match and format as numbers
  extract_stats <- function(lines, pattern) {
    m <- str_match(lines, pattern)
    m[!is.na(m[, 1]), -1] |>
      map_chr(\(x) str_replace(x, ",", "")) |>
      map_int(as.numeric)
  }

  # The raw porechop file is hard to work with because of a mix of newline and
  # carriage returns and ANSI color codes. We'll clean it up in two steps,
  # pausing to use one of the codes to extract the selected adapters from the
  # matching stats table.
  lines <-  file |>
    readr::read_file() |>
    str_replace_all("\r", "\n") |>
    I() |>
    readr::read_lines() |>
    str_subset("^$", negate = TRUE) |>
    trimws()
  selected_adapter_lines <- lines |>
    str_subset("\\033\\[32m") |>
    cli::ansi_strip()
  selected_adapter_sets <- selected_adapter_lines |>
    str_split_i(" {2,}", 1)
  lines <- lines |>
    cli::ansi_strip()

  # The table with sequence identities for known adapter sets can be found from
  # the header line up to the trimming message
  idx_start <- str_which(
    lines,
    "Set +%ID +%ID"
  )
  idx_end <- str_which(
    lines,
    "Trimming adapters from read ends"
  )
  adapter_matches <- tibble(
    line = lines[seq(idx_start + 1, idx_end - 1)]
  ) |>
    separate(line,
      into = c("adapter_set", "best_start_perc", "best_end_perc"),
      sep = " {2,}"
    ) |>
    mutate(
      across(starts_with("best_"), as.numeric),
      selected = adapter_set %in% selected_adapter_sets
    )

  # Adapter sequences used for trimming are between the message indicating that
  # adapter trimming is starting and the first progress indicator
  idx_start <- str_which(
    lines,
    "Trimming adapters from read ends"
  )
  idx_delta <- str_which(
    tail(lines, -idx_start),
    "[[:digit:],]+ / [[:digit:],]+ \\("
  )[1]
  trimming_sequences <- tibble(
    line = lines[seq(idx_start + 1, idx_start + idx_delta - 1)]
  ) |>
    separate(line, into = c("sequence_id", "sequence"), sep = ": +")

  # Line after "Loading reads" has the sample file name
  sample_file <- lines[str_which(lines, "Loading reads") + 1]

  # Number of reads processed: Parse from line like "16,174 reads loaded"
  reads_total <- extract_stats(lines,
    "(?<total>[[:digit:],]+) reads loaded"
  )[1]

  # Reads with adapters trimmed from start
  # Example line: "2,193 / 16,174 reads had adapters trimmed from their start
  # (22,313 bp removed)"
  stats_start <- extract_stats(lines,
    "(?<trimmed>[[:digit:],]+) / (?<total>[[:digit:],]+) reads had adapters trimmed from their start \\((?<bp>[[:digit:],]+) bp removed\\)"
  )

  # Reads with adapters trimmed from end: Parse similar to start stats
  stats_end <- extract_stats(lines,
    "(?<trimmed>[[:digit:],]+) / (?<total>[[:digit:],]+) reads had adapters trimmed from their end \\((?<bp>[[:digit:],]+) bp removed\\)"
  )

  # Reads split at middle adapters
  # Example line: "44 / 16,174 reads were split based on middle adapters"
  stats_middle <- extract_stats(lines,
    "(?<trimmed>[[:digit:],]+) / (?<total>[[:digit:],]+) reads were split based on middle adapters"
  )

  # All totals should be the same
  totals <- list(
    reads_total,
    stats_start["total"], stats_end["total"], stats_middle["total"]
  )

  if (length(unique(totals)) != 1L) {
    stop("Read totals don't match.")
  }

  results <- list(
    reads_total = reads_total |> unname(),
    reads_trimmed_start = stats_start["trimmed"] |> unname(),
    reads_trimmed_end = stats_end["trimmed"] |> unname(),
    reads_split = stats_middle["trimmed"] |> unname(),
    adapter_matches = adapter_matches,
    adapter_sequences = trimming_sequences
  )

  results
}
