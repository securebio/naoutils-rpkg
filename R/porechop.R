# Functions for working with Porechop and Porechop ABI log files

#' Read useful information from a Porechop log file
#'
#' @param file Path to Porechop (ABI) log file
#'
#' @return A list containing:
#'   * `file`: Character. Name of input file provided to Porechop
#'   * `reads_total`: Integer. Total number of reads processed
#'   * `reads_trimmed_start`: Integer. Number of reads with adapters trimmed from start
#'   * `reads_trimmed_end`: Integer. Number of reads with adapters trimmed from end
#'   * `reads_split`: Integer. Number of reads split based on middle adapters
#'   * `adapter_matches`: Tibble with columns: adapter_set, best_start_perc, best_end_perc, selected
#'   * `adapter_sequences`: Tibble with columns: sequence_id, sequence
#' @export
read_porechop_log <- function(file) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  PATTERN_READS_LOADED <- "(?<total>[[:digit:],]+) reads loaded"
  PATTERN_START_TRIM <- "(?<trimmed>[[:digit:],]+) / (?<total>[[:digit:],]+) reads had adapters trimmed from their start \\((?<bp>[[:digit:],]+) bp removed\\)"
  PATTERN_END_TRIM <- "(?<trimmed>[[:digit:],]+) / (?<total>[[:digit:],]+) reads had adapters trimmed from their end \\((?<bp>[[:digit:],]+) bp removed\\)"
  PATTERN_MIDDLE_SPLIT <- "(?<trimmed>[[:digit:],]+) / (?<total>[[:digit:],]+) reads were split based on middle adapters"
  PATTERN_ADAPTER_TABLE_HEADER <- "Set +%ID +%ID"
  PATTERN_TRIMMING_SECTION <- "Trimming adapters from read ends"
  PATTERN_PROGRESS_INDICATOR <- "[[:digit:],]+ / [[:digit:],]+ \\("
  PATTERN_NO_ADAPTERS <- "No adapters found"

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

  no_adapters_found <- any(str_detect(lines, PATTERN_NO_ADAPTERS))

  # The table with sequence identities for known adapter sets can be found from
  # the header line up to the trimming message (or "No adapters found" message)
  idx_adapter_table_start <- str_which(lines, PATTERN_ADAPTER_TABLE_HEADER)
  if (length(idx_adapter_table_start) == 0) {
    stop("Could not find adapter matching table header in log file. ",
         "Expected pattern: '", PATTERN_ADAPTER_TABLE_HEADER, "'")
  }

  if (no_adapters_found) {
    idx_adapter_table_end <- str_which(lines, PATTERN_NO_ADAPTERS)
  } else {
    idx_adapter_table_end <- str_which(lines, PATTERN_TRIMMING_SECTION)
    if (length(idx_adapter_table_end) == 0) {
      stop("Could not find adapter trimming section in log file. ",
           "Expected pattern: '", PATTERN_TRIMMING_SECTION, "'")
    }
  }

  adapter_matches <- tibble(
    line = lines[seq(idx_adapter_table_start + 1, idx_adapter_table_end - 1)]
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
  if (no_adapters_found) {
    trimming_sequences <- tibble(
      sequence_id = character(0),
      sequence = character(0)
    )
  } else {
    idx_trimming_start <- str_which(lines, PATTERN_TRIMMING_SECTION)
    # This pattern was already found above, so we know it exists

    idx_progress <- str_which(
      tail(lines, -idx_trimming_start),
      PATTERN_PROGRESS_INDICATOR
    )[1]

    if (is.na(idx_progress)) {
      stop("Could not find progress indicator after trimming section in log file. ",
           "Expected pattern: '", PATTERN_PROGRESS_INDICATOR, "'")
    }

    trimming_sequences <- tibble(
      line = lines[seq(idx_trimming_start + 1, idx_trimming_start + idx_progress - 1)]
    ) |>
      separate(line, into = c("sequence_id", "sequence"), sep = ": +")
  }

  # Line after "Loading reads" has the sample file name
  sample_file <- lines[str_which(lines, "Loading reads") + 1]

  # Number of reads processed: Parse from line like "16,174 reads loaded"
  reads_total <- extract_stats(lines, PATTERN_READS_LOADED)[1]
  if (is.na(reads_total) || length(reads_total) == 0) {
    stop("Could not find total reads count in log file. ",
         "Expected pattern: '", PATTERN_READS_LOADED, "'")
  }

  # Extract trimming statistics (or set to 0 if no adapters found)
  if (no_adapters_found) {
    stats_start <- c(trimmed = 0L, total = reads_total, bp = 0L)
    stats_end <- c(trimmed = 0L, total = reads_total, bp = 0L)
    stats_middle <- c(trimmed = 0L, total = reads_total)
  } else {
    # Reads with adapters trimmed from start
    # Example line: "2,193 / 16,174 reads had adapters trimmed from their start
    # (22,313 bp removed)"
    stats_start <- extract_stats(lines, PATTERN_START_TRIM)
    if (length(stats_start) == 0) {
      stop("Could not find start trimming statistics in log file. ",
           "Expected pattern: '", PATTERN_START_TRIM, "'")
    }

    stats_end <- extract_stats(lines, PATTERN_END_TRIM)
    if (length(stats_end) == 0) {
      stop("Could not find end trimming statistics in log file. ",
           "Expected pattern: '", PATTERN_END_TRIM, "'")
    }

    # Reads split at middle adapters
    # Example line: "44 / 16,174 reads were split based on middle adapters"
    stats_middle <- extract_stats(lines, PATTERN_MIDDLE_SPLIT)
    if (length(stats_middle) == 0) {
      stop("Could not find middle adapter splitting statistics in log file. ",
           "Expected pattern: '", PATTERN_MIDDLE_SPLIT, "'")
    }
  }

  # Validate that all totals match (only when adapters were found)
  if (!no_adapters_found) {
    totals <- list(
      reads_total,
      stats_start["total"], stats_end["total"], stats_middle["total"]
    )

    unique_totals <- unique(totals)
    if (length(unique_totals) != 1L) {
      stop(sprintf(
        "Read totals don't match across different sections. Found: %s",
        paste(sapply(totals, \(x) paste(names(x), x, sep = "=")), collapse = ", ")
      ))
    }
  }

  results <- list(
    file = sample_file,
    reads_total = reads_total |> unname(),
    reads_trimmed_start = stats_start["trimmed"] |> unname(),
    reads_trimmed_end = stats_end["trimmed"] |> unname(),
    reads_split = stats_middle["trimmed"] |> unname(),
    adapter_matches = adapter_matches,
    adapter_sequences = trimming_sequences
  )

  results
}
