
#' Read one or more Kraken2 sample reports
#'
#' This function can read one or more Kraken(2) sample reports, as well as
#' output of `kraken2-inspect`. Kraken2 reports with and without minimizer data
#' are supported.
#'
#' @details
#' The column names in the resulting data frame are "clade_fragments_percent",
#' "clade_fragments", "node_fragments", "minimizers_total",
#' "minimizers_distinct", "rank_code", "taxid", "scientific_name"; the
#' minimizer columns are only included if `minimizers = TRUE`.
#'
#' This function should work on Kraken and Kraken2 reports as well as on the
#' output of `kraken2-inspect`, but we've not tested it on original Kraken
#' reports. Currently there is no check as to whether the minimizer data is in
#' the reports or not; if `minimizers` is set incorrectly, the function will
#' return without error but the reports will be parsed incorrectly.
#'
#' - [Kraken2 docs](https://github.com/DerrickWood/kraken2/wiki/Manual#sample-report-output-format)
#' - [Kraken docs](https://ccb.jhu.edu/software/kraken/MANUAL.html#output-format)
#'
#' @param files Paths to files
#' @param minimizers Boolean specifying whether the files include minimizer
#'   data
#' @param id A string or NULL. If a string, create column with that name for
#'   storing the files or the names of the `files` vector
#' @param num_threads Parameter passed onto `readr::read_tsv()`
#'
#' @export
read_kraken_reports <- function(files,
                                minimizers,
                                id = "file",
                                num_threads = 1) {
  if (minimizers) {
    cts <- "nnnnncic"
    cns <- kraken_col_names[1:8]
  } else {
    cns <- kraken_col_names[c(1:3, 6:8)]
    cts <- "nnncic"
  }
  if (is.null(names(files))) {
    names(files) <- files
  }
  if (is.null(id)) {
    id <- rlang::zap()
  }
  # Note: We call map on the list of functions, instead of passing the function
  # list to `read_tsv()`, because `read_tsv()` can fail when the multi-file
  # input is too large. In future, consider checking if furrr is installed, and
  # if so use future_map.
  files |>
    map(\(x) readr::read_tsv(x,
      col_names = cns,
      col_types = cts,
      comment = "#", # Allows reading Kraken inspect.txt files
      trim_ws = FALSE,
      progress = FALSE,
      num_threads = num_threads
    )) |>
    list_rbind(names_to = id) |>
    # Store taxonomic level info, then clean up scientific name
    mutate(
      rank_level = scientific_name |> str_extract("^ *") |> str_length() |>
        (\(x) x %/% 2L)(),
      across(scientific_name, \(x) str_replace(x, "^ *", ""))
    )
}

kraken_col_names <- c(
  "clade_fragments_percent",
  "clade_fragments",
  "node_fragments",
  "minimizers_total", # Used when `minimizers = TRUE`
  "minimizers_distinct", # Used when `minimizers = TRUE`
  "rank_code",
  "taxid",
  "scientific_name",
  "rank_level", # created by `read_kraken_reports()`
  "taxonomy" # created by `kraken_add_taxonomy()`
)

#' Complete a set of Kraken reports by filling in zero counts (Not tested!)
#'
#' This function uses `tidyr::complete()` to add rows with zeros for taxa that
#' missing in some samples. Note: The resulting data frame will no longer be
#' compatible with `kraken_add_taxonomy()`.
#'
#' @param x A data frame of Kraken sample reports
#' @param explicit Boolean passed to `tidyr::complete()`
#'
#' @export
kraken_complete <- function(x, explicit = TRUE) {
  stopifnot(inherits(x, "data.frame"))
  # Fill columns are quantities for a given (sample, taxon) pair
  fill_cols <- head(kraken_col_names, 5) |>
    intersect(names(x))
  fill_list <- fill_cols |>
    rep_named(0) |>
    as.list()
  # Use nesting() to group columns denoting taxon properties
  nesting_cols <- tail(kraken_col_names, 5) |>
    intersect(names(x))
  # Assume other columns are unique to the specific microbiome sample
  sample_cols <- names(x) |> setdiff(fill_cols) |> setdiff(nesting_cols)
  x |>
    tidyr::complete(
      !!!syms(sample_cols),
      tidyr::nesting(!!!syms(nesting_cols)),
      fill = fill_list,
      explicit = explicit
    )
}

#' Remove nonstandard ranks from Kraken reports
#'
#' @param x A data frame of Kraken sample reports
#'
#' @export
kraken_remove_nonstandard_ranks <- function(x) {
  stopifnot(inherits(x, "data.frame"))
  stopifnot("rank_code" %in% names(x))
  x |>
    dplyr::filter(stringr::str_detect(rank_code, "[0-9]", negate = TRUE))
}

#' Add full taxonomy path to a Kraken report data frame
#'
#' This function requires that rows are in the order of the original Kraken
#' report, so cannot be used after `kraken_complete()`.
#'
#' You can apply this function to a data frame of Kraken reports for multiple
#' samples; however, doing so is typically inefficient due to needing to
#' recompute paths for taxa that appear in multiple samples. Instead, consider
#' generating taxonomy paths from the original database and joining these to
#' the kraken report data frame.
#'
#' @param x A data frame containing the fields "rank_code", "scientific_name",
#'   and "rank_level" with rows in the original order of the Kraken sample
#'   report(s), as returned by `read_kraken_reports()`.
#'
#' @export
kraken_add_taxonomy <- function(x) {
  stopifnot(inherits(x, "data.frame"))
  stopifnot(all(c("rank_code", "scientific_name", "rank_level") %in% names(x)))

  lvls <- x$rank_level

  get_clade_idx <- function(i) {
    # TODO: consider if get a speedup by testing if have any children (which can do by looking at the next row)
    next_non_child <- match(TRUE, seq_along(lvls) > i & lvls <= lvls[i], nomatch = length(lvls) + 1)
    seq.int(i, next_non_child - 1)
  }

  # Iterate through the nodes, adding the node's string to the tax path for it and its children
  x <- x |> add_column(taxonomy = "")
  for (i in seq.int(nrow(x))) {
    node_string <- str_c(x[[i, "rank_code"]], "__", x[[i, "scientific_name"]], ";")
    clade_idx <- get_clade_idx(i)
    x[clade_idx, "taxonomy"] <- str_c(x[clade_idx, ][["taxonomy"]], node_string)
  }

  x
}
