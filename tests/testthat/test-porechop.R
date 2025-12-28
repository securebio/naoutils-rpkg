test_that("read_porechop_log parses valid log file", {
  log_file <- test_path("fixtures", "porechop_sample.log")
  result <- read_porechop_log(log_file)

  # Check structure
  expect_type(result, "list")
  expect_named(
    result,
    c(
      "reads_total", "reads_trimmed_start",
      "reads_trimmed_end", "reads_split",
      "adapter_matches", "adapter_sequences"
    )
  )

  # Check types
  expect_type(result$reads_total, "integer")
  expect_type(result$reads_trimmed_start, "integer")
  expect_type(result$reads_trimmed_end, "integer")
  expect_type(result$reads_split, "integer")
  expect_s3_class(result$adapter_matches, "tbl_df")
  expect_s3_class(result$adapter_sequences, "tbl_df")
})

test_that("read_porechop_log validates input", {
  expect_error(
    read_porechop_log("nonexistent_file.log"),
    "File not found"
  )
})

test_that("read_porechop_log extracts correct statistics", {
  log_file <- test_path("fixtures", "porechop_sample.log")
  result <- read_porechop_log(log_file)

  # Check that all counts are non-negative
  expect_gte(result$reads_total, 0)
  expect_gte(result$reads_trimmed_start, 0)
  expect_gte(result$reads_trimmed_end, 0)
  expect_gte(result$reads_split, 0)

  # Check that trimmed counts don't exceed total
  expect_lte(result$reads_trimmed_start, result$reads_total)
  expect_lte(result$reads_trimmed_end, result$reads_total)
  expect_lte(result$reads_split, result$reads_total)
})

test_that("read_porechop_log adapter_matches has correct structure", {
  log_file <- test_path("fixtures", "porechop_sample.log")
  result <- read_porechop_log(log_file)

  # Check adapter_matches columns
  expect_true(
    all(c("adapter_set", "best_start_perc", "best_end_perc", "selected") %in%
      names(result$adapter_matches))
  )

  # Check column types
  expect_type(result$adapter_matches$adapter_set, "character")
  expect_type(result$adapter_matches$best_start_perc, "double")
  expect_type(result$adapter_matches$best_end_perc, "double")
  expect_type(result$adapter_matches$selected, "logical")

  # Check that there are multiple adapter sets
  expect_gt(nrow(result$adapter_matches), 0)

  # Check that percentages are in valid range (0-100)
  expect_true(all(result$adapter_matches$best_start_perc >= 0))
  expect_true(all(result$adapter_matches$best_start_perc <= 100))
  expect_true(all(result$adapter_matches$best_end_perc >= 0))
  expect_true(all(result$adapter_matches$best_end_perc <= 100))
})

test_that("read_porechop_log detects selected adapters", {
  log_file <- test_path("fixtures", "porechop_sample.log")
  result <- read_porechop_log(log_file)

  # Check that at least one adapter is selected
  expect_true(any(result$adapter_matches$selected))

  # Check that adapter_sequences is not empty
  expect_gt(nrow(result$adapter_sequences), 0)
})

test_that("read_porechop_log adapter_sequences has correct structure", {
  log_file <- test_path("fixtures", "porechop_sample.log")
  result <- read_porechop_log(log_file)

  # Check adapter_sequences columns
  expect_true(
    all(c("sequence_id", "sequence") %in% names(result$adapter_sequences))
  )

  # Check column types
  expect_type(result$adapter_sequences$sequence_id, "character")
  expect_type(result$adapter_sequences$sequence, "character")

  # Check that sequences are non-empty
  expect_true(all(nchar(result$adapter_sequences$sequence) > 0))
  expect_true(all(nchar(result$adapter_sequences$sequence_id) > 0))
})

test_that("read_porechop_log returns expected values for sample file", {
  log_file <- test_path("fixtures", "porechop_sample.log")
  result <- read_porechop_log(log_file)

  # These are the expected values from the MultiQC test file
  expect_equal(result$reads_total, 10000)
  expect_equal(result$reads_trimmed_start, 7100)
  expect_equal(result$reads_trimmed_end, 4849)
  expect_equal(result$reads_split, 7)

  # Check that specific adapters are selected
  selected_adapters <- result$adapter_matches |>
    dplyr::filter(selected) |>
    dplyr::pull(adapter_set)

  expect_true("SQK-NSK007" %in% selected_adapters)
  expect_true("SQK-MAP006 short" %in% selected_adapters)
})
