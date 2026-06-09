#!/usr/bin/env Rscript

###############################################################################
# Transparent, reproducible Q-methodology workflow around the R package qmethod
# by trippm@tripplab.com on 090626
#
# Purpose:
#   This script reads a project-specific Q-sort CSV file, validates the forced
#   distribution, runs qmethod solutions with 1, 2, 3, and 4 factors, exports
#   loadings, defining sorts, factor arrays, distinguishing/consensus statements,
#   and generates a strict LaTeX solution-comparison report.
#
# Expected input CSV format:
#
#   statement_id, statement_text, category, P1, P2, P3, P4, P5, P6, P7, P8
#
# Expected Q-sort design:
#
#   -6  
#   -5  -5  
#   -4  -4  -4  
#   -3  -3  -3  -3
#   -2  -2  -2  -2  -2  
#   -1  -1  -1  -1  -1  -1
#    0   0   0   0   0   0   0
#    1   1   1   1   1   1
#    2   2   2   2   2
#    3   3   3   3
#    4   4   4
#    5   5
#    6
#
# Important patch relative to the earlier version:
#   The qmethod package can return a non-table qdc component for a 1-factor
#   solution because distinguishing/consensus statements are not defined when
#   there is only one factor. Calling print(result) or summary(result) may then
#   crash through print.QmethodRes().
#
#   Therefore, this script NEVER directly calls print(result) or summary(result)
#   on the qmethod object. It extracts and writes components manually instead.
#
###############################################################################

###############################################################################
# 0. USER SETTINGS
###############################################################################

# Path to the input CSV file.
input_file <- "data/qsort_data.csv"

# Main output directory.
output_dir <- "outputs"

# Expected number of statements.
expected_n_statements <- 49L

# Expected participant columns.
expected_participants <- paste0("P", 1:8)

# Forced Q-sort distribution.
expected_distribution <- c(
  rep(-6, 1),
  rep(-5, 2),
  rep(-4, 3),
  rep(-3, 4),
  rep(-2, 5),
  rep(-1, 6),
  rep( 0, 7),
  rep( 1, 6),
  rep( 2, 5),
  rep( 3, 4),
  rep( 4, 3),
  rep( 5, 2),
  rep( 6, 1)
)

# Factor solutions to run.
factor_solutions <- 1:4

# qmethod settings.
#
# PCA is used as the default extraction method.
# Varimax is used for all solutions, including 1F, because qmethod accepts
# varimax as a standard rotation argument. For 1F, rotation has no substantive
# interpretive role, but using varimax avoids nonstandard-rotation warnings
# caused by passing rotation = "none".
extraction_method <- "PCA"
rotation_method <- "varimax"
correlation_method <- "pearson"

# Factor-loading significance threshold.
# For 49 statements:
#   SE = 1 / sqrt(49) = 0.143
#   p < 0.01 threshold = 2.58 / sqrt(49) = 0.369
loading_significance_threshold <- 2.58 / sqrt(expected_n_statements)

###############################################################################
# 1. PACKAGE CHECK
###############################################################################

required_packages <- c("qmethod")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", "),
      "\nInstall them with:\n",
      "install.packages(c(",
      paste(sprintf("\"%s\"", missing_packages), collapse = ", "),
      "))"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages(library(qmethod))

###############################################################################
# 2. OUTPUT DIRECTORIES
###############################################################################

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

subdirs <- c(
  "validation",
  "correlations",
  "factor_solutions",
  "factor_arrays",
  "distinguishing_consensus",
  "reports",
  "rds"
)

for (d in subdirs) {
  dir.create(file.path(output_dir, d), recursive = TRUE, showWarnings = FALSE)
}

###############################################################################
# 3. GENERAL HELPER FUNCTIONS
###############################################################################

# Write a CSV using UTF-8 encoding.
write_csv_utf8 <- function(x, path, row.names = FALSE) {
  utils::write.csv(
    x,
    file = path,
    row.names = row.names,
    fileEncoding = "UTF-8",
    na = ""
  )
}

# Capture object printing safely.
# This prevents a printing method from crashing the whole analysis.
safe_capture_print <- function(x) {
  out <- tryCatch(
    {
      capture.output(print(x))
    },
    error = function(e) {
      c(
        "Printing failed.",
        paste0("Error: ", conditionMessage(e))
      )
    }
  )
  out
}

# Extract a component from a qmethod result.
# qmethod objects are lists, but package versions may differ slightly.
# This function first tries by name, then by list index as fallback.
get_component <- function(result, name, index = NULL) {
  if (!is.null(result[[name]])) {
    return(result[[name]])
  }

  if (!is.null(index) && length(result) >= index) {
    return(result[[index]])
  }

  return(NULL)
}

# Convert object to data frame if possible.
safe_as_data_frame <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }

  out <- tryCatch(
    {
      as.data.frame(x, stringsAsFactors = FALSE)
    },
    error = function(e) {
      NULL
    }
  )

  return(out)
}

# Select the first nfactors numeric columns from a qmethod component.
# This is useful for loadings, z-scores, and rounded factor scores.
select_numeric_factor_columns <- function(x, nfactors, prefix = "F") {
  x <- safe_as_data_frame(x)

  if (is.null(x)) {
    return(NULL)
  }

  numeric_cols <- vapply(x, is.numeric, logical(1))

  if (sum(numeric_cols) < nfactors) {
    return(NULL)
  }

  out <- x[, which(numeric_cols)[seq_len(nfactors)], drop = FALSE]
  names(out) <- paste0(prefix, seq_len(nfactors))

  return(out)
}

# Escape text for LaTeX.
latex_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""

  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x <- gsub("~", "\\\\textasciitilde{}", x)

  return(x)
}

# Create a simple LaTeX table.
# This intentionally avoids HTML, Word, or PDF generation.
latex_table <- function(df, caption = NULL, label = NULL, align = NULL) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)

  if (is.null(align)) {
    align <- paste0(rep("l", ncol(df)), collapse = "")
  }

  lines <- character(0)

  lines <- c(lines, "\\begin{table}[htbp]")
  lines <- c(lines, "\\centering")

  if (!is.null(caption)) {
    lines <- c(lines, paste0("\\caption{", latex_escape(caption), "}"))
  }

  if (!is.null(label)) {
    lines <- c(lines, paste0("\\label{", latex_escape(label), "}"))
  }

  lines <- c(lines, paste0("\\begin{tabular}{", align, "}"))
  lines <- c(lines, "\\hline")

  header <- paste(latex_escape(names(df)), collapse = " & ")
  lines <- c(lines, paste0(header, " \\\\"))
  lines <- c(lines, "\\hline")

  for (i in seq_len(nrow(df))) {
    row <- paste(latex_escape(unlist(df[i, ], use.names = FALSE)), collapse = " & ")
    lines <- c(lines, paste0(row, " \\\\"))
  }

  lines <- c(lines, "\\hline")
  lines <- c(lines, "\\end{tabular}")
  lines <- c(lines, "\\end{table}")

  return(paste(lines, collapse = "\n"))
}

###############################################################################
# 4. Q-METHODOLOGY HELPER FUNCTIONS
###############################################################################

# Run qmethod robustly.
#
# The current qmethod interface is typically:
#   qmethod(dataset, nfactors, extraction = "PCA", rotation = "varimax",
#           forced = TRUE, distribution = NA, cor.method = "pearson",
#           silent = FALSE, ...)
#
# Older versions may expose nstat/nqsorts or differ slightly, so this function
# checks formal arguments before passing optional parameters.
run_qmethod_solution <- function(q_matrix, nfactors) {
  fmls <- names(formals(qmethod::qmethod))

  args <- list(
    dataset = as.data.frame(q_matrix),
    nfactors = nfactors
  )

  if ("extraction" %in% fmls) {
    args$extraction <- extraction_method
  }

  if ("rotation" %in% fmls) {
    args$rotation <- rotation_method
  }

  if ("forced" %in% fmls) {
    args$forced <- TRUE
  }

  if ("cor.method" %in% fmls) {
    args$cor.method <- correlation_method
  }

  if ("silent" %in% fmls) {
    args$silent <- TRUE
  }

  if ("nstat" %in% fmls) {
    args$nstat <- nrow(q_matrix)
  }

  if ("nqsorts" %in% fmls) {
    args$nqsorts <- ncol(q_matrix)
  }

  # Do not pass distribution unless a non-forced design is needed.
  # For forced = TRUE, qmethod calculates the distribution automatically
  # from the observed Q-sort data.
  result <- do.call(qmethod::qmethod, args)

  return(result)
}

# Convert qmethod's flagged matrix to a logical matrix with known row/column names.
as_logical_flag_matrix <- function(x, participants, nfactors) {
  row_names <- participants
  col_names <- paste0("F", seq_len(nfactors))

  if (is.null(x)) {
    out <- matrix(
      FALSE,
      nrow = length(row_names),
      ncol = nfactors,
      dimnames = list(row_names, col_names)
    )
    return(out)
  }

  x_df <- safe_as_data_frame(x)

  if (is.null(x_df) || ncol(x_df) < nfactors) {
    out <- matrix(
      FALSE,
      nrow = length(row_names),
      ncol = nfactors,
      dimnames = list(row_names, col_names)
    )
    return(out)
  }

  x_df <- x_df[, seq_len(nfactors), drop = FALSE]
  out <- as.matrix(x_df)

  # Force logical interpretation.
  storage.mode(out) <- "logical"

  rownames(out) <- row_names
  colnames(out) <- col_names

  return(out)
}

# Fallback flagging in case qmethod does not return usable flags.
#
# This follows a conservative simple rule:
#   - A Q-sort is defining only if it loads significantly on exactly one factor.
#   - If it loads significantly on more than one factor, it is confounded.
#   - If it loads significantly on no factor, it is non-defining.
#
# This fallback is not meant to replace qmethod's own qflag() when available.
fallback_flags_from_loadings <- function(loadings_numeric, threshold) {
  L <- as.matrix(loadings_numeric)
  abs_L <- abs(L)

  significant <- abs_L >= threshold

  out <- matrix(
    FALSE,
    nrow = nrow(significant),
    ncol = ncol(significant),
    dimnames = dimnames(significant)
  )

  for (i in seq_len(nrow(significant))) {
    sig_idx <- which(significant[i, ])

    if (length(sig_idx) >= 1) {
      out[i, sig_idx] <- TRUE
    }
  }

  return(out)
}

# Classify Q-sorts as defining, confounded, or non-defining.
classify_defining_sorts <- function(loadings_numeric, flagged_matrix) {
  participants <- rownames(loadings_numeric)
  factor_names <- colnames(loadings_numeric)

  flag_count <- rowSums(flagged_matrix, na.rm = TRUE)

  status <- ifelse(
    flag_count == 0,
    "non_defining",
    ifelse(flag_count == 1, "defining", "confounded")
  )

  defining_factor <- rep(NA_character_, length(participants))
  defining_loading <- rep(NA_real_, length(participants))
  loading_sign <- rep(NA_character_, length(participants))

  for (i in seq_along(participants)) {
    idx <- which(flagged_matrix[i, ])

    if (length(idx) == 1) {
      defining_factor[i] <- factor_names[idx]
      defining_loading[i] <- as.matrix(loadings_numeric)[i, idx]
      loading_sign[i] <- ifelse(defining_loading[i] >= 0, "positive", "negative")
    }

    if (length(idx) > 1) {
      defining_factor[i] <- paste(factor_names[idx], collapse = ";")
      defining_loading[i] <- NA_real_
      loading_sign[i] <- "confounded"
    }
  }

  out <- data.frame(
    participant = participants,
    status = status,
    defining_factor = defining_factor,
    defining_loading = defining_loading,
    loading_sign = loading_sign,
    stringsAsFactors = FALSE
  )

  return(out)
}

# Assign forced-distribution scores from z-scores.
#
# This is used only as a fallback if qmethod does not provide rounded
# factor scores in a clean numeric table.
assign_forced_scores_from_z <- function(z_vector, distribution_vector) {
  scores_sorted <- sort(distribution_vector, decreasing = TRUE)
  idx <- order(z_vector, decreasing = TRUE, na.last = TRUE)

  out <- rep(NA_real_, length(z_vector))
  out[idx] <- scores_sorted

  return(out)
}

# Build a clean factor array table from qmethod outputs.
#
# The output includes:
#   statement_id
#   statement_text
#   category
#   F1_z, F2_z, ...
#   F1_score, F2_score, ...
standardize_factor_array <- function(result, metadata, nfactors, distribution_vector) {
  # qmethod usually stores statement z-scores in result$zsc.
  zsc_raw <- get_component(result, "zsc", index = 5)

  z_df <- select_numeric_factor_columns(zsc_raw, nfactors, prefix = "F")

  if (is.null(z_df)) {
    stop(
      paste0("Could not extract statement z-scores for ", nfactors, "-factor solution."),
      call. = FALSE
    )
  }

  names(z_df) <- paste0("F", seq_len(nfactors), "_z")

  # qmethod usually stores rounded factor scores in result$zsc_n.
  zsc_n_raw <- get_component(result, "zsc_n", index = 6)
  score_df <- select_numeric_factor_columns(zsc_n_raw, nfactors, prefix = "F")

  if (!is.null(score_df)) {
    names(score_df) <- paste0("F", seq_len(nfactors), "_score")
  }

  # Fallback: assign forced scores by ranking z-scores.
  if (is.null(score_df)) {
    score_df <- as.data.frame(
      matrix(
        NA_real_,
        nrow = nrow(z_df),
        ncol = nfactors
      )
    )

    names(score_df) <- paste0("F", seq_len(nfactors), "_score")

    for (j in seq_len(nfactors)) {
      score_df[[j]] <- assign_forced_scores_from_z(
        z_vector = z_df[[j]],
        distribution_vector = distribution_vector
      )
    }
  }

  out <- cbind(
    metadata,
    z_df,
    score_df,
    stringsAsFactors = FALSE
  )

  return(out)
}

# Extract distinguishing / consensus table safely.
#
# For nfactors == 1, this is not applicable.
# For nfactors >= 2, qmethod normally provides result$qdc.
# The function handles unexpected structures safely.
extract_qdc_table <- function(result, nfactors) {
  if (nfactors < 2) {
    return(NULL)
  }

  qdc_raw <- get_component(result, "qdc", index = 8)

  qdc_df <- safe_as_data_frame(qdc_raw)

  if (is.null(qdc_df)) {
    return(NULL)
  }

  return(qdc_df)
}

# Count distinguishing and consensus statements from a qdc table.
#
# qmethod versions may differ in column names, so this function scans all
# character columns for likely labels.
count_qdc_labels <- function(qdc_df) {
  if (is.null(qdc_df) || nrow(qdc_df) == 0) {
    return(list(distinguishing = NA_integer_, consensus = NA_integer_))
  }

  char_cols <- vapply(qdc_df, is.character, logical(1))

  if (!any(char_cols)) {
    return(list(distinguishing = NA_integer_, consensus = NA_integer_))
  }

  txt <- unlist(qdc_df[, char_cols, drop = FALSE], use.names = FALSE)
  txt <- txt[!is.na(txt)]

  n_distinguishing <- sum(grepl("distinguish", txt, ignore.case = TRUE))
  n_consensus <- sum(grepl("consensus", txt, ignore.case = TRUE))

  return(list(
    distinguishing = n_distinguishing,
    consensus = n_consensus
  ))
}

# Write a manual, safe, text-only summary of one qmethod result.
#
# This replaces print(result) and summary(result), which can crash for 1F
# because qdc is not defined as a table for a one-factor solution.
write_safe_qmethod_summary <- function(
  result,
  nfactors,
  solution_label,
  loadings_export,
  defining_sorts,
  factor_array,
  path
) {
  lines <- character(0)

  lines <- c(lines, paste0("qmethod solution: ", solution_label))
  lines <- c(lines, paste0("Number of factors: ", nfactors))
  lines <- c(lines, paste0("Extraction method: ", extraction_method))
  lines <- c(lines, paste0("Rotation argument passed to qmethod: ", rotation_method))

  if (nfactors == 1) {
    lines <- c(
      lines,
      "",
      "Note: This is a 1-factor solution.",
      "Distinguishing and consensus statements are not applicable because there is only one factor.",
      "Direct print(result) and summary(result) are intentionally skipped to avoid qmethod print-method errors when qdc is not table-like."
    )
  } else {
    lines <- c(
      lines,
      "",
      "Note: This is a multi-factor solution.",
      "Distinguishing and consensus statements are applicable if qmethod returned a qdc table."
    )
  }

  lines <- c(lines, "", "Available qmethod result components:")
  lines <- c(lines, paste(names(result), collapse = ", "))

  lines <- c(lines, "", "Loadings:")
  lines <- c(lines, capture.output(print(loadings_export)))

  lines <- c(lines, "", "Defining-sort classification:")
  lines <- c(lines, capture.output(print(defining_sorts)))

  lines <- c(lines, "", "Factor array preview:")
  preview_cols <- c("statement_id", "statement_text", "category")
  z_cols <- grep("_z$", names(factor_array), value = TRUE)
  score_cols <- grep("_score$", names(factor_array), value = TRUE)
  preview_cols <- c(preview_cols, z_cols, score_cols)
  preview <- head(factor_array[, preview_cols, drop = FALSE], 10)
  lines <- c(lines, capture.output(print(preview)))

  f_char <- get_component(result, "f_char", index = 7)
  if (!is.null(f_char)) {
    lines <- c(lines, "", "Factor characteristics:")
    lines <- c(lines, safe_capture_print(f_char))
  }

  writeLines(lines, con = path, useBytes = TRUE)
}

###############################################################################
# 5. READ INPUT CSV
###############################################################################

if (!file.exists(input_file)) {
  stop(
    paste0("Input file not found: ", input_file),
    call. = FALSE
  )
}

raw_data <- utils::read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

###############################################################################
# 6. VALIDATE STRUCTURE
###############################################################################

validation_messages <- character(0)
validation_pass <- TRUE

required_metadata_cols <- c("statement_id", "statement_text", "category")

# Check metadata columns.
missing_metadata <- setdiff(required_metadata_cols, names(raw_data))

if (length(missing_metadata) > 0) {
  validation_pass <- FALSE
  validation_messages <- c(
    validation_messages,
    paste0(
      "Missing required metadata column(s): ",
      paste(missing_metadata, collapse = ", ")
    )
  )
}

# Check participant columns.
missing_participants <- setdiff(expected_participants, names(raw_data))

if (length(missing_participants) > 0) {
  validation_pass <- FALSE
  validation_messages <- c(
    validation_messages,
    paste0(
      "Missing expected participant column(s): ",
      paste(missing_participants, collapse = ", ")
    )
  )
}

# Check number of statements.
if (nrow(raw_data) != expected_n_statements) {
  validation_pass <- FALSE
  validation_messages <- c(
    validation_messages,
    paste0(
      "Invalid number of statements. Expected ",
      expected_n_statements,
      ", observed ",
      nrow(raw_data),
      "."
    )
  )
}

# Stop early if critical columns are missing.
if (!all(c(required_metadata_cols, expected_participants) %in% names(raw_data))) {
  validation_report <- data.frame(
    check = "CSV structure",
    status = "FAIL",
    message = paste(validation_messages, collapse = " | "),
    stringsAsFactors = FALSE
  )

  write_csv_utf8(
    validation_report,
    file.path(output_dir, "validation", "input_validation_report.csv")
  )

  stop(
    paste0(
      "Input validation failed before Q-sort checks.\n",
      paste(validation_messages, collapse = "\n")
    ),
    call. = FALSE
  )
}

###############################################################################
# 7. SEPARATE METADATA AND RAW Q-SORT DATA
###############################################################################

metadata <- raw_data[, required_metadata_cols, drop = FALSE]
q_data_raw <- raw_data[, expected_participants, drop = FALSE]

# Validate statement IDs.
if (any(is.na(metadata$statement_id) | metadata$statement_id == "")) {
  validation_pass <- FALSE
  validation_messages <- c(validation_messages, "Missing statement_id value(s).")
}

if (any(duplicated(metadata$statement_id))) {
  validation_pass <- FALSE
  validation_messages <- c(validation_messages, "Duplicated statement_id value(s).")
}

# Validate statement text.
if (any(is.na(metadata$statement_text) | metadata$statement_text == "")) {
  validation_pass <- FALSE
  validation_messages <- c(validation_messages, "Missing statement_text value(s).")
}

# Validate category.
if (any(is.na(metadata$category) | metadata$category == "")) {
  validation_pass <- FALSE
  validation_messages <- c(validation_messages, "Missing category value(s).")
}

###############################################################################
# 8. CONVERT Q-SORT DATA TO NUMERIC MATRIX
###############################################################################

q_data_numeric <- as.data.frame(
  lapply(q_data_raw, function(x) suppressWarnings(as.numeric(as.character(x)))),
  stringsAsFactors = FALSE
)

# Identify cells that failed numeric conversion.
non_numeric_cells <- is.na(q_data_numeric) & !is.na(q_data_raw)

if (any(non_numeric_cells)) {
  validation_pass <- FALSE
  validation_messages <- c(
    validation_messages,
    "Some Q-sort cells could not be converted to numeric values."
  )
}

# Identify missing Q-sort scores.
if (any(is.na(q_data_numeric))) {
  validation_pass <- FALSE
  validation_messages <- c(
    validation_messages,
    "Missing Q-sort score(s) detected."
  )
}

# Allowed scores.
allowed_scores <- sort(unique(expected_distribution))

invalid_score_mask <- !(as.matrix(q_data_numeric) %in% allowed_scores)

if (any(invalid_score_mask, na.rm = TRUE)) {
  validation_pass <- FALSE
  validation_messages <- c(
    validation_messages,
    paste0(
      "Invalid Q-sort score(s) detected. Allowed values are: ",
      paste(allowed_scores, collapse = ", "),
      "."
    )
  )
}

# Check integer scores.
integer_score_mask <- abs(
  as.matrix(q_data_numeric) - round(as.matrix(q_data_numeric))
) > .Machine$double.eps^0.5

if (any(integer_score_mask, na.rm = TRUE)) {
  validation_pass <- FALSE
  validation_messages <- c(
    validation_messages,
    "Non-integer Q-sort score(s) detected."
  )
}

# Create numeric Q matrix.
q_matrix <- as.matrix(q_data_numeric)
storage.mode(q_matrix) <- "numeric"

rownames(q_matrix) <- as.character(metadata$statement_id)
colnames(q_matrix) <- expected_participants

###############################################################################
# 9. VALIDATE FORCED DISTRIBUTION
###############################################################################

expected_counts <- as.integer(
  table(factor(expected_distribution, levels = allowed_scores))
)

distribution_check <- data.frame(
  score = allowed_scores,
  expected = expected_counts,
  stringsAsFactors = FALSE
)

for (participant in expected_participants) {
  observed_counts <- as.integer(
    table(factor(q_matrix[, participant], levels = allowed_scores))
  )

  distribution_check[[participant]] <- observed_counts
}

# Status by participant.
participant_distribution_status <- data.frame(
  participant = expected_participants,
  status = NA_character_,
  stringsAsFactors = FALSE
)

for (participant in expected_participants) {
  observed_counts <- as.integer(
    table(factor(q_matrix[, participant], levels = allowed_scores))
  )

  participant_distribution_status$status[
    participant_distribution_status$participant == participant
  ] <- ifelse(
    identical(observed_counts, expected_counts),
    "PASS",
    "FAIL"
  )
}

if (any(participant_distribution_status$status != "PASS")) {
  validation_pass <- FALSE
  validation_messages <- c(
    validation_messages,
    paste0(
      "Forced distribution failed for participant(s): ",
      paste(
        participant_distribution_status$participant[
          participant_distribution_status$status != "PASS"
        ],
        collapse = ", "
      )
    )
  )
}

# Status by score across participants.
distribution_check$status_by_score <- apply(
  distribution_check[, expected_participants, drop = FALSE],
  1,
  function(x) {
    if (all(x == distribution_check$expected[seq_along(x)])) {
      "PASS"
    } else {
      # This row-level check is less important than participant-level status.
      # It is kept only as a diagnostic.
      "CHECK_PARTICIPANT_COLUMNS"
    }
  }
)

###############################################################################
# 10. WRITE VALIDATION OUTPUTS
###############################################################################

write_csv_utf8(
  distribution_check,
  file.path(output_dir, "validation", "qsort_distribution_check.csv")
)

write_csv_utf8(
  participant_distribution_status,
  file.path(output_dir, "validation", "participant_distribution_status.csv")
)

validation_report <- data.frame(
  check = c(
    "required_metadata_columns",
    "participant_columns",
    "number_of_statements",
    "statement_ids",
    "statement_text",
    "category",
    "qsort_numeric",
    "qsort_missing_values",
    "qsort_allowed_scores",
    "qsort_integer_scores",
    "forced_distribution"
  ),
  status = c(
    ifelse(length(missing_metadata) == 0, "PASS", "FAIL"),
    ifelse(length(missing_participants) == 0, "PASS", "FAIL"),
    ifelse(nrow(raw_data) == expected_n_statements, "PASS", "FAIL"),
    ifelse(!any(is.na(metadata$statement_id) | metadata$statement_id == "" | duplicated(metadata$statement_id)), "PASS", "FAIL"),
    ifelse(!any(is.na(metadata$statement_text) | metadata$statement_text == ""), "PASS", "FAIL"),
    ifelse(!any(is.na(metadata$category) | metadata$category == ""), "PASS", "FAIL"),
    ifelse(!any(non_numeric_cells), "PASS", "FAIL"),
    ifelse(!any(is.na(q_data_numeric)), "PASS", "FAIL"),
    ifelse(!any(invalid_score_mask, na.rm = TRUE), "PASS", "FAIL"),
    ifelse(!any(integer_score_mask, na.rm = TRUE), "PASS", "FAIL"),
    ifelse(all(participant_distribution_status$status == "PASS"), "PASS", "FAIL")
  ),
  stringsAsFactors = FALSE
)

validation_report$message <- ""

if (length(validation_messages) > 0) {
  validation_report$message[validation_report$status == "FAIL"] <- paste(
    validation_messages,
    collapse = " | "
  )
}

write_csv_utf8(
  validation_report,
  file.path(output_dir, "validation", "input_validation_report.csv")
)

writeLines(
  c(
    "Q-methodology input validation report",
    paste0("Input file: ", input_file),
    paste0("Validation status: ", ifelse(validation_pass, "PASS", "FAIL")),
    "",
    if (length(validation_messages) == 0) {
      "No validation problems detected."
    } else {
      validation_messages
    }
  ),
  con = file.path(output_dir, "validation", "input_validation_report.txt"),
  useBytes = TRUE
)

# Stop before qmethod if validation failed.
if (!validation_pass) {
  stop(
    paste0(
      "Input validation failed. See: ",
      file.path(output_dir, "validation", "input_validation_report.txt")
    ),
    call. = FALSE
  )
}

###############################################################################
# 11. EXPORT PREPARED METADATA AND Q MATRIX
###############################################################################

write_csv_utf8(
  metadata,
  file.path(output_dir, "validation", "statement_metadata.csv")
)

write_csv_utf8(
  data.frame(statement_id = rownames(q_matrix), q_matrix, check.names = FALSE),
  file.path(output_dir, "validation", "q_matrix_numeric.csv")
)

###############################################################################
# 12. PARTICIPANT CORRELATION MATRIX
###############################################################################

# In Q methodology, participants/Q-sorts are correlated with each other.
# Since q_matrix columns are participants and rows are statements, cor(q_matrix)
# gives the participant-by-participant correlation matrix.
participant_correlation_matrix <- cor(
  q_matrix,
  method = correlation_method,
  use = "pairwise.complete.obs"
)

write_csv_utf8(
  data.frame(
    participant = rownames(participant_correlation_matrix),
    participant_correlation_matrix,
    check.names = FALSE
  ),
  file.path(output_dir, "correlations", "participant_correlation_matrix.csv")
)

###############################################################################
# 13. RUN QMETHOD SOLUTIONS
###############################################################################

solution_summaries <- list()
all_results <- list()

for (nf in factor_solutions) {
  solution_label <- paste0(nf, "F")

  message("Running qmethod solution: ", solution_label)

  ###########################################################################
  # 13.1 Run qmethod safely
  ###########################################################################

  result <- tryCatch(
    {
      run_qmethod_solution(q_matrix, nfactors = nf)
    },
    error = function(e) {
      stop(
        paste0(
          "qmethod failed for solution ",
          solution_label,
          ". Error: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  all_results[[solution_label]] <- result

  saveRDS(
    result,
    file = file.path(output_dir, "rds", paste0("qmethod_result_", solution_label, ".rds"))
  )

  ###########################################################################
  # 13.2 Extract loadings
  ###########################################################################

  # qmethod usually stores loadings in result$loa.
  loa_raw <- get_component(result, "loa", index = 3)
  loadings_numeric <- select_numeric_factor_columns(
    loa_raw,
    nfactors = nf,
    prefix = "F"
  )

  if (is.null(loadings_numeric)) {
    stop(
      paste0("Could not extract loadings for solution ", solution_label, "."),
      call. = FALSE
    )
  }

  rownames(loadings_numeric) <- expected_participants
  names(loadings_numeric) <- paste0("F", seq_len(nf))

  ###########################################################################
  # 13.3 Extract or reconstruct qmethod flags
  ###########################################################################

  flagged_raw <- get_component(result, "flagged", index = 4)

  flagged_matrix <- as_logical_flag_matrix(
    x = flagged_raw,
    participants = expected_participants,
    nfactors = nf
  )

  # If no usable flagging was returned, use the fallback threshold rule.
  if (all(flagged_matrix == FALSE)) {
    flagged_matrix <- fallback_flags_from_loadings(
      loadings_numeric = loadings_numeric,
      threshold = loading_significance_threshold
    )
  }

  flagged_df <- as.data.frame(flagged_matrix, stringsAsFactors = FALSE)
  names(flagged_df) <- paste0("flag_", names(flagged_df))

  loadings_export <- data.frame(
    participant = expected_participants,
    loadings_numeric,
    flagged_df,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  write_csv_utf8(
    loadings_export,
    file.path(output_dir, "factor_solutions", paste0("loadings_", solution_label, ".csv"))
  )

  ###########################################################################
  # 13.4 Classify defining, confounded, and non-defining sorts
  ###########################################################################

  defining_sorts <- classify_defining_sorts(
    loadings_numeric = loadings_numeric,
    flagged_matrix = flagged_matrix
  )

  write_csv_utf8(
    defining_sorts,
    file.path(output_dir, "factor_solutions", paste0("defining_sorts_", solution_label, ".csv"))
  )

  ###########################################################################
  # 13.5 Export factor characteristics if available
  ###########################################################################

  f_char <- get_component(result, "f_char", index = 7)

  if (!is.null(f_char)) {
    writeLines(
      safe_capture_print(f_char),
      con = file.path(output_dir, "factor_solutions", paste0("factor_characteristics_", solution_label, ".txt")),
      useBytes = TRUE
    )
  }

  ###########################################################################
  # 13.6 Build and export factor arrays
  ###########################################################################

  factor_array <- standardize_factor_array(
    result = result,
    metadata = metadata,
    nfactors = nf,
    distribution_vector = expected_distribution
  )

  write_csv_utf8(
    factor_array,
    file.path(output_dir, "factor_arrays", paste0("factor_arrays_", solution_label, ".csv"))
  )

  # Export one ranked table per factor.
  for (j in seq_len(nf)) {
    f_name <- paste0("F", j)
    score_col <- paste0(f_name, "_score")
    z_col <- paste0(f_name, "_z")

    ranked_array <- factor_array[
      order(factor_array[[score_col]], factor_array[[z_col]], decreasing = TRUE),
      c("statement_id", "statement_text", "category", z_col, score_col),
      drop = FALSE
    ]

    write_csv_utf8(
      ranked_array,
      file.path(
        output_dir,
        "factor_arrays",
        paste0("factor_array_", solution_label, "_", f_name, "_ranked.csv")
      )
    )
  }

  ###########################################################################
  # 13.7 Export distinguishing and consensus statements where applicable
  ###########################################################################

  qdc_export <- NULL
  n_distinguishing <- NA_integer_
  n_consensus <- NA_integer_

  if (nf >= 2) {
    qdc_table <- extract_qdc_table(
      result = result,
      nfactors = nf
    )

    if (!is.null(qdc_table)) {
      # If qdc has one row per statement, attach metadata.
      if (nrow(qdc_table) == nrow(metadata)) {
        qdc_export <- cbind(metadata, qdc_table, stringsAsFactors = FALSE)
      } else {
        qdc_export <- qdc_table
      }

      write_csv_utf8(
        qdc_export,
        file.path(
          output_dir,
          "distinguishing_consensus",
          paste0("distinguishing_consensus_", solution_label, ".csv")
        )
      )

      qdc_counts <- count_qdc_labels(qdc_export)
      n_distinguishing <- qdc_counts$distinguishing
      n_consensus <- qdc_counts$consensus
    } else {
      writeLines(
        "Distinguishing / consensus table could not be extracted from qmethod result.",
        con = file.path(
          output_dir,
          "distinguishing_consensus",
          paste0("distinguishing_consensus_", solution_label, "_not_available.txt")
        ),
        useBytes = TRUE
      )
    }
  } else {
    writeLines(
      c(
        "Not applicable for a 1-factor solution.",
        "There are no between-factor comparisons when only one factor is extracted."
      ),
      con = file.path(
        output_dir,
        "distinguishing_consensus",
        "distinguishing_consensus_1F_not_applicable.txt"
      ),
      useBytes = TRUE
    )
  }

  ###########################################################################
  # 13.8 Write safe text summary for this qmethod solution
  ###########################################################################

  write_safe_qmethod_summary(
    result = result,
    nfactors = nf,
    solution_label = solution_label,
    loadings_export = loadings_export,
    defining_sorts = defining_sorts,
    factor_array = factor_array,
    path = file.path(
      output_dir,
      "factor_solutions",
      paste0("qmethod_safe_summary_", solution_label, ".txt")
    )
  )

  ###########################################################################
  # 13.9 Solution diagnostics
  ###########################################################################

  # Eigenvalues can be approximated from rotated loadings as column sums of
  # squared factor loadings. This is mainly a comparative diagnostic here.
  eigenvalues <- colSums(as.matrix(loadings_numeric)^2, na.rm = TRUE)
  variance_percent <- 100 * eigenvalues / ncol(q_matrix)
  total_variance_percent <- sum(variance_percent, na.rm = TRUE)

  n_defining_total <- sum(defining_sorts$status == "defining", na.rm = TRUE)
  n_confounded <- sum(defining_sorts$status == "confounded", na.rm = TRUE)
  n_non_defining <- sum(defining_sorts$status == "non_defining", na.rm = TRUE)

  factor_names <- paste0("F", seq_len(nf))

  defining_per_factor <- sapply(
    factor_names,
    function(f) {
      sum(
        defining_sorts$status == "defining" &
          defining_sorts$defining_factor == f,
        na.rm = TRUE
      )
    }
  )

  factors_with_two_or_more_defining <- sum(defining_per_factor >= 2)

  warning_flags <- character(0)

  if (nf > 1 && any(defining_per_factor == 0)) {
    warning_flags <- c(warning_flags, "at_least_one_factor_without_defining_sorts")
  }

  if (nf > 1 && any(defining_per_factor == 1)) {
    warning_flags <- c(warning_flags, "single_sort_factor_present")
  }

  if (n_confounded > 0) {
    warning_flags <- c(warning_flags, "confounded_sorts_present")
  }

  if (n_non_defining > 0) {
    warning_flags <- c(warning_flags, "non_defining_sorts_present")
  }

  if (nf == 4) {
    warning_flags <- c(warning_flags, "four_factor_solution_cautious_with_only_8_participants")
  }

  if (nf == 1 && n_defining_total >= 6) {
    warning_flags <- c(warning_flags, "one_factor_solution_plausible_if_interpretively_coherent")
  }

  if (length(warning_flags) == 0) {
    warning_flags <- "none"
  }

  solution_summaries[[solution_label]] <- data.frame(
    solution = solution_label,
    nfactors = nf,
    extraction = extraction_method,
    rotation_argument = rotation_method,
    rotation_note = ifelse(
      nf == 1,
      "Rotation has no substantive role for a one-factor solution.",
      "Varimax rotation used."
    ),
    total_variance_percent = round(total_variance_percent, 2),
    eigenvalues = paste(round(eigenvalues, 3), collapse = "; "),
    variance_percent_by_factor = paste(round(variance_percent, 2), collapse = "; "),
    defining_sorts_total = n_defining_total,
    confounded_sorts = n_confounded,
    non_defining_sorts = n_non_defining,
    defining_sorts_by_factor = paste(
      paste0(names(defining_per_factor), "=", defining_per_factor),
      collapse = "; "
    ),
    factors_with_two_or_more_defining_sorts = factors_with_two_or_more_defining,
    distinguishing_statements = ifelse(
      is.na(n_distinguishing),
      "NA",
      as.character(n_distinguishing)
    ),
    consensus_statements = ifelse(
      is.na(n_consensus),
      "NA",
      as.character(n_consensus)
    ),
    warning_flags = paste(warning_flags, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

###############################################################################
# 14. EXPORT SOLUTION COMPARISON TABLE
###############################################################################

solution_comparison <- do.call(rbind, solution_summaries)

write_csv_utf8(
  solution_comparison,
  file.path(output_dir, "reports", "solution_comparison_table.csv")
)

###############################################################################
# 15. GENERATE STRICT LATEX REPORT
###############################################################################

distribution_latex_table <- data.frame(
  score = allowed_scores,
  expected_count = expected_counts,
  stringsAsFactors = FALSE
)

validation_latex_table <- validation_report[, c("check", "status"), drop = FALSE]

participant_status_latex_table <- participant_distribution_status

solution_latex_table <- solution_comparison[, c(
  "solution",
  "total_variance_percent",
  "defining_sorts_total",
  "confounded_sorts",
  "non_defining_sorts",
  "defining_sorts_by_factor",
  "distinguishing_statements",
  "consensus_statements",
  "warning_flags"
), drop = FALSE]

generated_files <- data.frame(
  file = c(
    file.path(output_dir, "validation", "input_validation_report.csv"),
    file.path(output_dir, "validation", "input_validation_report.txt"),
    file.path(output_dir, "validation", "qsort_distribution_check.csv"),
    file.path(output_dir, "validation", "participant_distribution_status.csv"),
    file.path(output_dir, "validation", "statement_metadata.csv"),
    file.path(output_dir, "validation", "q_matrix_numeric.csv"),
    file.path(output_dir, "correlations", "participant_correlation_matrix.csv"),
    file.path(output_dir, "factor_solutions", "loadings_1F.csv"),
    file.path(output_dir, "factor_solutions", "loadings_2F.csv"),
    file.path(output_dir, "factor_solutions", "loadings_3F.csv"),
    file.path(output_dir, "factor_solutions", "loadings_4F.csv"),
    file.path(output_dir, "factor_solutions", "defining_sorts_1F.csv"),
    file.path(output_dir, "factor_solutions", "defining_sorts_2F.csv"),
    file.path(output_dir, "factor_solutions", "defining_sorts_3F.csv"),
    file.path(output_dir, "factor_solutions", "defining_sorts_4F.csv"),
    file.path(output_dir, "factor_arrays", "factor_arrays_1F.csv"),
    file.path(output_dir, "factor_arrays", "factor_arrays_2F.csv"),
    file.path(output_dir, "factor_arrays", "factor_arrays_3F.csv"),
    file.path(output_dir, "factor_arrays", "factor_arrays_4F.csv"),
    file.path(output_dir, "distinguishing_consensus", "distinguishing_consensus_1F_not_applicable.txt"),
    file.path(output_dir, "distinguishing_consensus", "distinguishing_consensus_2F.csv"),
    file.path(output_dir, "distinguishing_consensus", "distinguishing_consensus_3F.csv"),
    file.path(output_dir, "distinguishing_consensus", "distinguishing_consensus_4F.csv"),
    file.path(output_dir, "reports", "solution_comparison_table.csv")
  ),
  stringsAsFactors = FALSE
)

latex_lines <- c(
  "\\documentclass[11pt]{article}",
  "\\usepackage[utf8]{inputenc}",
  "\\usepackage[T1]{fontenc}",
  "\\usepackage[margin=1in]{geometry}",
  "\\usepackage{array}",
  "\\usepackage{longtable}",
  "\\usepackage{booktabs}",
  "\\usepackage{hyperref}",
  "\\setlength{\\parindent}{0pt}",
  "\\setlength{\\parskip}{6pt}",
  "",
  "\\title{Q-Methodology Solution Comparison Report}",
  "\\author{Automated qmethod workflow}",
  paste0("\\date{", latex_escape(as.character(Sys.Date())), "}"),
  "",
  "\\begin{document}",
  "\\maketitle",
  "",
  "\\section*{Input summary}",
  paste0("Input file: \\texttt{", latex_escape(input_file), "}\\\\"),
  paste0("Number of statements: ", expected_n_statements, "\\\\"),
  paste0("Number of Q-sorts: ", length(expected_participants), "\\\\"),
  paste0("Participants: \\texttt{", latex_escape(paste(expected_participants, collapse = ", ")), "}\\\\"),
  paste0("Extraction method: \\texttt{", latex_escape(extraction_method), "}\\\\"),
  paste0("Rotation argument passed to \\texttt{qmethod}: \\texttt{", latex_escape(rotation_method), "}\\\\"),
  paste0("Correlation method: \\texttt{", latex_escape(correlation_method), "}\\\\"),
  paste0("Significant loading threshold used for fallback diagnostics: $\\pm ", round(loading_significance_threshold, 3), "$\\\\"),
  "",
  "\\section*{Validation summary}",
  latex_table(
    validation_latex_table,
    caption = "Input validation checks.",
    label = "tab:input-validation",
    align = "ll"
  ),
  "",
  "\\section*{Forced Q-sort distribution}",
  latex_table(
    distribution_latex_table,
    caption = "Expected forced distribution for the 49-statement Q-sort.",
    label = "tab:forced-distribution",
    align = "rr"
  ),
  "",
  "\\section*{Participant distribution status}",
  latex_table(
    participant_status_latex_table,
    caption = "Forced-distribution validation status by participant.",
    label = "tab:participant-distribution-status",
    align = "ll"
  ),
  "",
  "\\section*{Solution comparison}",
  latex_table(
    solution_latex_table,
    caption = "Comparison of 1-, 2-, 3-, and 4-factor qmethod solutions.",
    label = "tab:solution-comparison",
    align = "lrrrrllll"
  ),
  "",
  "\\section*{Interpretive caution}",
  "The table above is a diagnostic aid. The final retained solution should not be selected only by explained variance. It should also be evaluated by the number of defining Q-sorts per factor, the presence of confounded or non-defining sorts, the coherence of factor arrays, and the substantive interpretability of distinguishing and consensus statements.",
  "",
  "A one-factor solution is acceptable when it represents a coherent dominant shared viewpoint and additional factors do not add stable or interpretable perspectives. Multi-factor solutions are preferable only when each retained factor is statistically and substantively defensible.",
  "",
  "For the one-factor solution, distinguishing and consensus statements are not applicable because there are no between-factor comparisons.",
  "",
  "\\section*{Generated output files}",
  "\\begin{longtable}{p{0.95\\textwidth}}",
  "\\toprule",
  "\\textbf{File} \\\\",
  "\\midrule",
  paste0(latex_escape(generated_files$file), " \\\\", collapse = "\n"),
  "\\bottomrule",
  "\\end{longtable}",
  "",
  "\\end{document}"
)

latex_report_path <- file.path(output_dir, "reports", "solution_comparison_report.tex")

writeLines(
  latex_lines,
  con = latex_report_path,
  useBytes = TRUE
)

###############################################################################
# 16. FINAL MESSAGE
###############################################################################

message("Q-methodology workflow completed successfully.")
message("Solution comparison CSV: ", file.path(output_dir, "reports", "solution_comparison_table.csv"))
message("Strict LaTeX report: ", latex_report_path)

  
