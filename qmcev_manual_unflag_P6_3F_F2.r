#!/usr/bin/env Rscript

###############################################################################
# Transparent, reproducible Q-methodology workflow around the R package qmethod
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
#   -6  -5  -5  -4  -4  -4  -3  -3  -3  -3
#   -2  -2  -2  -2  -2  -1  -1  -1  -1  -1  -1
#    0   0   0   0   0   0   0
#    1   1   1   1   1   1
#    2   2   2   2   2
#    3   3   3   3
#    4   4   4
#    5   5
#    6
#
# Important methodological patch:
#
#   In the automatic 3-factor solution, P6 was flagged as defining F2 with a
#   negative loading. For interpretive clarity, this script applies a documented
#   manual-flagging override:
#
#       P6 is manually unflagged from F2 in the 3-factor solution.
#
#   The script then recalculates the 3F factor z-scores, rounded factor scores,
#   factor characteristics, and distinguishing/consensus statements using the
#   revised flag matrix through qmethod::qzscores().
#
#   The extraction and rotation are NOT rerun. Only the post-rotation flagging
#   and downstream statement calculations are adjusted.
#
###############################################################################

###############################################################################
# 0. USER SETTINGS
###############################################################################

# Path to the input CSV file.
input_file <- "data/qsort_data.csv"

# Main output directory.
output_dir <- "outputs_manual_unflag_P6_3F_F2"

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
extraction_method <- "PCA"
rotation_method <- "varimax"
correlation_method <- "pearson"

# Significant loading threshold.
#
# For 49 statements:
#   SE = 1 / sqrt(49) = 0.143
#   p < 0.01 threshold = 2.58 / sqrt(49) = 0.369
loading_significance_threshold <- 2.58 / sqrt(expected_n_statements)

# Manual-flagging override switch.
#
# Set to TRUE to unflag P6 from F2 in the 3-factor solution.
apply_manual_unflag_P6_F2_in_3F <- TRUE

# Manual-flagging label used in output filenames.
manual_3F_suffix <- "manual_unflag_P6_F2"

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

  return(out)
}

# Extract a component from a qmethod result.
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

  result <- do.call(qmethod::qmethod, args)

  return(result)
}

# Convert qmethod's flagged matrix to a logical matrix.
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

  storage.mode(out) <- "logical"

  rownames(out) <- row_names
  colnames(out) <- col_names

  return(out)
}

# Fallback flagging based on significant loadings.
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

# Apply project-specific manual flagging overrides.
#
# Current override:
#   In the 3F solution, P6 is manually unflagged from F2 because it was
#   automatically flagged with a negative loading. This supports interpretation
#   of factors as shared positive viewpoints.
apply_manual_flag_overrides <- function(flagged_matrix, solution_label) {
  manual_notes <- character(0)

  if (apply_manual_unflag_P6_F2_in_3F && solution_label == "3F") {
    if ("P6" %in% rownames(flagged_matrix) && "F2" %in% colnames(flagged_matrix)) {
      old_value <- flagged_matrix["P6", "F2"]
      flagged_matrix["P6", "F2"] <- FALSE

      manual_notes <- c(
        manual_notes,
        paste0(
          "Manual override applied in 3F: P6/F2 flag changed from ",
          old_value,
          " to FALSE. P6 is treated as non-defining for F2."
        )
      )
    } else {
      manual_notes <- c(
        manual_notes,
        "Manual override requested in 3F, but P6 or F2 was not found in the flag matrix."
      )
    }
  }

  return(list(
    flagged_matrix = flagged_matrix,
    manual_notes = manual_notes
  ))
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

# Recalculate z-scores and factor scores using the current flag matrix.
#
# This is essential for manual flagging:
#   If flags change, factor arrays must be recalculated from qzscores().
recalculate_qzscores_with_flags <- function(q_matrix, nfactors, loadings_numeric, flagged_matrix) {
  result_z <- qmethod::qzscores(
    dataset = as.data.frame(q_matrix),
    nfactors = nfactors,
    loa = as.data.frame(loadings_numeric),
    flagged = as.data.frame(flagged_matrix),
    forced = TRUE
  )

  return(result_z)
}

# Assign forced-distribution scores from z-scores.
#
# This is only a fallback if qmethod does not provide rounded scores.
assign_forced_scores_from_z <- function(z_vector, distribution_vector) {
  scores_sorted <- sort(distribution_vector, decreasing = TRUE)
  idx <- order(z_vector, decreasing = TRUE, na.last = TRUE)

  out <- rep(NA_real_, length(z_vector))
  out[idx] <- scores_sorted

  return(out)
}

# Build a clean factor array table from qmethod or qzscores output.
standardize_factor_array <- function(result, metadata, nfactors, distribution_vector) {
  zsc_raw <- get_component(result, "zsc", index = 5)

  z_df <- select_numeric_factor_columns(zsc_raw, nfactors, prefix = "F")

  if (is.null(z_df)) {
    stop(
      paste0("Could not extract statement z-scores for ", nfactors, "-factor solution."),
      call. = FALSE
    )
  }

  names(z_df) <- paste0("F", seq_len(nfactors), "_z")

  zsc_n_raw <- get_component(result, "zsc_n", index = 6)
  score_df <- select_numeric_factor_columns(zsc_n_raw, nfactors, prefix = "F")

  if (!is.null(score_df)) {
    names(score_df) <- paste0("F", seq_len(nfactors), "_score")
  }

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

# Calculate or extract distinguishing/consensus statements.
#
# For nfactors == 1:
#   Not applicable.
#
# For nfactors >= 2:
#   Recalculate using qmethod::qdc() from the z-scores and factor characteristics
#   after manual flagging.
extract_or_calculate_qdc_table <- function(q_matrix, result_for_statements, nfactors) {
  if (nfactors < 2) {
    return(NULL)
  }

  zsc <- get_component(result_for_statements, "zsc", index = 5)
  f_char <- get_component(result_for_statements, "f_char", index = 7)

  if (is.null(zsc) || is.null(f_char)) {
    return(NULL)
  }

  sed <- NULL

  if (is.list(f_char) && !is.null(f_char$sd_dif)) {
    sed <- f_char$sd_dif
  } else if (is.list(f_char) && length(f_char) >= 3) {
    sed <- f_char[[3]]
  }

  if (is.null(sed)) {
    return(NULL)
  }

  qdc_df <- tryCatch(
    {
      as.data.frame(
        qmethod::qdc(
          dataset = as.data.frame(q_matrix),
          nfactors = nfactors,
          zsc = zsc,
          sed = sed
        ),
        stringsAsFactors = FALSE
      )
    },
    error = function(e) {
      NULL
    }
  )

  return(qdc_df)
}

# Count distinguishing and consensus statements from a qdc table.
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
write_safe_qmethod_summary <- function(
  result_original,
  result_for_statements,
  nfactors,
  solution_label,
  loadings_export,
  defining_sorts,
  factor_array,
  manual_notes,
  path
) {
  lines <- character(0)

  lines <- c(lines, paste0("qmethod solution: ", solution_label))
  lines <- c(lines, paste0("Number of factors: ", nfactors))
  lines <- c(lines, paste0("Extraction method: ", extraction_method))
  lines <- c(lines, paste0("Rotation argument passed to qmethod: ", rotation_method))

  if (length(manual_notes) > 0) {
    lines <- c(lines, "", "Manual flagging notes:")
    lines <- c(lines, manual_notes)
  } else {
    lines <- c(lines, "", "Manual flagging notes:")
    lines <- c(lines, "No manual flagging overrides applied.")
  }

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
      "Distinguishing and consensus statements are applicable if qmethod::qdc() returned a table."
    )
  }

  lines <- c(lines, "", "Available original qmethod result components:")
  lines <- c(lines, paste(names(result_original), collapse = ", "))

  lines <- c(lines, "", "Loadings and flags:")
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

  f_char <- get_component(result_for_statements, "f_char", index = 7)

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

if (any(is.na(metadata$statement_id) | metadata$statement_id == "")) {
  validation_pass <- FALSE
  validation_messages <- c(validation_messages, "Missing statement_id value(s).")
}

if (any(duplicated(metadata$statement_id))) {
  validation_pass <- FALSE
  validation_messages <- c(validation_messages, "Duplicated statement_id value(s).")
}

if (any(is.na(metadata$statement_text) | metadata$statement_text == "")) {
  validation_pass <- FALSE
  validation_messages <- c(validation_messages, "Missing statement_text value(s).")
}

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

non_numeric_cells <- is.na(q_data_numeric) & !is.na(q_data_raw)

if (any(non_numeric_cells)) {
  validation_pass <- FALSE
  validation_messages <- c(
    validation_messages,
    "Some Q-sort cells could not be converted to numeric values."
  )
}

if (any(is.na(q_data_numeric))) {
  validation_pass <- FALSE
  validation_messages <- c(
    validation_messages,
    "Missing Q-sort score(s) detected."
  )
}

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

distribution_check$status_by_score <- "PASS"

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
  # 13.1 Run qmethod normally
  ###########################################################################

  result_original <- tryCatch(
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

  all_results[[solution_label]] <- result_original

  saveRDS(
    result_original,
    file = file.path(output_dir, "rds", paste0("qmethod_result_", solution_label, "_automatic.rds"))
  )

  ###########################################################################
  # 13.2 Extract loadings from the original qmethod run
  ###########################################################################

  loa_raw <- get_component(result_original, "loa", index = 3)

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
  # 13.3 Extract automatic flags
  ###########################################################################

  flagged_raw <- get_component(result_original, "flagged", index = 4)

  flagged_matrix_automatic <- as_logical_flag_matrix(
    x = flagged_raw,
    participants = expected_participants,
    nfactors = nf
  )

  if (all(flagged_matrix_automatic == FALSE)) {
    flagged_matrix_automatic <- fallback_flags_from_loadings(
      loadings_numeric = loadings_numeric,
      threshold = loading_significance_threshold
    )
  }

  ###########################################################################
  # 13.4 Apply manual flagging overrides
  ###########################################################################

  override_result <- apply_manual_flag_overrides(
    flagged_matrix = flagged_matrix_automatic,
    solution_label = solution_label
  )

  flagged_matrix <- override_result$flagged_matrix
  manual_notes <- override_result$manual_notes

  manual_override_applied <- length(manual_notes) > 0 &&
    any(grepl("Manual override applied", manual_notes))

  ###########################################################################
  # 13.5 Recalculate z-scores and factor scores if manual flags changed
  ###########################################################################

  if (manual_override_applied) {
    message("Manual flagging override applied for ", solution_label, ". Recalculating z-scores and factor scores.")

    result_for_statements <- recalculate_qzscores_with_flags(
      q_matrix = q_matrix,
      nfactors = nf,
      loadings_numeric = loadings_numeric,
      flagged_matrix = flagged_matrix
    )

    saveRDS(
      result_for_statements,
      file = file.path(output_dir, "rds", paste0("qmethod_result_", solution_label, "_", manual_3F_suffix, ".rds"))
    )
  } else {
    result_for_statements <- result_original
  }

  ###########################################################################
  # 13.6 Export loadings and final flags
  ###########################################################################

  flagged_df <- as.data.frame(flagged_matrix, stringsAsFactors = FALSE)
  names(flagged_df) <- paste0("flag_", names(flagged_df))

  flagged_auto_df <- as.data.frame(flagged_matrix_automatic, stringsAsFactors = FALSE)
  names(flagged_auto_df) <- paste0("auto_flag_", names(flagged_auto_df))

  loadings_export <- data.frame(
    participant = expected_participants,
    loadings_numeric,
    flagged_auto_df,
    flagged_df,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  loadings_filename <- if (manual_override_applied) {
    paste0("loadings_", solution_label, "_", manual_3F_suffix, ".csv")
  } else {
    paste0("loadings_", solution_label, ".csv")
  }

  write_csv_utf8(
    loadings_export,
    file.path(output_dir, "factor_solutions", loadings_filename)
  )

  # Also preserve automatic loadings/flags for transparency if manual override was used.
  if (manual_override_applied) {
    loadings_export_auto_only <- data.frame(
      participant = expected_participants,
      loadings_numeric,
      flagged_auto_df,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    write_csv_utf8(
      loadings_export_auto_only,
      file.path(output_dir, "factor_solutions", paste0("loadings_", solution_label, "_automatic.csv"))
    )
  }

  ###########################################################################
  # 13.7 Classify defining, confounded, and non-defining sorts
  ###########################################################################

  defining_sorts <- classify_defining_sorts(
    loadings_numeric = loadings_numeric,
    flagged_matrix = flagged_matrix
  )

  defining_filename <- if (manual_override_applied) {
    paste0("defining_sorts_", solution_label, "_", manual_3F_suffix, ".csv")
  } else {
    paste0("defining_sorts_", solution_label, ".csv")
  }

  write_csv_utf8(
    defining_sorts,
    file.path(output_dir, "factor_solutions", defining_filename)
  )

  if (manual_override_applied) {
    defining_sorts_auto <- classify_defining_sorts(
      loadings_numeric = loadings_numeric,
      flagged_matrix = flagged_matrix_automatic
    )

    write_csv_utf8(
      defining_sorts_auto,
      file.path(output_dir, "factor_solutions", paste0("defining_sorts_", solution_label, "_automatic.csv"))
    )
  }

  ###########################################################################
  # 13.8 Export factor characteristics
  ###########################################################################

  f_char <- get_component(result_for_statements, "f_char", index = 7)

  fchar_filename <- if (manual_override_applied) {
    paste0("factor_characteristics_", solution_label, "_", manual_3F_suffix, ".txt")
  } else {
    paste0("factor_characteristics_", solution_label, ".txt")
  }

  if (!is.null(f_char)) {
    writeLines(
      safe_capture_print(f_char),
      con = file.path(output_dir, "factor_solutions", fchar_filename),
      useBytes = TRUE
    )
  }

  ###########################################################################
  # 13.9 Build and export factor arrays
  ###########################################################################

  factor_array <- standardize_factor_array(
    result = result_for_statements,
    metadata = metadata,
    nfactors = nf,
    distribution_vector = expected_distribution
  )

  factor_array_filename <- if (manual_override_applied) {
    paste0("factor_arrays_", solution_label, "_", manual_3F_suffix, ".csv")
  } else {
    paste0("factor_arrays_", solution_label, ".csv")
  }

  write_csv_utf8(
    factor_array,
    file.path(output_dir, "factor_arrays", factor_array_filename)
  )

  for (j in seq_len(nf)) {
    f_name <- paste0("F", j)
    score_col <- paste0(f_name, "_score")
    z_col <- paste0(f_name, "_z")

    ranked_array <- factor_array[
      order(factor_array[[score_col]], factor_array[[z_col]], decreasing = TRUE),
      c("statement_id", "statement_text", "category", z_col, score_col),
      drop = FALSE
    ]

    ranked_filename <- if (manual_override_applied) {
      paste0("factor_array_", solution_label, "_", manual_3F_suffix, "_", f_name, "_ranked.csv")
    } else {
      paste0("factor_array_", solution_label, "_", f_name, "_ranked.csv")
    }

    write_csv_utf8(
      ranked_array,
      file.path(output_dir, "factor_arrays", ranked_filename)
    )
  }

  ###########################################################################
  # 13.10 Export distinguishing and consensus statements
  ###########################################################################

  qdc_export <- NULL
  n_distinguishing <- NA_integer_
  n_consensus <- NA_integer_

  if (nf >= 2) {
    qdc_table <- extract_or_calculate_qdc_table(
      q_matrix = q_matrix,
      result_for_statements = result_for_statements,
      nfactors = nf
    )

    if (!is.null(qdc_table)) {
      if (nrow(qdc_table) == nrow(metadata)) {
        qdc_export <- cbind(metadata, qdc_table, stringsAsFactors = FALSE)
      } else {
        qdc_export <- qdc_table
      }

      qdc_filename <- if (manual_override_applied) {
        paste0("distinguishing_consensus_", solution_label, "_", manual_3F_suffix, ".csv")
      } else {
        paste0("distinguishing_consensus_", solution_label, ".csv")
      }

      write_csv_utf8(
        qdc_export,
        file.path(output_dir, "distinguishing_consensus", qdc_filename)
      )

      qdc_counts <- count_qdc_labels(qdc_export)
      n_distinguishing <- qdc_counts$distinguishing
      n_consensus <- qdc_counts$consensus
    } else {
      qdc_not_available_filename <- if (manual_override_applied) {
        paste0("distinguishing_consensus_", solution_label, "_", manual_3F_suffix, "_not_available.txt")
      } else {
        paste0("distinguishing_consensus_", solution_label, "_not_available.txt")
      }

      writeLines(
        "Distinguishing / consensus table could not be extracted or recalculated.",
        con = file.path(output_dir, "distinguishing_consensus", qdc_not_available_filename),
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
  # 13.11 Write safe text summary for this qmethod solution
  ###########################################################################

  safe_summary_filename <- if (manual_override_applied) {
    paste0("qmethod_safe_summary_", solution_label, "_", manual_3F_suffix, ".txt")
  } else {
    paste0("qmethod_safe_summary_", solution_label, ".txt")
  }

  write_safe_qmethod_summary(
    result_original = result_original,
    result_for_statements = result_for_statements,
    nfactors = nf,
    solution_label = solution_label,
    loadings_export = loadings_export,
    defining_sorts = defining_sorts,
    factor_array = factor_array,
    manual_notes = manual_notes,
    path = file.path(output_dir, "factor_solutions", safe_summary_filename)
  )

  ###########################################################################
  # 13.12 Solution diagnostics
  ###########################################################################

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

  if (manual_override_applied) {
    warning_flags <- c(warning_flags, "manual_flagging_override_applied")
  }

  if (length(warning_flags) == 0) {
    warning_flags <- "none"
  }

  solution_name_for_report <- if (manual_override_applied) {
    paste0(solution_label, "_", manual_3F_suffix)
  } else {
    solution_label
  }

  solution_summaries[[solution_name_for_report]] <- data.frame(
    solution = solution_name_for_report,
    nfactors = nf,
    extraction = extraction_method,
    rotation_argument = rotation_method,
    rotation_note = ifelse(
      nf == 1,
      "Rotation has no substantive role for a one-factor solution.",
      "Varimax rotation used."
    ),
    manual_flagging = ifelse(
      manual_override_applied,
      paste(manual_notes, collapse = " | "),
      "none"
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
  "manual_flagging",
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
    file.path(output_dir, "factor_solutions", "loadings_3F_manual_unflag_P6_F2.csv"),
    file.path(output_dir, "factor_solutions", "loadings_3F_automatic.csv"),
    file.path(output_dir, "factor_solutions", "loadings_4F.csv"),
    file.path(output_dir, "factor_solutions", "defining_sorts_1F.csv"),
    file.path(output_dir, "factor_solutions", "defining_sorts_2F.csv"),
    file.path(output_dir, "factor_solutions", "defining_sorts_3F_manual_unflag_P6_F2.csv"),
    file.path(output_dir, "factor_solutions", "defining_sorts_3F_automatic.csv"),
    file.path(output_dir, "factor_solutions", "defining_sorts_4F.csv"),
    file.path(output_dir, "factor_arrays", "factor_arrays_1F.csv"),
    file.path(output_dir, "factor_arrays", "factor_arrays_2F.csv"),
    file.path(output_dir, "factor_arrays_3F_manual_unflag_P6_F2.csv"),
    file.path(output_dir, "factor_arrays", "factor_arrays_4F.csv"),
    file.path(output_dir, "distinguishing_consensus", "distinguishing_consensus_1F_not_applicable.txt"),
    file.path(output_dir, "distinguishing_consensus", "distinguishing_consensus_2F.csv"),
    file.path(output_dir, "distinguishing_consensus", "distinguishing_consensus_3F_manual_unflag_P6_F2.csv"),
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
  "\\section*{Manual flagging note}",
  "For the 3-factor solution, P6 was automatically flagged as defining F2 with a negative loading. The manually adjusted workflow unflags P6 from F2 and treats P6 as non-defining. The 3F factor arrays, factor characteristics, and distinguishing/consensus statements are recalculated from the revised flag matrix.",
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
    caption = "Comparison of 1-, 2-, 3-, and 4-factor qmethod solutions. The 3F solution is reported after the documented manual flagging override.",
    label = "tab:solution-comparison",
    align = "lrrrrlllll"
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
message("Manual unflagging P6/F2 in 3F enabled: ", apply_manual_unflag_P6_F2_in_3F)
message("Solution comparison CSV: ", file.path(output_dir, "reports", "solution_comparison_table.csv"))
message("Strict LaTeX report: ", latex_report_path)

