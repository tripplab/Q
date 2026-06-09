#!/usr/bin/env Rscript

###############################################################################
# MCEV Project-specific Q-methodology workflow around qmethod by trippm@tripplab.com
#
# Workflow:
# 1. Read CSV
# 2. Validate 49 statements and forced distribution
# 3. Separate metadata from numeric Q-sort matrix
# 4. Run qmethod for 1F, 2F, 3F, and 4F
# 5. Export loadings
# 6. Export factor arrays
# 7. Export distinguishing / consensus statements where applicable
# 8. Generate a strict LaTeX solution-comparison report
#
# Expected input CSV columns:
# statement_id, statement_text, category, P1, P2, P3, P4, P5, P6, P7, P8
###############################################################################

###############################################################################
# 0. USER SETTINGS
###############################################################################

input_file <- "Qsort_results.csv"
output_dir <- "qresults"

expected_n_statements <- 49L
expected_participants <- paste0("P", 1:8)

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

factor_solutions <- 1:4
extraction_method <- "PCA"
rotation_method_multi_factor <- "varimax"
correlation_method <- "pearson"

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
      paste(sprintf('\"%s\"', missing_packages), collapse = ", "),
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
# 3. HELPER FUNCTIONS
###############################################################################

write_csv <- function(x, path, row.names = FALSE) {
  utils::write.csv(
    x,
    file = path,
    row.names = row.names,
    fileEncoding = "UTF-8",
    na = ""
  )
}

get_component <- function(result, name, index) {
  if (!is.null(result[[name]])) {
    return(result[[name]])
  }
  if (length(result) >= index) {
    return(result[[index]])
  }
  return(NULL)
}

safe_as_data_frame <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  as.data.frame(x, stringsAsFactors = FALSE)
}

select_numeric_factor_columns <- function(x, nfactors) {
  x <- safe_as_data_frame(x)
  if (is.null(x)) {
    return(NULL)
  }

  numeric_cols <- vapply(x, is.numeric, logical(1))

  if (sum(numeric_cols) < nfactors) {
    stop(
      "Could not identify enough numeric factor columns in qmethod output.",
      call. = FALSE
    )
  }

  out <- x[, which(numeric_cols)[seq_len(nfactors)], drop = FALSE]
  names(out) <- paste0("F", seq_len(nfactors))
  out
}

as_logical_matrix <- function(x, nrow_expected, ncol_expected, row_names, col_names) {
  if (is.null(x)) {
    out <- matrix(
      FALSE,
      nrow = nrow_expected,
      ncol = ncol_expected,
      dimnames = list(row_names, col_names)
    )
    return(out)
  }

  x <- as.data.frame(x, stringsAsFactors = FALSE)

  if (ncol(x) < ncol_expected) {
    out <- matrix(
      FALSE,
      nrow = nrow_expected,
      ncol = ncol_expected,
      dimnames = list(row_names, col_names)
    )
    return(out)
  }

  x <- x[, seq_len(ncol_expected), drop = FALSE]
  out <- as.matrix(x)

  storage.mode(out) <- "logical"

  rownames(out) <- row_names
  colnames(out) <- col_names

  out
}

fallback_flags_from_loadings <- function(loadings_numeric, threshold) {
  abs_loadings <- abs(as.matrix(loadings_numeric))
  significant <- abs_loadings >= threshold

  out <- matrix(
    FALSE,
    nrow = nrow(significant),
    ncol = ncol(significant),
    dimnames = dimnames(significant)
  )

  for (i in seq_len(nrow(significant))) {
    sig_idx <- which(significant[i, ])

    if (length(sig_idx) == 1) {
      out[i, sig_idx] <- TRUE
    } else if (length(sig_idx) > 1) {
      max_idx <- which.max(abs_loadings[i, ])
      out[i, sig_idx] <- TRUE

      if (length(sig_idx) > 1) {
        out[i, sig_idx] <- TRUE
      }
    }
  }

  out
}

classify_defining_sorts <- function(loadings_numeric, flagged_matrix) {
  participant <- rownames(loadings_numeric)
  factor_names <- colnames(loadings_numeric)

  flag_count <- rowSums(flagged_matrix, na.rm = TRUE)

  status <- ifelse(
    flag_count == 0,
    "non_defining",
    ifelse(flag_count == 1, "defining", "confounded")
  )

  defining_factor <- rep(NA_character_, length(participant))
  defining_loading <- rep(NA_real_, length(participant))
  loading_sign <- rep(NA_character_, length(participant))

  for (i in seq_along(participant)) {
    idx <- which(flagged_matrix[i, ])

    if (length(idx) == 1) {
      defining_factor[i] <- factor_names[idx]
      defining_loading[i] <- loadings_numeric[i, idx]
      loading_sign[i] <- ifelse(defining_loading[i] >= 0, "positive", "negative")
    } else if (length(idx) > 1) {
      defining_factor[i] <- paste(factor_names[idx], collapse = ";")
      defining_loading[i] <- NA_real_
      loading_sign[i] <- "confounded"
    }
  }

  data.frame(
    participant = participant,
    status = status,
    defining_factor = defining_factor,
    defining_loading = defining_loading,
    loading_sign = loading_sign,
    stringsAsFactors = FALSE
  )
}

assign_forced_scores_from_z <- function(z_vector, distribution_vector) {
  scores_sorted <- sort(distribution_vector, decreasing = TRUE)
  idx <- order(z_vector, decreasing = TRUE, na.last = TRUE)

  out <- rep(NA_real_, length(z_vector))
  out[idx] <- scores_sorted

  out
}

standardize_factor_array <- function(result, metadata, nfactors, distribution_vector) {
  zsc <- get_component(result, "zsc", 5)
  zsc_n <- get_component(result, "zsc_n", 6)

  z_df <- select_numeric_factor_columns(zsc, nfactors)
  names(z_df) <- paste0("F", seq_len(nfactors), "_z")

  score_df <- NULL

  if (!is.null(zsc_n)) {
    candidate_scores <- safe_as_data_frame(zsc_n)
    numeric_cols <- vapply(candidate_scores, is.numeric, logical(1))

    if (sum(numeric_cols) >= nfactors) {
      score_df <- candidate_scores[, which(numeric_cols)[seq_len(nfactors)], drop = FALSE]
      names(score_df) <- paste0("F", seq_len(nfactors), "_score")
    }
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
      score_df[[j]] <- assign_forced_scores_from_z(z_df[[j]], distribution_vector)
    }
  }

  out <- cbind(
    metadata,
    z_df,
    score_df,
    stringsAsFactors = FALSE
  )

  out
}

extract_qdc_table <- function(result, q_matrix, nfactors) {
  if (nfactors < 2) {
    return(NULL)
  }

  qdc_table <- get_component(result, "qdc", 8)

  if (!is.null(qdc_table)) {
    return(as.data.frame(qdc_table, stringsAsFactors = FALSE))
  }

  zsc <- get_component(result, "zsc", 5)
  f_char <- get_component(result, "f_char", 7)

  sed <- NULL

  if (!is.null(f_char)) {
    if (is.list(f_char) && length(f_char) >= 3) {
      sed <- as.data.frame(f_char[[3]])
    }
  }

  if (!is.null(zsc) && !is.null(sed)) {
    return(as.data.frame(qmethod::qdc(q_matrix, nfactors, zsc, sed)))
  }

  return(NULL)
}

latex_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""

  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x <- gsub("~", "\\\\textasciitilde{}", x)

  x
}

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
    row <- paste(latex_escape(df[i, ]), collapse = " & ")
    lines <- c(lines, paste0(row, " \\\\"))
  }

  lines <- c(lines, "\\hline")
  lines <- c(lines, "\\end{tabular}")
  lines <- c(lines, "\\end{table}")

  paste(lines, collapse = "\n")
}

run_qmethod_solution <- function(q_matrix, nfactors) {
  rotation <- ifelse(nfactors == 1, "none", rotation_method_multi_factor)

  fmls <- names(formals(qmethod::qmethod))

  args <- list(
    dataset = as.data.frame(q_matrix),
    nfactors = nfactors,
    rotation = rotation,
    forced = TRUE,
    cor.method = correlation_method
  )

  if ("extraction" %in% fmls) {
    args$extraction <- extraction_method
  }

  if ("distribution" %in% fmls) {
    args$distribution <- NULL
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

  do.call(qmethod::qmethod, args)
}

###############################################################################
# 4. READ INPUT CSV
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
# 5. VALIDATE STRUCTURE
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

if ("statement_id" %in% names(raw_data)) {
  if (any(is.na(raw_data$statement_id) | raw_data$statement_id == "")) {
    validation_pass <- FALSE
    validation_messages <- c(validation_messages, "Missing statement_id value(s).")
  }

  if (any(duplicated(raw_data$statement_id))) {
    validation_pass <- FALSE
    validation_messages <- c(validation_messages, "Duplicated statement_id value(s).")
  }
}

if ("statement_text" %in% names(raw_data)) {
  if (any(is.na(raw_data$statement_text) | raw_data$statement_text == "")) {
    validation_pass <- FALSE
    validation_messages <- c(validation_messages, "Missing statement_text value(s).")
  }
}

if ("category" %in% names(raw_data)) {
  if (any(is.na(raw_data$category) | raw_data$category == "")) {
    validation_pass <- FALSE
    validation_messages <- c(validation_messages, "Missing category value(s).")
  }
}

if (!all(c(required_metadata_cols, expected_participants) %in% names(raw_data))) {
  validation_report <- data.frame(
    check = "CSV structure",
    status = "FAIL",
    message = paste(validation_messages, collapse = " | "),
    stringsAsFactors = FALSE
  )

  write_csv(
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
# 6. SEPARATE METADATA AND Q-SORT MATRIX
###############################################################################

metadata <- raw_data[, required_metadata_cols, drop = FALSE]

q_data_raw <- raw_data[, expected_participants, drop = FALSE]

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

integer_score_mask <- abs(as.matrix(q_data_numeric) - round(as.matrix(q_data_numeric))) > .Machine$double.eps^0.5

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
# 7. VALIDATE FORCED DISTRIBUTION
###############################################################################

expected_counts <- as.integer(table(factor(expected_distribution, levels = allowed_scores)))

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

distribution_check$status <- apply(
  distribution_check[, expected_participants, drop = FALSE],
  1,
  function(x) {
    if (all(x == distribution_check$expected)) {
      "PASS"
    } else {
      "FAIL"
    }
  }
)

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

write_csv(
  distribution_check,
  file.path(output_dir, "validation", "qsort_distribution_check.csv")
)

write_csv(
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

write_csv(
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
# 8. EXPORT PREPARED METADATA AND Q-MATRIX
###############################################################################

write_csv(
  metadata,
  file.path(output_dir, "validation", "statement_metadata.csv")
)

write_csv(
  data.frame(statement_id = rownames(q_matrix), q_matrix, check.names = FALSE),
  file.path(output_dir, "validation", "q_matrix_numeric.csv")
)

###############################################################################
# 9. PARTICIPANT CORRELATION MATRIX
###############################################################################

participant_correlation_matrix <- cor(
  q_matrix,
  method = correlation_method,
  use = "pairwise.complete.obs"
)

write_csv(
  data.frame(
    participant = rownames(participant_correlation_matrix),
    participant_correlation_matrix,
    check.names = FALSE
  ),
  file.path(output_dir, "correlations", "participant_correlation_matrix.csv")
)

###############################################################################
# 10. RUN QMETHOD SOLUTIONS
###############################################################################

solution_summaries <- list()
all_results <- list()

for (nf in factor_solutions) {
  solution_label <- paste0(nf, "F")

  message("Running qmethod solution: ", solution_label)

  result <- run_qmethod_solution(q_matrix, nfactors = nf)
  all_results[[solution_label]] <- result

  saveRDS(
    result,
    file = file.path(output_dir, "rds", paste0("qmethod_result_", solution_label, ".rds"))
  )

  captured_summary <- capture.output({
    print(result)
    cat("\n\n")
    suppressWarnings(print(summary(result)))
  })

  writeLines(
    captured_summary,
    con = file.path(output_dir, "factor_solutions", paste0("qmethod_summary_", solution_label, ".txt")),
    useBytes = TRUE
  )

  ###########################################################################
  # 10.1 EXPORT LOADINGS
  ###########################################################################

  loa_raw <- get_component(result, "loa", 3)
  loa_df <- safe_as_data_frame(loa_raw)

  if (is.null(loa_df)) {
    stop(paste0("Could not extract loadings for solution ", solution_label), call. = FALSE)
  }

  loadings_numeric <- select_numeric_factor_columns(loa_df, nf)
  rownames(loadings_numeric) <- expected_participants

  flagged_raw <- get_component(result, "flagged", 4)
  flagged_matrix <- as_logical_matrix(
    flagged_raw,
    nrow_expected = length(expected_participants),
    ncol_expected = nf,
    row_names = expected_participants,
    col_names = paste0("F", seq_len(nf))
  )

  if (all(flagged_matrix == FALSE)) {
    flagged_matrix <- fallback_flags_from_loadings(
      loadings_numeric,
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

  write_csv(
    loadings_export,
    file.path(output_dir, "factor_solutions", paste0("loadings_", solution_label, ".csv"))
  )

  ###########################################################################
  # 10.2 EXPORT DEFINING SORTS
  ###########################################################################

  defining_sorts <- classify_defining_sorts(
    loadings_numeric = loadings_numeric,
    flagged_matrix = flagged_matrix
  )

  write_csv(
    defining_sorts,
    file.path(output_dir, "factor_solutions", paste0("defining_sorts_", solution_label, ".csv"))
  )

  ###########################################################################
  # 10.3 EXPORT FACTOR CHARACTERISTICS IF AVAILABLE
  ###########################################################################

  f_char <- get_component(result, "f_char", 7)

  if (!is.null(f_char)) {
    f_char_capture <- capture.output(print(f_char))

    writeLines(
      f_char_capture,
      con = file.path(output_dir, "factor_solutions", paste0("factor_characteristics_", solution_label, ".txt")),
      useBytes = TRUE
    )
  }

  ###########################################################################
  # 10.4 EXPORT FACTOR ARRAYS
  ###########################################################################

  factor_array <- standardize_factor_array(
    result = result,
    metadata = metadata,
    nfactors = nf,
    distribution_vector = expected_distribution
  )

  write_csv(
    factor_array,
    file.path(output_dir, "factor_arrays", paste0("factor_arrays_", solution_label, ".csv"))
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

    write_csv(
      ranked_array,
      file.path(output_dir, "factor_arrays", paste0("factor_array_", solution_label, "_", f_name, "_ranked.csv"))
    )
  }

  ###########################################################################
  # 10.5 EXPORT DISTINGUISHING / CONSENSUS STATEMENTS
  ###########################################################################

  qdc_export <- NULL
  n_distinguishing <- NA_integer_
  n_consensus <- NA_integer_

  if (nf >= 2) {
    qdc_table <- extract_qdc_table(
      result = result,
      q_matrix = q_matrix,
      nfactors = nf
    )

    if (!is.null(qdc_table)) {
      if (nrow(qdc_table) == nrow(metadata)) {
        qdc_export <- cbind(metadata, qdc_table, stringsAsFactors = FALSE)
      } else {
        qdc_export <- qdc_table
      }

      write_csv(
        qdc_export,
        file.path(output_dir, "distinguishing_consensus", paste0("distinguishing_consensus_", solution_label, ".csv"))
      )

      if ("dist.and.cons" %in% names(qdc_export)) {
        n_distinguishing <- sum(
          grepl("Distinguishes", qdc_export$dist.and.cons),
          na.rm = TRUE
        )

        n_consensus <- sum(
          qdc_export$dist.and.cons == "Consensus",
          na.rm = TRUE
        )
      }
    } else {
      writeLines(
        "Distinguishing / consensus table could not be extracted from qmethod result.",
        con = file.path(output_dir, "distinguishing_consensus", paste0("distinguishing_consensus_", solution_label, "_not_available.txt")),
        useBytes = TRUE
      )
    }
  } else {
    writeLines(
      "Not applicable for a 1-factor solution.",
      con = file.path(output_dir, "distinguishing_consensus", "distinguishing_consensus_1F_not_applicable.txt"),
      useBytes = TRUE
    )
  }

  ###########################################################################
  # 10.6 SOLUTION DIAGNOSTICS
  ###########################################################################

  eigenvalues <- colSums(as.matrix(loadings_numeric)^2, na.rm = TRUE)
  variance_percent <- 100 * eigenvalues / ncol(q_matrix)
  total_variance_percent <- sum(variance_percent, na.rm = TRUE)

  unique_defining <- defining_sorts$status == "defining"
  n_defining_total <- sum(unique_defining, na.rm = TRUE)
  n_confounded <- sum(defining_sorts$status == "confounded", na.rm = TRUE)
  n_non_defining <- sum(defining_sorts$status == "non_defining", na.rm = TRUE)

  defining_per_factor <- sapply(
    paste0("F", seq_len(nf)),
    function(f) sum(defining_sorts$defining_factor == f & defining_sorts$status == "defining", na.rm = TRUE)
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
    rotation = ifelse(nf == 1, "none", rotation_method_multi_factor),
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
    distinguishing_statements = ifelse(is.na(n_distinguishing), "NA", as.character(n_distinguishing)),
    consensus_statements = ifelse(is.na(n_consensus), "NA", as.character(n_consensus)),
    warning_flags = paste(warning_flags, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

###############################################################################
# 11. EXPORT SOLUTION COMPARISON TABLE
###############################################################################

solution_comparison <- do.call(rbind, solution_summaries)

write_csv(
  solution_comparison,
  file.path(output_dir, "reports", "solution_comparison_table.csv")
)

###############################################################################
# 12. GENERATE STRICT LATEX REPORT
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
    file.path(output_dir, "validation", "qsort_distribution_check.csv"),
    file.path(output_dir, "validation", "q_matrix_numeric.csv"),
    file.path(output_dir, "correlations", "participant_correlation_matrix.csv"),
    file.path(output_dir, "factor_solutions", "loadings_1F.csv"),
    file.path(output_dir, "factor_solutions", "loadings_2F.csv"),
    file.path(output_dir, "factor_solutions", "loadings_3F.csv"),
    file.path(output_dir, "factor_solutions", "loadings_4F.csv"),
    file.path(output_dir, "factor_arrays", "factor_arrays_1F.csv"),
    file.path(output_dir, "factor_arrays", "factor_arrays_2F.csv"),
    file.path(output_dir, "factor_arrays", "factor_arrays_3F.csv"),
    file.path(output_dir, "factor_arrays", "factor_arrays_4F.csv"),
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
  paste0("Rotation for multi-factor solutions: \\texttt{", latex_escape(rotation_method_multi_factor), "}\\\\"),
  paste0("Correlation method: \\texttt{", latex_escape(correlation_method), "}\\\\"),
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
# 13. FINAL MESSAGE
###############################################################################

message("Q-methodology workflow completed successfully.")
message("Solution comparison CSV: ", file.path(output_dir, "reports", "solution_comparison_table.csv"))
message("Strict LaTeX report: ", latex_report_path)


  
