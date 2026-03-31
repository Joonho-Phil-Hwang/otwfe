# =============================================================================
# otwfe_file.R
#
# Large-scale CSV panel data analysis without loading the full file into memory.
# Uses the otwfe_init / otwfe_update / otwfe_finalize pipeline internally.
#
# Key properties:
#   - Full dataset is never materialized in R memory
#   - Adaptive chunked reading with automatic unit-boundary detection
#   - Automatic T_support detection via early-termination time-column scan
#   - Automatic time remapping to {1,...,T} (supports arbitrary integer / year values)
#   - Algebraically exact: matches plm to machine precision
#
# Prerequisite: file must be sorted by id_col (non-decreasing)
# =============================================================================

# =============================================================================
# Internal helper functions
# =============================================================================

# --------------------------------------------------------------------------
# .csv_col_names(): read column names from CSV header
# --------------------------------------------------------------------------
.csv_col_names <- function(path, sep) {
  names(data.table::fread(path, sep = sep, nrows = 0L, showProgress = FALSE))
}

# --------------------------------------------------------------------------
# .csv_read_chunk_seq(): read next n_rows rows from an open sequential connection
#
# con        : open text connection (file position maintained across calls)
# header_line: header row string (prepended so fread can infer column names)
# sep        : CSV delimiter
# n_rows     : maximum rows to read
# returns    : data.frame or NULL (EOF)
#
# Sequential reading avoids the O(skip) re-scan overhead of fread(skip=N),
# making total I/O O(n) instead of O(n^2 / chunk_size).
# --------------------------------------------------------------------------
.csv_read_chunk_seq <- function(con, header_line, sep, n_rows) {
  lines <- readLines(con, n = n_rows)
  if (length(lines) == 0L) return(NULL)
  chunk <- data.table::fread(
    input        = paste(c(header_line, lines), collapse = "\n"),
    sep          = sep,
    showProgress = FALSE
  )
  if (nrow(chunk) == 0L) return(NULL)
  as.data.frame(chunk)
}

# --------------------------------------------------------------------------
# .detect_time_levels(): detect unique calendar times via early-termination scan
#
# Optimization: maintains a file connection (no re-scanning) + early exit
#   - Reads scan_chunk lines at a time, extracts only the time column
#   - Stops after stable_rounds consecutive chunks with no new time values
#   - For sorted panels with small T, terminates after a few hundred thousand rows
#   - Worst case (time values concentrated at end of file): full scan
# --------------------------------------------------------------------------
.detect_time_levels <- function(path, time_col, sep, verbose,
                                 scan_chunk    = 1e5L,
                                 stable_rounds = 5L) {
  if (verbose) cat("  Scanning time column...\n")

  # Get column index from header
  header_line <- readLines(path, n = 1L)
  col_names   <- strsplit(header_line, sep, fixed = TRUE)[[1L]]
  col_idx     <- which(col_names == time_col)
  if (length(col_idx) == 0L)
    stop(sprintf("Column '%s' not found in file.", time_col))

  # Maintain file connection — no re-scanning from beginning
  con <- file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  readLines(con, n = 1L)  # skip header

  seen         <- integer(0)
  stable       <- 0L
  rows_scanned <- 0L

  repeat {
    lines <- readLines(con, n = as.integer(scan_chunk))
    if (length(lines) == 0L) break
    rows_scanned <- rows_scanned + length(lines)

    # Extract col_idx-th field from each line
    vals <- suppressWarnings(as.integer(
      vapply(strsplit(lines, sep, fixed = TRUE),
             function(x) if (length(x) >= col_idx) x[[col_idx]] else NA_character_,
             character(1L))
    ))
    vals <- vals[!is.na(vals)]

    new_vals <- setdiff(unique(vals), seen)
    if (length(new_vals) == 0L) {
      stable <- stable + 1L
      if (stable >= stable_rounds) break  # early exit
    } else {
      seen   <- sort(c(seen, new_vals))
      stable <- 0L
    }
  }

  if (verbose && rows_scanned < 1e7)
    cat(sprintf("  Early exit after scanning %s rows\n",
                format(rows_scanned, big.mark = ",")))

  sort(seen)
}

# --------------------------------------------------------------------------
# .find_unit_boundary(): index of the last row of the last complete unit
#
# Returns: integer index, or NA if the entire chunk is one unit
# --------------------------------------------------------------------------
.find_unit_boundary <- function(df, id_col) {
  ids     <- df[[id_col]]
  n       <- length(ids)
  if (n == 0L) return(NA_integer_)
  last_id <- ids[n]
  non_last <- which(ids != last_id)
  if (length(non_last) == 0L) return(NA_integer_)
  as.integer(max(non_last))
}

# =============================================================================
# Main function: otwfe_file()
# =============================================================================
#'
#' Two-Way Fixed Effects regression on a large CSV panel file
#'
#' @param path       Path to CSV file (must be sorted by id_col)
#' @param id_col     Name of the unit identifier column
#' @param time_col   Name of the calendar time column
#' @param y_col      Name of the dependent variable column
#' @param x_cols     Character vector of covariate column names
#' @param chunk_size Maximum number of rows per chunk (default 1,000,000)
#' @param sep        CSV delimiter (default ",")
#' @param verbose    Print progress messages
#' @return An object of class "otwfe"
#'
otwfe_file <- function(path,
                       id_col,
                       time_col,
                       y_col,
                       x_cols,
                       chunk_size = 1e6L,
                       sep        = ",",
                       verbose    = TRUE) {

  chunk_size <- as.integer(chunk_size)
  t_total    <- proc.time()
  sep72      <- strrep("=", 72)
  sep_m      <- strrep("-", 72)

  # -----------------------------------------------------------------------
  # Step 0: Input validation
  # -----------------------------------------------------------------------
  if (!file.exists(path))
    stop(sprintf("File not found: '%s'", path))

  col_names <- .csv_col_names(path, sep)
  for (col in c(id_col, time_col, y_col, x_cols))
    if (!col %in% col_names)
      stop(sprintf("Column not found in file: '%s'", col))

  if (verbose) {
    cat(sprintf("\n%s\n", sep72))
    cat(sprintf("  otwfe_file: %s\n", basename(path)))
    cat(sprintf("  id=%s  time=%s  y=%s  x=(%s)\n",
                id_col, time_col, y_col, paste(x_cols, collapse = ", ")))
    cat(sprintf("  chunk_size = %s rows\n", format(chunk_size, big.mark = ",")))
    cat(sprintf("%s\n", sep72))
  }

  # -----------------------------------------------------------------------
  # Step 1: T_support detection (early-termination time column scan)
  # -----------------------------------------------------------------------
  if (verbose) cat("\n[Step 1] Detecting T_support\n")
  t1           <- proc.time()
  time_levels  <- .detect_time_levels(path, time_col, sep, verbose)
  T_support    <- length(time_levels)
  time_remap   <- setNames(seq_along(time_levels), as.character(time_levels))
  baseline_idx <- 1L   # after remapping, time = 1 is the baseline

  if (verbose)
    cat(sprintf("  T_support = %d  |  time range: %s to %s  [%.1f sec]\n",
                T_support,
                time_levels[1L], time_levels[T_support],
                (proc.time() - t1)[["elapsed"]]))

  # -----------------------------------------------------------------------
  # Step 2: Build initialization chunk
  #   Accumulate chunks until all T calendar times are covered,
  #   ensuring warm-up can always span the full T_support.
  #
  #   Open a single sequential file connection here and reuse it through
  #   Step 4 — eliminates the O(skip) re-scan overhead of fread(skip=N).
  # -----------------------------------------------------------------------
  if (verbose) cat("\n[Step 2] Building initialization chunk (covering all T periods)\n")
  t2        <- proc.time()

  # Open once; on.exit ensures the connection is closed even on error
  con         <- file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  header_line <- readLines(con, n = 1L)   # consume header row

  init_df       <- NULL
  n_init_chunks <- 0L
  rows_read     <- 0L   # used only for verbose progress messages below

  repeat {
    raw <- .csv_read_chunk_seq(con, header_line, sep, chunk_size)
    if (is.null(raw) || nrow(raw) == 0L) break

    rows_read     <- rows_read + nrow(raw)
    n_init_chunks <- n_init_chunks + 1L

    raw[[time_col]] <- time_remap[as.character(raw[[time_col]])]
    init_df <- if (is.null(init_df)) raw else rbind(init_df, raw)

    times_found <- length(unique(init_df[[time_col]]))
    if (verbose)
      cat(sprintf("  Chunk %d accumulated: %s rows,  time coverage %d/%d\n",
                  n_init_chunks,
                  format(nrow(init_df), big.mark = ","),
                  times_found, T_support))

    if (times_found >= T_support) break
    if (nrow(raw) < chunk_size)    break   # EOF
  }
  rm(raw)

  if (is.null(init_df))
    stop("No data could be read from the file.")

  times_in_init <- length(unique(init_df[[time_col]]))
  if (times_in_init < T_support)
    warning(sprintf(
      "Initialization chunk covers only %d/%d time periods. Warm-up quality may be reduced.",
      times_in_init, T_support))

  # Detect unit boundary: incomplete last unit is carried over
  bnd <- .find_unit_boundary(init_df, id_col)
  if (is.na(bnd)) {
    first_chunk <- init_df
    carry_over  <- NULL
  } else {
    first_chunk <- init_df[seq_len(bnd),            , drop = FALSE]
    carry_over  <- init_df[(bnd + 1L):nrow(init_df), , drop = FALSE]
  }
  rm(init_df); invisible(gc())

  if (verbose)
    cat(sprintf("  -> First chunk: %s rows  |  carry-over: %s rows  [%.1f sec]\n",
                format(nrow(first_chunk), big.mark = ","),
                format(if (!is.null(carry_over)) nrow(carry_over) else 0L,
                       big.mark = ","),
                (proc.time() - t2)[["elapsed"]]))

  # -----------------------------------------------------------------------
  # Step 3: State initialization (warm-up selection + first chunk)
  # -----------------------------------------------------------------------
  if (verbose) cat("\n[Step 3] State initialization\n")
  handle <- otwfe_init(
    x_cols        = x_cols,
    time_col      = time_col,
    id_col        = id_col,
    y_col         = y_col,
    T_support     = T_support,
    baseline_time = baseline_idx,
    verbose       = verbose
  )
  handle <- otwfe_update(handle, first_chunk)
  rm(first_chunk); invisible(gc())

  # -----------------------------------------------------------------------
  # Step 4: Adaptive chunked processing with carry-over
  # -----------------------------------------------------------------------
  if (verbose) cat("\n[Step 4] Processing chunks\n")
  chunk_idx <- n_init_chunks

  repeat {
    raw <- .csv_read_chunk_seq(con, header_line, sep, chunk_size)
    if (is.null(raw) || nrow(raw) == 0L) break

    rows_read <- rows_read + nrow(raw)
    chunk_idx <- chunk_idx + 1L
    is_eof    <- nrow(raw) < chunk_size

    raw[[time_col]] <- time_remap[as.character(raw[[time_col]])]

    # Sort-order check: last id in carry_over <= first id in new chunk
    if (!is.null(carry_over) && nrow(raw) > 0L) {
      if (carry_over[[id_col]][nrow(carry_over)] > raw[[id_col]][1L])
        stop(sprintf(
          "Sort order violation at chunk %d: file must be sorted by '%s'.",
          chunk_idx, id_col))
    }

    combined <- if (is.null(carry_over)) raw else rbind(carry_over, raw)
    rm(raw); invisible(gc())

    if (is_eof) {
      carry_over <- NULL
      handle     <- otwfe_update(handle, combined)
      rm(combined); invisible(gc())
      break
    }

    bnd <- .find_unit_boundary(combined, id_col)

    if (is.na(bnd)) {
      carry_over <- combined
    } else {
      carry_over <- combined[(bnd + 1L):nrow(combined), , drop = FALSE]
      complete   <- combined[seq_len(bnd),              , drop = FALSE]
      rm(combined)
      handle     <- otwfe_update(handle, complete)
      rm(complete); invisible(gc())
    }

    if (verbose && chunk_idx %% 5L == 0L)
      cat(sprintf("  Chunk %d done | %s rows processed\n",
                  chunk_idx, format(rows_read, big.mark = ",")))
  }

  if (!is.null(carry_over) && nrow(carry_over) > 0L) {
    handle <- otwfe_update(handle, carry_over)
    rm(carry_over); invisible(gc())
  }

  # -----------------------------------------------------------------------
  # Step 5: Finalize
  # -----------------------------------------------------------------------
  if (verbose) cat("\n[Step 5] Finalizing\n")
  result <- otwfe_finalize(handle)

  t_elapsed <- (proc.time() - t_total)[["elapsed"]]
  if (verbose) {
    S <- result$state
    cat(sprintf("\n%s\n", sep_m))
    cat(sprintf("  Done: N = %s units,  n = %s obs,  p = %d\n",
                format(S$N, big.mark = ","),
                format(S$n, big.mark = ","),
                S$p))
    cat(sprintf("  Total elapsed: %.1f sec  (%.1f min)\n",
                t_elapsed, t_elapsed / 60))
    cat(sprintf("  State size: %.4f MB\n",
                as.numeric(object.size(S)) / 1e6))
    cat(sprintf("\n  Estimates (x covariates):\n"))
    cat(sprintf("    %-10s  %8s  %8s  %8s  %7s\n",
                "", "coef", "SE", "SE(CR)", "t"))
    cat(sprintf("    %s\n", strrep("-", 52)))
    theta <- S$theta_hat[x_cols]
    se_cl <- sqrt(diag(S$sigma2_hat * S$inv_dotZtZ)[x_cols])
    se_cr <- sqrt(diag(S$Vcr_hat)[x_cols])
    for (xn in x_cols)
      cat(sprintf("    %-10s  %8.4f  %8.4f  %8.4f  %7.2f\n",
                  xn, theta[xn], se_cl[xn], se_cr[xn],
                  theta[xn] / se_cr[xn]))
    cat(sprintf("    %s\n", strrep("-", 52)))
    cat("    * SE    : based on classical (homoskedastic) variance\n")
    cat("    * SE(CR): based on HC0 cluster-robust variance (Arellano 1987)\n")
    cat(sprintf("%s\n", sep_m))
  }

  result$time_levels_original <- time_levels
  result$time_remap            <- time_remap
  result
}
