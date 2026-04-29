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
# .csv_read_chunk(): read n_rows data rows starting after skip_rows data rows
#
# fread with skip=(skip_rows+1) + header=FALSE reads:
#   - skip_rows+1 total lines (1 header + skip_rows data rows)
#   - then reads the next n_rows data rows as col_names-named columns
#
# This is far faster than the previous readLines + paste + fread(string)
# approach: fread's C parser scans at ~50M rows/sec vs R's readLines ~1M/sec.
# --------------------------------------------------------------------------
.csv_read_chunk <- function(path, skip_rows, n_rows, sep, col_names) {
  # tryCatch: fread errors when skip >= total file lines (past-EOF condition)
  chunk <- tryCatch(
    data.table::fread(
      file         = path,
      skip         = skip_rows + 1L,   # +1 accounts for the header row
      nrows        = n_rows,
      sep          = sep,
      header       = FALSE,
      col.names    = col_names,
      showProgress = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(chunk) || nrow(chunk) == 0L) return(NULL)
  data.table::setDF(chunk)   # in-place class change, no copy
}

# --------------------------------------------------------------------------
# .detect_time_levels(): detect unique calendar times via early-termination scan
#
# Uses fread(select=time_col) so only the needed column is parsed.
# Early exit after stable_rounds consecutive scan-chunks with no new time values.
# --------------------------------------------------------------------------
.detect_time_levels <- function(path, time_col, sep, verbose,
                                 scan_chunk    = 1e5L,
                                 stable_rounds = 5L,
                                 col_names     = NULL) {
  if (verbose) cat("  Scanning time column...\n")

  if (is.null(col_names))
    col_names <- .csv_col_names(path, sep)

  # When header=FALSE, fread ignores col.names for select lookup — use index.
  time_col_idx <- which(col_names == time_col)
  if (length(time_col_idx) == 0L)
    stop(sprintf("Column '%s' not found in file.", time_col))

  seen         <- integer(0)
  stable       <- 0L
  rows_scanned <- 0L

  repeat {
    chunk <- tryCatch(
      data.table::fread(
        file         = path,
        skip         = rows_scanned + 1L,   # +1 accounts for header row
        nrows        = as.integer(scan_chunk),
        sep          = sep,
        header       = FALSE,
        select       = time_col_idx,        # integer index works with header=FALSE
        showProgress = FALSE
      ),
      error = function(e) NULL   # handles skip >= total file lines
    )
    if (is.null(chunk) || nrow(chunk) == 0L) break

    n_read       <- nrow(chunk)
    rows_scanned <- rows_scanned + n_read

    vals     <- suppressWarnings(as.integer(chunk[[1L]]))   # only column selected
    vals     <- vals[!is.na(vals)]
    new_vals <- setdiff(unique(vals), seen)

    if (length(new_vals) == 0L) {
      stable <- stable + 1L
      if (stable >= stable_rounds) break
    } else {
      seen   <- sort(c(seen, new_vals))
      stable <- 0L
    }

    if (n_read < as.integer(scan_chunk)) break  # EOF: no more rows
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
#' @param path       Path to CSV file
#' @param id_col     Name of the unit identifier column
#' @param time_col   Name of the calendar time column
#' @param y_col      Name of the dependent variable column
#' @param x_cols     Character vector of covariate column names
#' @param chunk_size Maximum number of rows per chunk (default 1,000,000)
#' @param sep        CSV delimiter (default ",")
#' @param auto_sort  If TRUE, sort the file by id_col before processing.
#'                   Required when the file is not already sorted by id_col
#'                   (e.g., year-sorted panel data). The sorted copy is written
#'                   to a temporary file and cleaned up automatically.
#'                   Requires loading the full file into memory; use FALSE for
#'                   files that are already sorted or too large to fit in RAM.
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
                       auto_sort  = FALSE,
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

  # -----------------------------------------------------------------------
  # auto_sort: read full file, sort by id_col, write to temp CSV, update path
  # -----------------------------------------------------------------------
  tmp_sorted <- NULL
  if (auto_sort) {
    if (verbose) cat("  [auto_sort] Sorting file by", id_col, "...\n")
    t_sort <- proc.time()
    dt     <- data.table::fread(path, sep = sep, showProgress = FALSE)
    data.table::setkeyv(dt, id_col)
    tmp_sorted <- tempfile(fileext = ".csv")
    data.table::fwrite(dt, tmp_sorted)
    rm(dt); invisible(gc())
    path <- tmp_sorted
    on.exit(unlink(tmp_sorted), add = TRUE)
    if (verbose)
      cat(sprintf("  [auto_sort] Done in %.1f sec\n",
                  (proc.time() - t_sort)[["elapsed"]]))
  }

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
  time_levels  <- .detect_time_levels(path, time_col, sep, verbose,
                                       col_names = col_names)
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
  #   Accumulate chunks until all T calendar times are covered.
  #   rows_file tracks how many data rows have been read from the file.
  # -----------------------------------------------------------------------
  if (verbose) cat("\n[Step 2] Building initialization chunk (covering all T periods)\n")
  t2            <- proc.time()

  init_df       <- NULL
  rows_file     <- 0L   # data rows consumed from file so far
  n_init_chunks <- 0L

  repeat {
    raw <- .csv_read_chunk(path, rows_file, chunk_size, sep, col_names)
    if (is.null(raw) || nrow(raw) == 0L) break

    rows_file     <- rows_file + nrow(raw)
    n_init_chunks <- n_init_chunks + 1L

    raw[[time_col]] <- time_remap[as.character(raw[[time_col]])]
    init_df <- if (is.null(init_df)) raw else rbind(init_df, raw)
    rm(raw)

    times_found <- length(unique(init_df[[time_col]]))
    if (verbose)
      cat(sprintf("  Chunk %d accumulated: %s rows,  time coverage %d/%d\n",
                  n_init_chunks,
                  format(nrow(init_df), big.mark = ","),
                  times_found, T_support))

    if (times_found >= T_support) break
    if (nrow(init_df) < n_init_chunks * chunk_size) break  # EOF
  }

  if (is.null(init_df))
    stop("No data could be read from the file.")

  times_in_init <- length(unique(init_df[[time_col]]))
  if (times_in_init < T_support)
    warning(sprintf(
      "Initialization chunk covers only %d/%d time periods. Warm-up quality may be reduced.",
      times_in_init, T_support))

  # Find unit boundary: the last complete unit in init_df
  bnd <- .find_unit_boundary(init_df, id_col)
  if (is.na(bnd)) {
    first_chunk <- init_df
    carry_over  <- NULL
  } else {
    first_chunk <- init_df[seq_len(bnd),             , drop = FALSE]
    carry_over  <- init_df[(bnd + 1L):nrow(init_df), , drop = FALSE]
  }
  rm(init_df); invisible(gc())

  n_carry <- if (!is.null(carry_over)) nrow(carry_over) else 0L
  if (verbose)
    cat(sprintf("  -> First chunk: %s rows  |  carry-over: %s rows  [%.1f sec]\n",
                format(nrow(first_chunk), big.mark = ","),
                format(n_carry, big.mark = ","),
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
  #   Each iteration: read chunk_size new rows from file, prepend carry_over,
  #   find the last complete unit boundary, process complete units, and
  #   carry the tail forward.
  #   rows_file advances by nrow(raw) each iteration (not by the boundary),
  #   so the next fread(skip=rows_file) picks up exactly at the file position
  #   after the current chunk — the carry_over rows are NOT re-read from disk.
  # -----------------------------------------------------------------------
  if (verbose) cat("\n[Step 4] Processing chunks\n")
  chunk_idx <- n_init_chunks

  repeat {
    raw <- .csv_read_chunk(path, rows_file, chunk_size, sep, col_names)
    if (is.null(raw) || nrow(raw) == 0L) {
      if (!is.null(carry_over) && nrow(carry_over) > 0L)
        handle <- otwfe_update(handle, carry_over)
      break
    }

    chunk_idx <- chunk_idx + 1L
    rows_file <- rows_file + nrow(raw)
    is_eof    <- nrow(raw) < chunk_size

    raw[[time_col]] <- time_remap[as.character(raw[[time_col]])]
    combined <- if (is.null(carry_over)) raw else rbind(carry_over, raw)
    rm(raw)

    if (is_eof) {
      handle <- otwfe_update(handle, combined)
      rm(combined); invisible(gc())
      break
    }

    bnd <- .find_unit_boundary(combined, id_col)

    if (is.na(bnd)) {
      carry_over <- combined
    } else {
      complete   <- combined[seq_len(bnd),              , drop = FALSE]
      carry_over <- combined[(bnd + 1L):nrow(combined), , drop = FALSE]
      rm(combined)
      handle     <- otwfe_update(handle, complete)
      rm(complete); invisible(gc())
    }

    if (verbose && chunk_idx %% 5L == 0L)
      cat(sprintf("  Chunk %d done | %s rows processed\n",
                  chunk_idx, format(rows_file, big.mark = ",")))
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
