#' Utility functions for the standartox package

#' Read binary vector function
#'
#' Reads a binary vector and returns the deserialized data object. Supports both
#' RDS (compressed R objects) and FST (Fast Serialization Table) formats.
#'
#' @param vec raw vector; A binary raw vector containing serialized data.
#' @param type character; The format type of the binary data. Must be either
#'   \code{"rds"} (default) for R serialized objects or \code{"fst"} for
#'   Fast Serialization Table format.
#'
#' @return For \code{type = "rds"}, returns the deserialized R object.
#'   For \code{type = "fst"}, returns a \code{data.table} object.
#'
#' @details
#' This function is used internally to deserialize binary data downloaded from
#' web sources. The RDS format uses gzip compression and is suitable for any R
#' object. The FST format is optimized for data frames and tables.
#'
#' The connection is automatically closed using \code{on.exit()} for RDS files.
#' For FST files, a temporary file is created and automatically cleaned up.
#'
#' @author Andreas Scharmueller \email{andschar@@protonmail.com}
#' @noRd
#'
read_bin_vec = function(vec, type = c('rds', 'fst')) {
  # Input validation
  if (!is.raw(vec)) {
    stop("'vec' must be a raw vector")
  }
  type = match.arg(type)

  if (type == 'rds') {
    # from https://stackoverflow.com/questions/58135794/read-binary-vector/58136567#58136567
    con = gzcon(rawConnection(vec))
    res = readRDS(con)
    on.exit(close(con))
  }
  if (type == 'fst') {
    tmp = tempfile()
    on.exit(unlink(tmp), add = TRUE)
    writeBin(vec, tmp)
    res = fst::read_fst(tmp, as.data.table = TRUE)
  }
  res
}

#' Convert CAS number between hyphenated and non-hyphenated formats
#'
#' Converts Chemical Abstracts Service (CAS) Registry Numbers between hyphenated
#' (e.g., "50-00-0") and non-hyphenated (e.g., "50000") formats. The function
#' automatically detects the input format and converts to the opposite format.
#'
#' @param cas character; A vector of CAS numbers in either hyphenated or
#'   non-hyphenated format. The function automatically detects which format
#'   is provided.
#'
#' @return A character vector of CAS numbers in the opposite format from the input.
#'   If input contains hyphens, returns non-hyphenated format. If input lacks
#'   hyphens, returns hyphenated format (XXX-XX-X).
#'
#' @details
#' CAS numbers are unique identifiers for chemical substances. The standard format
#' includes hyphens (e.g., "50-00-0" for formaldehyde), but some databases store
#' them without hyphens. This function facilitates conversion between these formats.
#'
#' The function assumes valid CAS number structure where the last three characters
#' represent the check digit (1 char) and the two preceding digits.
#'
#' @author Andreas Scharmueller \email{andschar@@protonmail.com}
#' @noRd
#'
cas_conv = function(cas) {
  # Input validation
  if (!is.character(cas)) {
    stop("'cas' must be a character vector")
  }
  if (length(cas) == 0) {
    return(character(0))
  }

  if (any(grepl('-', cas))) {
    gsub('-', '', cas)
  } else {
    # Ensure CAS number has at least 4 characters for proper formatting
    if (any(nchar(cas) < 4)) {
      warning("Some CAS numbers are too short for proper formatting")
    }
    paste(substr(cas, 1, nchar(cas)-3),
          substr(cas, nchar(cas)-2, nchar(cas)-1),
          substr(cas, nchar(cas), nchar(cas)),
          sep = '-')
  }
}

#' Compose and display query parameter message
#'
#' Creates a formatted message displaying query parameters for Standartox database
#' queries. The function filters out NULL parameters, truncates long values,
#' and formats them in a user-friendly way.
#'
#' @param body list; A named list where each element represents a query parameter
#'   and its value(s). NULL values are automatically removed before display.
#'
#' @return NULL (invisible). The function is called for its side effect of
#'   printing a formatted message to the console via \code{message()}.
#'
#' @details
#' This function is used internally to provide user feedback about the parameters
#' being used in a Standartox query. Key features:
#' \itemize{
#'   \item Filters out NULL parameters
#'   \item Collapses multiple values with commas
#'   \item Truncates long parameter strings (>= 80 characters) to 70 characters
#'   \item Formats output with parameter names and values
#' }
#'
#' @author Andreas Scharmueller \email{andschar@@protonmail.com}
#' @noRd
#'
stx_message = function(body) {
  # Input validation
  if (!is.list(body)) {
    stop("'body' must be a list")
  }

  body2 = body[ sapply(body, function(x) !is.null(x)) ]
  body2 = lapply(body2, paste0, collapse = ', ')
  body2 = lapply(body2, function(x) {
    if (nchar(x) >= 80) {
      paste0(substr(paste0(x, collapse = ''), 1, 70), '...[truncated]')
    } else {
      x
    }
  })
  msg = paste0(paste0(names(body2), ': ', unlist(body2)),
               collapse = '\n')
  message('Standartox query running..\nParameters:\n', msg)
  invisible(NULL)
}

#' Remove columns that contain only NA values
#'
#' Filters out columns from a data frame or list where all values are NA.
#' This is useful for cleaning up query results or data tables that may
#' contain completely missing columns.
#'
#' @param dt data.frame, data.table, or list; The data structure to filter.
#'   Works with any object that can be passed to \code{Filter()}.
#'
#' @return An object of the same type as the input, but with all-NA columns
#'   removed. If the input is a data frame or data.table, returns a data frame
#'   or data.table with fewer columns. If all columns are NA, returns an
#'   empty structure of the same type.
#'
#' @details
#' This function uses R's \code{Filter()} function to remove columns where
#' all values are NA. It's particularly useful after subsetting data or
#' performing joins where certain columns may become entirely missing.
#'
#' @examples
#' \dontrun{
#' # Create a data frame with some all-NA columns
#' df <- data.frame(
#'   a = c(1, 2, 3),
#'   b = c(NA, NA, NA),
#'   c = c("x", "y", "z"),
#'   d = c(NA, NA, NA)
#' )
#' # Remove all-NA columns
#' rm_col_na(df)
#' # Returns only columns 'a' and 'c'
#' }
#'
#' @author Andreas Scharmueller \email{andschar@@protonmail.com}
#' @noRd
#'
rm_col_na = function(dt) {
  Filter(function(x) !all(is.na(x)), dt)
}

#' Calculate geometric mean
#'
#' Computes the geometric mean of a numeric vector. The geometric mean is
#' particularly useful for data that spans several orders of magnitude or
#' for concentrations in ecotoxicological studies.
#'
#' @param x numeric vector; A vector of positive values for which to compute
#'   the geometric mean. Negative values will result in NaN.
#' @param na.rm logical; If \code{TRUE} (default), NA values are removed
#'   before calculation. If \code{FALSE}, NA values will propagate.
#' @param zero.propagate logical; If \code{FALSE} (default), zero values are
#'   excluded from the calculation. If \code{TRUE}, the presence of any zero
#'   will cause the function to return 0.
#'
#' @return A numeric value representing the geometric mean. Returns \code{NaN}
#'   if any values are negative. Returns 0 if \code{zero.propagate = TRUE} and
#'   any zero values are present.
#'
#' @details
#' The geometric mean is calculated as \code{exp(mean(log(x)))} for positive values.
#' This implementation handles zeros and negative values explicitly:
#' \itemize{
#'   \item Negative values: Returns NaN with a warning
#'   \item Zero values: Either returns 0 (if \code{zero.propagate = TRUE}) or
#'         excludes them from calculation (default)
#' }
#'
#' When \code{zero.propagate = FALSE}, the denominator uses the original length
#' of the vector, not just the count of positive values.
#'
#' @references
#' \url{https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in}
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' gm_mean(c(2, 8))  # Returns 4
#'
#' # With zeros (default behavior - zeros excluded)
#' gm_mean(c(1, 2, 0, 4))
#'
#' # With zeros (propagate zeros)
#' gm_mean(c(1, 2, 0, 4), zero.propagate = TRUE)  # Returns 0
#'
#' # With NA values
#' gm_mean(c(2, 8, NA), na.rm = TRUE)
#' }
#'
#' @author Andreas Scharmueller \email{andschar@@protonmail.com}
#' @noRd
#'
# function to calculate the geometric mean
# https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
gm_mean = function(x, na.rm = TRUE, zero.propagate = FALSE){
  # Input validation
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector")
  }

  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }
}

#' Calculate geometric standard deviation
#'
#' Computes the geometric standard deviation of a numeric vector. The geometric
#' standard deviation is the multiplicative analog of the standard deviation and
#' is useful for log-normally distributed data such as toxicological concentrations.
#'
#' @param x numeric vector; A vector of strictly positive values. Non-positive
#'   values will result in NA with a warning.
#' @param na.rm logical; If \code{TRUE} (default), NA values are removed before
#'   calculation. If \code{FALSE} and NAs are present, the function returns NA.
#' @param sqrt.unbiased logical; If \code{TRUE} (default), uses the unbiased
#'   estimator with n-1 degrees of freedom. If \code{FALSE}, uses the maximum
#'   likelihood estimator with n degrees of freedom.
#'
#' @return A numeric value representing the geometric standard deviation.
#'   Returns \code{NA} with a warning if any values are non-positive or if
#'   \code{na.rm = FALSE} and NA values are present.
#'
#' @details
#' The geometric standard deviation is calculated as \code{exp(sd(log(x)))}.
#' This implementation:
#' \itemize{
#'   \item Requires all values to be strictly positive (> 0)
#'   \item Handles NA values according to the \code{na.rm} parameter
#'   \item Provides both unbiased (n-1) and maximum likelihood (n) estimators
#'   \item Returns NA with a warning for non-positive values
#' }
#'
#' The geometric standard deviation is always >= 1. Values close to 1 indicate
#' low relative variability, while larger values indicate higher variability.
#'
#' @references
#' Based on \code{EnvStats::geoSD()}
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' gm_sd(c(2, 4, 8))
#'
#' # With NA values
#' gm_sd(c(2, 4, NA, 8), na.rm = TRUE)
#'
#' # Using maximum likelihood estimator
#' gm_sd(c(2, 4, 8), sqrt.unbiased = FALSE)
#'
#' # Non-positive values produce NA with warning
#' gm_sd(c(2, 0, 4))  # Returns NA
#' }
#'
#' @author Andreas Scharmüller \email{andschar@@protonmail.com}
#' @noRd
#'
gm_sd = function (x, na.rm = TRUE, sqrt.unbiased = TRUE) {
  # after EnvStats::geoSD()
  if (!is.vector(x, mode = "numeric") || is.factor(x))
    stop("'x' must be a numeric vector")
  wna = which(is.na(x))
  if (length(wna) != 0) {
    if (na.rm)
      x = x[-wna]
    else return(NA)
  }
  if (any(x <= 0)) {
    warning("Non-positive values in 'x'")
    return(NA_real_)
  } else {
    sd.log = stats::sd(log(x))
    if (!sqrt.unbiased) {
      n = length(x)
      sd.log = sqrt((n - 1)/n) * sd.log
    }
    exp(sd.log)
  }
}

#' Flag outliers using the interquartile range (IQR) method
#'
#' Identifies outliers in a numeric vector using the Tukey fence method based on
#' the interquartile range (IQR). Values beyond a multiple of the IQR from the
#' quartiles are flagged as outliers.
#'
#' @param x numeric vector; The data to check for outliers.
#' @param lim numeric; The multiplier for the IQR to define outlier boundaries.
#'   Default is 1.5 (standard Tukey method). Common alternatives are 3 for
#'   "far outliers" or other values for custom sensitivity.
#' @param na.rm logical; If \code{TRUE} (default), NA values are ignored when
#'   calculating quartiles and IQR. If \code{FALSE}, presence of NA values will
#'   cause the function to fail.
#' @param ... Additional arguments passed to \code{stats::quantile()} and
#'   \code{stats::IQR()}.
#'
#' @return A logical vector of the same length as \code{x}, where \code{TRUE}
#'   indicates an outlier and \code{FALSE} indicates a normal value. NA values
#'   in the input are preserved as NA in the output (when \code{na.rm = TRUE}).
#'
#' @details
#' The function implements the standard Tukey fence method for outlier detection:
#' \itemize{
#'   \item Lower fence: Q1 - (lim × IQR)
#'   \item Upper fence: Q3 + (lim × IQR)
#'   \item Values outside these fences are flagged as outliers
#' }
#'
#' The default \code{lim = 1.5} corresponds to the standard boxplot definition
#' of outliers. Using \code{lim = 3} identifies only extreme outliers.
#'
#' This function uses \code{data.table::fifelse()} for efficient vectorized
#' conditional evaluation, which properly handles NA values.
#'
#' @examples
#' \dontrun{
#' # Basic usage with default IQR multiplier (1.5)
#' x <- c(1, 2, 3, 4, 5, 100)
#' flag_outliers(x)  # Returns: FALSE FALSE FALSE FALSE FALSE TRUE
#'
#' # With NA values
#' x_na <- c(1, 2, 3, NA, 5, 100)
#' flag_outliers(x_na, na.rm = TRUE)
#'
#' # Using stricter criterion for far outliers
#' flag_outliers(x, lim = 3)
#'
#' # In a data.table context
#' library(data.table)
#' dt <- data.table(concentration = c(0.1, 0.2, 0.15, 10, 0.18))
#' dt[, is_outlier := flag_outliers(concentration)]
#' }
#'
#' @author Andreas Scharmueller \email{andschar@@protonmail.com}
#' @noRd
#'
flag_outliers = function(x, lim = 1.5, na.rm = TRUE, ...) {
  # Input validation
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector")
  }
  if (!is.numeric(lim) || length(lim) != 1 || lim <= 0) {
    stop("'lim' must be a single positive number")
  }

  qnt = stats::quantile(x, probs = c(.25, .75), na.rm = na.rm, ...)
  H = lim * stats::IQR(x, na.rm = na.rm, ...)
  data.table::fifelse(x < qnt[1] - H | x > qnt[2] + H, TRUE, FALSE)
}
