sampler_tests <- function(key_vars, burnin, sample_size, thin,
                          n_chains = 1, progressbar = TRUE,
                          alpha_init, lambda_init) {
  if (!is.logical(progressbar)) stop("progessbar must be either TRUE or FALSE.")
  if (!is.numeric(n_chains) | round(n_chains) != n_chains | n_chains <= 0) {
    stop("n_chains must be a positive integer.")
  }
  if (!is.numeric(burnin) | round(burnin) != burnin | burnin < 0) {
    stop("burnin must be a non-negative integer.")
  }
  if (!is.numeric(sample_size) | round(sample_size) != sample_size | sample_size <= 0) {
    stop("sample_size must be a positive integer.")
  }
  if (!is.numeric(thin) | round(thin) != thin | thin <= 0) {
    stop("thin must be a positive integer.")
  }
  if (!is.null(alpha_init)) {
    if (!is.numeric(alpha_init)) {
      stop("alpha_init must be a numeric vector.")
    }
    if (length(alpha_init) != length(key_vars)) {
      stop("alpha_init must be the same length as key_vars.")
    }
  }
  if (!is.null(lambda_init)) {
    if (!is.numeric(lambda_init) | !all(round(lambda_init) == lambda_init) | any(lambda_init <= 0)) {
      stop("lambda_init must be a vector of integer labels.")
    }
    if (length(lambda_init) != nrow(data)) {
      stop("The length of lambda_init must be equal to the number of rows in data.")
    }
  }
}

#' Title
#'
#' @param data
#' @param key_vars
#'
#' @importFrom magrittr `%>%`
#' @return
#'
#' @examples
data_tests <- function(data, key_vars) {
  if (!is.character(key_vars)) stop("Argument key_vars must be a character vector.")
  if ("RLid" %in% key_vars) stop("key_vars/data cannot have a variable named 'RLid'.")
  if ("linkage" %in% key_vars) stop("key_vars/data cannot have a variable named 'linkage'.")
  if (is.data.frame(data)) {
    if (!all(key_vars %in% names(data))) stop("A value in key_vars is not present as a column in data.")

    # converting to strings and handling NAs
    data %>%
      dplyr::mutate_at(key_vars, as.character) %>%
      dplyr::mutate_at(key_vars, ~ifelse(is.na(.), "", .)) %>%
      return

  } else {
    for (df in data) {
      if (!all(key_vars %in% names(df))) {
        stop("One or more of the data.frames in data does not include all of the columns in key_vars.")
      }
    }
    # combining all data.frames, converting to strings, and handling NAs
    data %>%
      lapply(function(df) dplyr::mutate_at(df, key_vars, as.character)) %>%
      Reduce(dplyr::bind_rows, .) %>%
      dplyr::mutate_at(key_vars, ~ifelse(is.na(.), "", .)) %>%
      return
  }
}
