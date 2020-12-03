#' Calculates empirical frequencies of each value in a collection of records.
#'
#' \code{emp_dist} takes a data.frame of records \code{data} and \code{key_vars}
#' a vector of colmun names, then returns a list of data.frames (one for each
#' key variable), each containing a column of values and a column of
#' corresponding empirical frequencies.
#'
#' @importFrom magrittr `%>%`
#'
#' @param data     A data.frame. The data containing the records.
#' @param key_vars A character vector. The names of the columns from \code{data}
#'                 to use as key variables.
#'
#' @return A list of data.frames. Each list contains an index of values each key
#'   variable can take on and the corresponding empirical frequency.
#'
emp_dist <- function(data, key_vars){
  emp_freqs <- list()
  for (key_var in key_vars) {
    emp_freq <- data %>%
      dplyr::select(key_var) %>%
      table %>%
      prop.table %>%
      data.frame
    names(emp_freq) <- c(paste0(key_var), "freq")
    emp_freqs[[key_var]] <- emp_freq
  }
  return(emp_freqs)
}
