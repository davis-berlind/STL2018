#' Tests for validity of \code{lambda_pior} argument.
#'
#' Both \code{\link{recursiveRL}} and \code{\link{regressionRL}} require an
#' \code{lambda_pior} argument that encodes the prior on the linkage
#' structure \eqn{\Lambda}. \code{lambda_prior_tests} runs tests to make sure
#' \code{lambda_pior} is a valid data structure.
#'
#' @section Details:
#' \code{lambda_prior} must include an element named "prior" that contains a
#' string specifying the prior to be used for the linkage structure. Currently
#' the values 'uniform' and 'PYP' are supported. If 'prior' = 'PYP', then
#' \code{lambda_prior} must also contain two numeric arguments 'nu' and
#' 'sigma'. See \code{\link{pyp_mean}} and \code{\link{pyp_var}} for more
#' detail.
#'
#' @param lambda_prior A list. A list containing the type of prior to use for
#'                     the linkage structure and the corresponding parameters
#'                     for the prior.
#'
#' @return Exception and error message if \code{lambda_prior} is invalid.
#' @export
#'
#' @examples
#'
#' uni_prior <- list(prior = "uniform")
#' pyp_params <- pyp_solve(mn = 450, var = 500, N = 500)
#' pyp_prior <- list(prior = "PYP",
#'                   nu = pyp_params$nu,
#'                   sigma = pyp_params$sigma)
#'
#' lapply(list(uni_prior, pyp_prior), lambda_prior_tests)
#'
lambda_prior_tests <- function(lambda_prior){
  implemented <- c("uniform", "PYP")
  if (!is.list(lambda_prior)) {
    stop("lambda_prior must be a list.")
  } else if (is.null(lambda_prior[["prior"]])) {
    stop("lambda_prior must contain an element named 'prior'.")
  } else if (!is.character(lambda_prior[["prior"]])) {
    stop("lambda_prior[['prior']] must must be a string specifying the prior.")
  } else if (!lambda_prior[["prior"]] %in% implemented) {
    msg <- paste(lambda_prior[["prior"]],
                 "is not currently implemented for lambda_prior. Please choose a value from",
                 paste(implemented, collapse = ", "),
                 "for lambda_prior[['prior']].")
  } else {
    if (lambda_prior[["prior"]] == "PYP") {
      if (is.null(lambda_prior[["nu"]])){
        stop("If lambda_prior[['prior']] = 'PYP', then lambda_prior must contain an element named 'nu'.")
      } else if (is.null(lambda_prior[["sigma"]])){
        stop("If lambda_prior[['prior']] = 'PYP', then lambda_prior must contain an element named 'sigma'.")
      }
      sigma <- lambda_prior[["sigma"]]
      nu <- lambda_prior[["nu"]]
      if (!is.numeric(sigma) | sigma >= 1) {
        stop("For PYP prior, lambda_prior[['sigma']] must be a number less than 1.")
      }
      if (!is.numeric(nu)) {
        stop("For PYP prior,lambda_prior[['nu']] must be a number.")
      }
      if ((sigma < 1 & sigma >= 0) & nu <= -sigma){
        stop("For PYP prior, if lambda_prior[['sigma']] is in [0,1), lambda_prior[['nu']] must be greater than negative sigma.")
      }
      if (sigma < 0 & nu <= abs(sigma)){
        stop("For PYP prior, if lambda_prior[['sigma']] is less than zero, lambda_prior[['nu']] must be greater than m*|sigma| for some integer m.")
      }
    }
  }
}

generate_lambda_prior <- function(){
  # TODO: implement
  next
}

#' Tests for validity of \code{alpha_pior} argument.
#'
#' Both \code{\link{recursiveRL}} and \code{\link{regressionRL}} require an
#' \code{alpha_pior} argument that encodes the beta prior on the distortion
#' probabilites. \code{alpha_prior_tests} runs tests to make sure
#' \code{alpha_pior} is a valid data structure.
#'
#' @section Details:
#' \code{alpha_pior} must be a list of lists. Each list in \code{alpha_pior}
#' must be named after an element of \code{key_vars}, and must contain two
#' strictly positive numeric elements \code{alpha} and \code{beta}.
#' \code{alpha_prior_tests} runs tests to make sure \code{alpha_prior} is
#' structured according to these standards.
#'
#' @family Prior Functions.
#'
#' @param key_vars    A character vector. The names of the fields in the record
#'                    data.
#' @param alpha_prior A list of lists. A list of beta distribution shape
#'                    parameters representing the prior for the distortion
#'                    probabilities.
#'
#' @return Exception and error message if \code{alpha_prior} is invalid.
#'
#' @export
#'
#' @examples
#' # creating uniform prior for RL500
#'
#' key_vars <- c("fname_c1", "lname_c1", "by", "bm", "bd")
#'
#' alpha_prior <- list(fname_c1 = list(alpha = 1, beta = 1),
#'                     lname_c1 = list(alpha = 1, beta = 1),
#'                     by       = list(alpha = 1, beta = 1),
#'                     bm       = list(alpha = 1, beta = 1),
#'                     bd       = list(alpha = 1, beta = 1))
#'
#' alpha_prior_tests(key_vars, alpha_prior)
#'
alpha_prior_tests <- function(alpha_prior, key_vars){
  if (!is.list(alpha_prior)){
    stop("alpha_prior must be a list.")
  } else {
    for (key_var in key_vars){
      if (!key_var %in% names(alpha_prior)){
        stop(paste("There is no element in alpha_prior named", key_var))
      }
      else if (!is.list(alpha_prior[[key_var]])){
        stop(paste(key_var, "in alpha_prior is not a list object."))
      }
      else if (!is.list(alpha_prior[[key_var]])){
        stop(paste(key_var, "in alpha_prior is not a list object."))
      }
      else if (!all(c("alpha", "beta") %in% names(alpha_prior[[key_var]]))){
        msg <- paste(key_var, "in alpha_prior must contain two numeric",
                     "elements called alpha and beta.")
        stop(msg)
      } else {
        for (arg in c("alpha", "beta")){
          val <- alpha_prior[[key_var]][[arg]]
          if (!is.numeric(val) | val <= 0){
            stop(paste("For",key_var,arg,"must be a strictly postive number."))
          }
        }
      }
    }
  }
}

#' Generate prior for distortion probability \eqn{\alpha}.
#'
#' Given a vector of key variable names, \code{generate_alpha_prior} prompts
#' the user to enter the parameters of the prior on the distortion
#' probability. The list that \code{generate_alpha_prior} can then directly
#' be used as the \code{alpha_prior} argument in either
#' \code{\link{recursiveRL}} or \code{\link{regressionRL}}.
#'
#' @family Prior Functions.
#'
#' @param key_vars A character vector. The names of the fields in the record
#'                 data.
#'
#' @return A formatted list encoding the prior for \eqn{\alpha}.
#' @export
#'
#' @examples
#' key_vars <- c("fname_c1", "lname_c1", "by", "bm", "bd")
#' alpha_prior <- generate_alpha_prior(key_vars)
#'
generate_alpha_prior <- function(key_vars){
  alpha_prior <- list()
  while (TRUE) {
    m <- readline("Use uniform prior for every field in key_vars? [y/n] ")
    if (m %in% c("y", "n")) break
  }
  for (key_var in key_vars){
    if (m == "y"){
      alpha_prior[[key_var]] <- list(alpha = 1, beta = 1)
    } else {
      while (TRUE) {
        alpha <- readline(paste("For", key_var, "enter a number greater than zero for alpha: "))
        if (!is.na(suppressWarnings(as.numeric(alpha))) & alpha > 0) break
      }
      while (TRUE) {
        beta <- readline(paste("For", key_var, "enter a number greater than zero for beta: "))
        if (!is.na(suppressWarnings(as.numeric(beta))) & beta > 0) break
      }
      alpha_prior[[key_var]] <- list(alpha = as.numeric(alpha), beta = as.numeric(beta))
    }
  }
  return(alpha_prior)
}

#' Log density of \code{beta} given a prior.
#'
#' @param beta       A numeric vector. The sample of \eqn{\beta} to
#'                   calculate the log density for
#' @param beta_prior A list. A list encoding the prior for \eqn{\beta}.
#'
#' @return Log density of \code{beta}.
#'
log_beta_prior_dens <- function(beta, beta_prior){
  if (beta_prior[["prior"]] == "flat") return(0)
}

#' Log density of \code{sigma_y} given a prior.
#'
#' @param sigma_y       A numeric vector. The sample of \eqn{\sigma^_y} to
#'                      calculate the log density for.
#' @param sigma_y_prior A list. A list encoding the prior for \eqn{\sigma^_y}.
#'
#' @return Log density of \code{sigma_y}.
#'
log_sigma_y_prior_dens <- function(sigma_y, sigma_y_prior){
  if (sigma_y_prior[["prior"]] == "flat") return(0)
}

#' Log density of \code{sigma_x} given a prior.
#'
#' @param sigma_x       A numeric vector. The sample of
#'                      \eqn{\Sigma^_{x|\widetilde{x}} to calculate the log
#'                      density for.
#' @param sigma_x_prior A list. A list encoding the prior for
#'                      \eqn{\Sigma^_{x|\widetilde{x}}.
#'
#' @return Log density of \code{sigma_x}.
#'
log_sigma_x_prior_dens <- function(sigma_x, sigma_x_prior){
  if (sigma_x_prior[["prior"]] == "flat") return(0)
}

#' Log density of \code{sigma_data} given a prior.
#'
#' @param sigma_data       A numeric vector. The sample of
#'                         \eqn{\Sigma^_{\widetilde{x}} to calculate the log
#'                         density for.
#' @param sigma_data_prior A list. A list encoding the prior for
#'                         \eqn{\Sigma^_{\widetilde{x}}.
#'
#' @return Log density of \code{sigma_data}.
#'
log_sigma_data_prior_dens <- function(sigma_data, sigma_data_prior){
  if (sigma_data_prior[["prior"]] == "flat") return(0)
}
