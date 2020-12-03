#' List the implemented proposal distributions for the distortion probability
#' \eqn{\alpha}.
#'
#' Returns a list of the proposal distributions that are currently
#' implemented for generating new values of  the distortion probability
#' \eqn{\alpha} in \code{\link{alpha_metropolis}}.
#'
#' @param verbose A logical. If TRUE, the list of implemented proposal
#'                distributions will be printed
#'
#' @return A character vector of the implemented proposal distributions.
#' @export
#'
get_alpha_proposals <- function(verbose = TRUE){
  implemented <- c("RRW")
  description <- c("reflected random walk")
  if (verbose) cat(paste0("\n",implemented, ": ", description))
  else return(implemented)
}

#' Tests for validity of \code{alpha_proposal} argument.
#'
#' \code{alpha_proposal_tests} performs tests to make sure \code{alpha_proposal}
#' is properly formatted for use in either \code{\link{recursiveRL}} or
#' \code{\link{regressionRL}}.
#'
#' @section Details:
#' \code{alpha_proposal} must be either a list containing a string element
#' named \code{proposal}, which specifies the type of proposal distribution to
#' use for the distortion probabilities \eqn{\alpha}, and a non-negative numeric
#' element named \code{drift}, which controls the correlation of the proposals.
#' Otherwise, \code{alpha_proposal} can be a list of lists, with each sublist being
#' named after a value in \code{key_vars} and each containing separate \code{proposal}
#' and \code{drift} parameters. Combinations are also allowed, so that if some value
#' in \code{key_vars} does not appear in \code{alpha_proposal},
#' \code{alpha_proposal[["proposal"]]} and \code{alpha_proposal[["drift"]]} will be used
#' for the missing values.
#'
#' @seealso \code{\link{get_alpha_proposals}} for currently implemented values of
#' \code{proposal}.
#'
#' @importFrom dplyr setequal
#'
#' @param alpha_proposal A list. A list containing the proposal distribution to
#'                       use in the Metropolis-Hastings step of
#'                       \code{\link{recursiveRL}} and \code{\link{regressionRL}}.
#' @param key_vars       A character vector. The names of the fields in the record
#'                       data.
#'
#' @return Exception and error message if \code{alpha_proposal} is invalid.
#'
#' @export
#'
#' @examples
#' # creating reflected random walk proposal for RL500
#' key_vars <- c("fname_c1", "lname_c1", "by", "bm", "bd")
#'
#' # one proposal for all fields
#' alpha_proposal <- list(proposal = "RRW", drift = 0.1)
#' alpha_proposal_tests(alpha_proposal, key_vars)
#'
#' # different proposal for each field
#' alpha_proposal <- list(fname_c1 = list(proposal = "RRW", drift = 0.1),
#'                        lname_c1 = list(proposal = "RRW", drift = 0.2),
#'                        by       = list(proposal = "RRW", drift = 0.3),
#'                        bm       = list(proposal = "RRW", drift = 0.4),
#'                        bd       = list(proposal = "RRW", drift = 0.5))
#' alpha_proposal_tests(alpha_proposal, key_vars)
#'
alpha_proposal_tests <- function(alpha_proposal, key_vars){
  params <- c("proposal", "drift")
  lnames <- names(alpha_proposal)
  if (!is.list(alpha_proposal)) {
    stop("alpha_proposal must be a list specifying the proposal distribution to use.")
  }
  if (!setequal(key_vars, lnames) & !all(params %in% lnames)) {
    stop(paste("alpha_proposal must either contain two elements 'proposal' and 'drift'",
               "or contain lists with names corresponding to every element of key_vars"))
  }
  if (all(params %in% lnames)) {
    proposal <- alpha_proposal[["proposal"]]
    drift <- alpha_proposal[["drift"]]
    if (!is.numeric(drift) | drift < 0) {
      stop("Drift parameters must be non-negative numbers.")
    }
    if (!is.character(proposal) | !proposal %in% get_alpha_proposals(verbose = FALSE)) {
      msg <- paste(proposal, "is not an implemented proposal distribution for alpha_metropolis.",
                   "See get_alpha_proposals() for a list of available proposal distributions.")
      stop(msg)
    }
  }
  if (any(key_vars %in% lnames)) {
    for (key_var in intersect(key_vars, lnames)) {
      proposal <- alpha_proposal[[key_var]][["proposal"]]
      drift <- alpha_proposal[[key_var]][["drift"]]
      if (!is.numeric(drift) | drift < 0) {
        stop("Drift parameters must be non-negative numbers.")
      }
      if (!is.character(proposal) | !proposal %in% get_alpha_proposals(verbose = FALSE)) {
        msg <- paste(proposal, "is not an implemented proposal distribution for alpha_metropolis.",
                     "See get_alpha_proposals() for a list of available proposal distributions.")
        stop(msg)
      }
    }
  }
}

# while (TRUE) {
#   proposal <- readline(paste0("Specify a proposal distribution for ", key_var,"? [y/n] "))
#   if (proposal == "y"){
#     while (TRUE) {
#       proposal <- readline(paste("Choose a proposal distribution from the list above:", get_alpha_proposals()))
#       if (proposal %in% get_alpha_proposals(verbose = FALSE)) break
#     }
#     alpha_prior[[key_var]][["proposal"]] <- proposal
#     break
#   }
#   else if (proposal == "n") break
# }
# while (TRUE) {
#   drift <- readline(paste0("Specify a drift parameter for ", key_var,"? [y/n] "))
#   if (drift == "y"){
#     while (TRUE) {
#       if (!is.null(alpha_prior[[key_var]][["proposal"]])){
#         if (alpha_prior[[key_var]][["proposal"]] == "RRW"){
#           drift <- readline("Enter a number in [0,1]: ")
#           if (!is.na(suppressWarnings(as.numeric(drift))) & drift >= 0 & drift <= 1) break
#         }
#       } else{
#         drift <- readline("Enter a non-negative number: ")
#         if (!is.na(suppressWarnings(as.numeric(drift))) & drift >= 0) break
#       }
#     }
#     alpha_prior[[key_var]][["drift"]] <- as.numeric(drift)
#     break
#   }
#   else if (drift == "n") break
# }
