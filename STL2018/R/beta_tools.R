#' Optimization function for backsolving shape parameters for a beta prior.
#'
#' @param par Parameters to optimize over.
#' @param mn  A float. The mean condition of the beta prior. Must be strictly
#'            positive.
#' @param var A float. The variance condition of the beta prior. Must be
#'            strictly positive.
#'
#' @return A vector of residuals.
#'
beta_optim <- function(par, mn, var){
  r <- rep(NA, length(par))
  r[1] <- -mn + abs(par[1]) / (abs(par[1]) + abs(par[2]))
  r[2] <- -var + abs(par[1])*abs(par[2]) / ((abs(par[1]) + abs(par[2]))^2 * (abs(par[1]) + abs(par[2]) + 1))
  return(r)
}

#' Solve for the shape parameters of beta distribtution given mean and variance.
#'
#' \code{beta_solve} backsolves for \eqn{\alpha} and \eqn{\beta} of a beta
#' distribution given a mean and variance.
#'
#' @param mn  A float. The mean condition of the beta prior. Must be strictly
#'            positive.
#' @param var A float. The variance condition of the beta prior. Must be
#'            strictly positive.
#'
#' @return A list. \eqn{(\alpha, \beta)}
#' @export
#'
#' @examples
#' beta_solve(mn = 2.6, var = 3.4)
#'
beta_solve <- function(mn, var){
  if (mn < 0) stop("Mean must be greater than zero")
  if (var < 0) stop("Variance must be greater than zero")
  soln <- BB::BBsolve(par = c(0.5, 0.5), fn = function(par) beta_optim(par, mn, var))
  return(list(alpha = abs(soln$par[1]), beta = abs(soln$par[2])))
}
