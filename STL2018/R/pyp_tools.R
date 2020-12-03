#' Calculate \eqn{\log \frac{\Gamma(a + b)}{\Gamma(a)}}.
#'
#' @family PYP Functions
#'
#' @param a,b A float. \eqn{a > 0} and \eqn{b > -a}.
#'
#' @return \eqn{\log \frac{\Gamma(a + b)}{\Gamma(a)}}
#'
lgratio <- function(a,b){
  lgamma(a+b) - lgamma(a)
}

#' PYP function parameter tests.
#'
#' @family PYP Functions
#'
#' @param nu,sigma A float. PYP parameters such that \eqn{\sigma \in [0,1)} and
#'                 \eqn{\nu > -\sigma} or \eqn{\sigma < 0} and \eqn{\nu > m|\sigma|}
#'                 for some integer \eqn{m}.
#' @param N        An integer. The number of records.
#'
#' @return Warning if broken parameters are provided.
#'
pyp_warnings <- function(nu, sigma, N){
  if (!is.numeric(N) | round(N) != N | N <= 0) stop("N must be a positive integer.")
  if (!is.numeric(sigma) | sigma >= 1) stop("sigma must be a number less than 1.")
  if (!is.numeric(nu)) stop("nu must be a number.")
  if ((sigma < 1 & sigma >= 0) & nu <= -sigma){
    stop("If sigma is in [0,1), nu must be greater than negative sigma.")
  }
  if (sigma < 0 & nu <= abs(sigma)){
    stop("If sigma is less than zero, nu must be greater than m*|sigma| for some integer m.")
  }
}

#' Calculate the mean of a Pittman-Yor process.
#'
#' Given \eqn{N} entities (in this case records), \code{pyp_mean} calculates the
#' expected number of clusters for a Pittman-Yor process with parameters
#' \eqn{\sigma} and \eqn{\nu}, \deqn{\frac{\nu}{\sigma} \left[\frac{(\nu +
#' \sigma)_{N\uparrow}}{\nu_{N\uparrow}} - 1\right]} where \eqn{x_{s\uparrow} =
#' \frac{\Gamma(x+s)}{\Gamma(x)}}.
#'
#' @family PYP Functions
#'
#' @param nu,sigma A float. PYP parameters such that \eqn{\sigma \in [0,1)} and
#'                 \eqn{\nu > -\sigma} or \eqn{\sigma < 0} and \eqn{\nu > m|\sigma|}
#'                 for some integer \eqn{m}.
#' @param N        An integer. The number of records.
#'
#' @return Mean of PYP prior.
#' @export
#'
pyp_mean <- function(nu, sigma, N){
  pyp_warnings(nu, sigma, N)
  mn <- (nu / sigma)*(exp(lgratio(nu + sigma, N) - lgratio(nu, N)) - 1)
  return(mn)
}

#' Calculate the variance of a Pittman-Yor process.
#'
#' Given \eqn{N} entities (in this case records), \code{pyp_var} calculates the
#' variance of the number of clusters for a Pittman-Yor process with parameters
#' \eqn{\sigma} and \eqn{\nu},
#' \deqn{\frac{\nu(\nu+\sigma)}{\sigma^2}\frac{(\nu+2\sigma)_{N\uparrow}}{\nu_N\uparrow}
#' - \left[\frac{\nu}{\sigma}\frac{(\nu +
#' \sigma)_{N\uparrow}}{\nu_{N\uparrow}}\right]^2 - \frac{\nu}{\sigma}\frac{(\nu
#' + \sigma)_{N\uparrow}}{\nu_{N\uparrow}}} where \eqn{x_{s\uparrow} =
#' \frac{\Gamma(x+s)}{\Gamma(x)}}.
#'
#' @family PYP Functions
#'
#' @param nu,sigma A float. The PYP parameters such that \eqn{\sigma \in [0,1)}
#'                 and \eqn{\nu > - \sigma} or \eqn{\sigma < 0} and
#'                 \eqn{\nu > m|\sigma|} for some integer \eqn{m}.
#' @param N        An integer. The number of records.
#'
#' @return Variance of PYP prior.
#' @export
#'
pyp_var <- function(nu, sigma, N){
  pyp_warnings(nu, sigma, N)
  p1 <- exp(lgratio(nu + 2*sigma, N) - lgratio(nu, N))
  p2 <- exp(lgratio(nu + sigma, N) - lgratio(nu, N))
  ret <- (nu*(nu + sigma) / sigma^2)*p1 - ((nu /sigma)*p2)^2 - (nu / sigma)*p2
  return(ret)
}

#' Optimization function for backsolving parameters of a PYP prior.
#'
#' Gives non-linear equations for \eqn{\nu} and \eqn{\sigma} in terms of a known
#' mean, variance, and \eqn{N}.
#'
#' @family PYP Functions
#'
#' @param par Parameters to optimize over.
#' @param mn  An integer. The mean condition of the PYP prior.
#' @param var A float. The variance condition of the PYP prior. Must be strictly
#'            positive.
#' @param N   An integer. The number of records.
#'
#' @return A vector of residuals.
#'
pyp_optim <- function(par, mn, var, N){
  r <- rep(NA, length(par))
  sigma <- min(par[1], 0.999)
  nu <- max(par[2], -sigma + 0.001)

  r[1] <- -mn + pyp_mean(nu, sigma, N)
  r[2] <- -var + pyp_var(nu, sigma, N)

  return(r)
}

#' Solve for parameters a PYP prior given a mean, variance, and \eqn{N}.
#'
#' \code{pyp_solve} finds the parameters \eqn{\nu} and \eqn{\sigma} that give a
#' PYP prior given the mean \code{mn}, variance \code{var}, for a given number
#' of entities \eqn{N}.
#'
#' @param mn  An integer. The mean condition of the PYP prior.
#' @param var A float. The variance condition of the PYP prior. Must be
#'            strictly positive.
#' @param N   An integer. The number of records.
#'
#' @family PYP Functions
#'
#' @return A list. \eqn{(\sigma, \nu).}
#' @export
#'
#' @examples
#' pyp_solve(mn = 450, var = 500, N = 500)
pyp_solve <- function(mn, var, N) {
  if (!is.numeric(mn) | mn <= 0) stop("mn must be a positive number.")
  if (!is.numeric(var) | var <= 0) stop("var must be a positive number.")
  if (!is.numeric(N) | round(N) != N | N <= 0) stop("N must be a positive integer.")
  soln <- BB::BBsolve(par = c(0.5, 0.5), fn = function(par) pyp_optim(par, mn, var, N))
  sigma <- min(0.999, soln$par[1])
  if (soln$par[1] < 0){
    nu <- round(soln$par[2] / abs(sigma)) * abs(sigma)
    return(list(sigma = sigma, nu = nu))
  }
  else{
    return(list(sigma = sigma, nu = soln$par[2]))
  }
}
