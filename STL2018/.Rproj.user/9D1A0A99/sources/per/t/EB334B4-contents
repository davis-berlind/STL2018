#' Metropolis-Hastings update of distortion probability for \eqn{\ell}th field.
#'
#' Implementation of the Metropolis-Hastings update step for the distortion
#' probability \eqn{\alpha_\ell}. The default proposal distribution is a
#' reflected random walk.
#'
#' @family Metropolis-Hastings Update Functions
#'
#' @param current  A float. The current value of \eqn{\alpha_\ell} in the Markov
#'                 chain.
#' @param linkage  An integer vector. The cluster labels for each record.
#' @param field    A character vector. The values for the \eqn{\ell^{th}} field of
#'                 the records.
#' @param thetas_l A matrix. The first column of \code{thetas_l} must contain
#'                 all the possible values that the \eqn{\ell^{th}} field can
#'                 take (i.e. the unique values of \code{field}), and the second
#'                 column must contain the corresponding frequencies of each value
#'                 in the whole set of records.
#' @param a        A float. The \eqn{\alpha} shape parameter of the beta prior
#'                 for \eqn{\alpha_\ell}.
#' @param b        A float. The \eqn{\beta} shape parameter of the beta prior
#'                 for \eqn{\alpha_\ell}.
#' @param proposal A string. Specifies the proposal distribution to be used in
#'                 generating a new value for \eqn{\alpha_\ell}. Currently only
#'                 accepts "RRW" for reflected random walk.
#' @param drift    A float. A control parameter for tuning the variance of the
#'                 proposal distribution. larger values will lead to less
#'                 correlated proposals.
#'
#' @return The next value of \eqn{\alpha_\ell} in the Markov chain.
#'
alpha_metropolis <- function(current, linkage, field, thetas_l, a, b, proposal, drift){
  if (proposal == "RRW"){
    if (drift > 1) stop("Drift parameter must be in [0,1] for reflected random walk proposal.")

    # draw new alpha_l from proposal distribution
    new <- runif(1, min = current - drift, max = current + drift)
    new <- ifelse(new > 1, 2 - new, abs(new))

    # calculate acceptance probability
    log.r <- sum(log_distortion_cluster(linkage, field, new, thetas_l),
                 log(dbeta(new, a, b)),
                 -log_distortion_cluster(linkage, field, current, thetas_l),
                 -log(dbeta(current, a, b)))
  }
  # metropolis update step
  if(log(runif(1)) < log.r){
    return(new)
  } else {
    return(current)
  }
}

#' Metropolis-Hastings update of regression coefficients \eqn{\beta}.
#'
#' Implementation of the Metropolis-Hastings update step for the
#' regression coefficients \eqn{\beta}.
#'
#' @param linkage       An integer vector. The cluster labels for each record.
#' @param model_data    A data.frame. The model.frame of the regression. The
#'                      first column must contain the outcome variable \eqn{y}
#'                      and the remaining columns must contain the model
#'                      covariates \eqn{X}.
#' @param beta          A numeric vector. The current value of the regression
#'                      coefficients. Can also be a single float in the case of
#'                      univariate regression.
#' @param sigma_data    A matrix. The variance-covariance matrix of the
#'                      covariates. Must be symetric positive-definite. Can also
#'                      be a strictly positive float in the case of univariate
#'                      regression.
#' @param sigma_y       A float. The variance of the outcome variable \eqn{y}.
#'                      Must be strictly positive.
#' @param sigma_x       A matrix. The variance-covariance matrix of the
#'                      distortion of the covariates. Must be symetric
#'                      positive-definite. Can also be a strictly positive float
#'                      in the case of univariate regression.
#' @param beta_prior    A list. A list encoding the prior distribution of
#'                      \eqn{\beta}.
#' @param beta_proposal A list. A list encoding the proposal distribution to use
#'                      for \eqn{beta} in the Metropolis-Hastings step.
#'
#' @return An updated vector of the regression coefficients.
#'
beta_metropolis <- function(linkage, model_data, beta, sigma_data, sigma_y, sigma_x,
                            beta_prior, beta_proposal){
  if (beta_proposal[["proposal"]] == "MVN") {
    if (length(beta_proposal[["drift"]]) == 1) {
      S <- diag(rep(beta_proposal[["drift"]], length(beta)))
    } else {
      S <- diag(beta_proposal[["drift"]])
    }
    beta_new <- mvnfast::rmvn(n = 1, mu = beta, sigma = S)
    log.r <- sum(log_joint_lm_cluster(linkage, model_data, beta_new, sigma_data, sigma_y, sigma_x),
                 -log_joint_lm_cluster(linkage, model_data, beta, sigma_data, sigma_y, sigma_x))
    log.r <- log.r + log_beta_prior_dens(beta_new, beta_prior) - log_beta_prior_dens(beta, beta_prior)
    if(log(runif(1)) < log.r){
      return(beta_new)
    } else {
      return(beta)
    }
  } else if (beta_proposal[["proposal"]] == "normal") {
    for (i in 1:length(beta)) {
      if (length(beta_proposal[["drift"]]) == length(beta)) {
        s <- beta_proposal[["drift"]][i]
      } else {
        s <- beta_proposal[["drift"]]
      }
      beta_new <- beta
      beta_new[i] <- rnorm(n = 1, mean = beta[i], sd = s)
      log.r <- sum(log_joint_lm_cluster(linkage, model_data, beta_new, sigma_data, sigma_y, sigma_x),
                   -log_joint_lm_cluster(linkage, model_data, beta, sigma_data, sigma_y, sigma_x))
      log.r <- log.r + log_beta_prior_dens(beta_new, beta_prior) - log_beta_prior_dens(beta, beta_prior)
      if(log(runif(1)) < log.r){
        beta <- beta_new
      }
    }
    return(beta)
  }
}

#' Metropolis-Hastings update of the variance of \eqn{y}, \eqn{\sigma^2_y}.
#'
#' Implementation of the Metropolis-Hastings update step for the variance
#' of \eqn{y}, \eqn{\sigma^2_y}.
#'
#' @param linkage          An integer vector. The cluster labels for each
#'                         record.
#' @param model_data       A data.frame. The model.frame of the regression. The
#'                         first column must contain the outcome variable
#'                         \eqn{y} and the remaining columns must contain the
#'                         model covariates \eqn{X}.
#' @param beta             A numeric vector. The value of the regression
#'                         coefficients. Can also be a single float in the case
#'                         of univariate regression.
#' @param sigma_data       A matrix. The variance-covariance matrix of the
#'                         covariates. Must be symetric positive-definite. Can
#'                         also be a strictly positive float in the case of
#'                         univariate regression.
#' @param sigma_y          A float. The current variance of the outcome variable
#'                         \eqn{y}. Must be strictly positive.
#' @param sigma_x          A matrix. The variance-covariance matrix of the
#'                         distortion of the covariates. Must be symetric
#'                         positive-definite. Can also be a strictly positive
#'                         float in the case of univariate regression.
#' @param sigma_y_prior    A list. A list encoding the prior distribution of
#'                         \eqn{\sigma^2_y}.
#' @param sigma_y_proposal A list. A list encoding the proposal distribution to use
#'                         for \eqn{\sigma^2_y} in the Metropolis-Hastings step.
#'
#' @return An updated value of \eqn{\sigma^2_y}.
#'
sigma_y_metropolis <- function(linkage, model_data, beta, sigma_data, sigma_y, sigma_x,
                               sigma_y_prior, sigma_y_proposal){
  if (sigma_y_proposal[["proposal"]] == "folded-normal") {
    sigma_y_new <- VGAM::rfoldnorm(1, mean = sigma_y, sd = sigma_y_proposal[["drift"]])
    log.r <- sum(log_joint_lm_cluster(linkage, model_data, beta, sigma_data, sigma_y_new, sigma_x),
                 -log_joint_lm_cluster(linkage, model_data, beta, sigma_data, sigma_y, sigma_x),
                 log_sigma_y_prior_dens(sigma_y_new, sigma_y_prior),
                 -log_sigma_y_prior_dens(sigma_y, sigma_y_prior))
  }
  if(log(runif(1)) < log.r){
    return(sigma_y_new)
  } else {
    return(sigma_y)
  }
}

#' Metropolis-Hastings update of \eqn{\Sigma_{x|\widetilde{x}}}.
#'
#' Implementation of the Metropolis-Hastings update step for the covariate
#' distortion variance \eqn{\Sigma_{x|\widetilde{x}}}.
#'
#' @param linkage       An integer vector. The cluster labels for each record.
#' @param model_data    A data.frame. The model.frame of the regression. The
#'                      first column must contain the outcome variable \eqn{y}
#'                      and the remaining columns must contain the model
#'                      covariates \eqn{X}.
#' @param beta          A numeric vector. The value of the regression
#'                      coefficients. Can also be a single float in the case of
#'                      univariate regression.
#' @param sigma_data    A matrix. The variance-covariance matrix of the
#'                      covariates. Must be symetric positive-definite. Can also
#'                      be a strictly positive float in the case of univariate
#'                      regression.
#' @param sigma_y       A float. The variance of the outcome variable
#'                      \eqn{y}. Must be strictly positive.
#' @param sigma_x       A matrix. The current variance-covariance matrix of the
#'                      distortion of the covariates. Must be symetric
#'                      positive-definite. Can also be a strictly positive float
#'                      in the case of univariate regression.
#' @param sigma_x_prior A list. A list encoding the prior distribution of
#'                      \eqn{\Sigma_{x|\widetilde{x}}}.
#' @param sigma_x_proposal A list. A list encoding the proposal distribution to use
#'                         for \eqn{\Sigma_{x|\widetilde{x}}} in the
#'                         Metropolis-Hastings step.
#'
#' @return An updated value of \eqn{\Sigma_{x|\widetilde{x}}}.
#'
sigma_x_metropolis <- function(linkage, model_data, beta, sigma_data, sigma_y, sigma_x,
                               sigma_x_prior, sigma_x_proposal){
  P <- sqrt(length(sigma_x))
  if (sigma_x_proposal[["proposal"]] == "folded-normal") {
    if (length(sigma_x_proposal[["drift"]]) == 1) {
      drift <- rep(sigma_x_proposal[["drift"]], P)
    } else {
      drift <- sigma_x_proposal[["drift"]]
    }
    for (i in 1:P) {
      sigma_x_new <- sigma_x
      sigma_x_new[i,i] <- VGAM::rfoldnorm(1, mean = sigma_x[i,i], sd = drift)
      log.r <- sum(log_joint_lm_cluster(linkage, model_data, beta, sigma_data, sigma_y, sigma_x_new),
                   -log_joint_lm_cluster(linkage, model_data, beta, sigma_data, sigma_y, sigma_x),
                   log_sigma_x_prior_dens(sigma_x_new, sigma_x_prior),
                   -log_sigma_x_prior_dens(sigma_x, sigma_x_prior))
      if(log(runif(1)) < log.r){
        sigma_x <- sigma_x_new
      }
    }
    return(sigma_x)
  }
  if (sigma_x_proposal[["proposal"]] == "IW") {
    drift <- sigma_x_proposal[["drift"]]
    sigma_x_new <- CholWishart::rInvWishart(n = 1, df = drift, Sigma = (drift - P - 1) * sigma_x)[,,1]
    log.r <- sum(log_joint_lm_cluster(linkage, model_data, beta, sigma_data, sigma_y, sigma_x_new),
                 -log_joint_lm_cluster(linkage, model_data, beta, sigma_data, sigma_y, sigma_x),
                 CholWishart::dInvWishart(sigma_x,
                                          df = drift,
                                          Sigma = sigma_x_new * (drift - P - 1),
                                          log = TRUE),
                 -CholWishart::dInvWishart(sigma_x_new, df = drift,
                                           Sigma = sigma_x * (drift - P - 1),
                                           log = TRUE),
                 log_sigma_x_prior_dens(sigma_x_new, sigma_x_prior),
                 -log_sigma_x_prior_dens(sigma_x, sigma_x_prior))
    if(log(runif(1)) < log.r){
      return(sigma_x_new)
    } else {
      return(sigma_x)
    }
  }
}

#' Metropolis-Hastings update of the variance-covariance of \eqn{X}.
#'
#' Implementation of the Metropolis-Hastings update step for the
#' variance-covariance of the covariates \eqn{\Sigma_{\widetilde{x}}}.
#'
#' @param linkage       An integer vector. The cluster labels for each record.
#' @param model_data    A data.frame. The model.frame of the regression. The
#'                      first column must contain the outcome variable \eqn{y}
#'                      and the remaining columns must contain the model
#'                      covariates \eqn{X}.
#' @param beta          A numeric vector. The value of the regression
#'                      coefficients. Can also be a single float in the case of
#'                      univariate regression.
#' @param sigma_data    A matrix. The current variance-covariance matrix of the
#'                      covariates. Must be symetric positive-definite. Can also
#'                      be a strictly positive float in the case of univariate
#'                      regression.
#' @param sigma_y       A float. The variance of the outcome variable
#'                      \eqn{y}. Must be strictly positive.
#' @param sigma_x       A matrix. The variance-covariance matrix of the
#'                      distortion of the covariates. Must be symetric
#'                      positive-definite. Can also be a strictly positive float
#'                      in the case of univariate regression.
#' @param sigma_data_prior    A list. A list encoding the prior distribution of
#'                            \eqn{\Sigma_{\widetilde{x}}}.
#' @param sigma_data_proposal A list. A list encoding the proposal distribution to
#'                            use for \eqn{\Sigma_{\widetilde{x}}} in the
#'                            Metropolis-Hastings step.
#'
#' @return An updated value of \eqn{\Sigma_{\widetilde{x}}}.
#'
sigma_data_metropolis <- function(linkage, model_data, beta, sigma_data, sigma_y, sigma_x,
                                  sigma_data_prior, sigma_data_proposal){
  P <- sqrt(length(sigma_data))
  if (sigma_data_proposal[["proposal"]] == "folded-normal") {
    if (length(sigma_data_proposal[["drift"]]) == 1) {
      drift <- rep(sigma_data_proposal[["drift"]], P)
    } else {
      drift <- sigma_data_proposal[["drift"]]
    }
    for (i in 1:P) {
      sigma_data_new <- sigma_data
      sigma_data_new[i,i] <- VGAM::rfoldnorm(1, mean = sigma_data[i,i], sd = drift)
      log.r <- sum(log_joint_lm_cluster(linkage, model_data, beta, sigma_data_new, sigma_y, sigma_x),
                   -log_joint_lm_cluster(linkage, model_data, beta, sigma_data, sigma_y, sigma_x),
                   log_sigma_data_prior_dens(sigma_data_new, sigma_data_prior),
                   -log_sigma_data_prior_dens(sigma_data, sigma_data_prior))
      if(log(runif(1)) < log.r){
        sigma_data <- sigma_data_new
      }
    }
    return(sigma_data)
  }
  if (sigma_data_proposal[["proposal"]] == "IW") {
    drift <- sigma_data_proposal[["drift"]]
    sigma_data_new <- CholWishart::rInvWishart(n = 1, df = drift, Sigma = (drift - P - 1) * sigma_data)[,,1]
    log.r <- sum(log_joint_lm_cluster(linkage, model_data, beta, sigma_data_new, sigma_y, sigma_x),
                 -log_joint_lm_cluster(linkage, model_data, beta, sigma_data, sigma_y, sigma_x),
                 CholWishart::dInvWishart(sigma_data,
                                          df = drift,
                                          Sigma = sigma_data_new * (drift - P - 1) ,
                                          log = TRUE),
                 -CholWishart::dInvWishart(sigma_data_new,
                                           df = drift,
                                           Sigma = sigma_data * (drift - P - 1),
                                           log = TRUE),
                 log_sigma_x_prior_dens(sigma_data_new, sigma_data_prior),
                 -log_sigma_x_prior_dens(sigma_data, sigma_data_prior))
    if(log(runif(1)) < log.r){
      return(sigma_data_new)
    } else {
      return(sigma_data)
    }
  }
}

