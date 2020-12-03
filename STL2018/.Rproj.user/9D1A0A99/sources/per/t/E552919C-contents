#' Likelihood for the \eqn{\ell}th field of the hit-and-miss model.
#'
#' \code{hitmiss} gives the likelihood of the \eqn{\ell^{th}} field of a record
#' taking the value \code{v} given that the true value for the field is
#' \code{v_true}, the empirical frequency of \code{v} in all of the records is
#' \eqn{\theta_{v,\ell}}, and the distortion probability for the \eqn{\ell^{th}}
#' field is \eqn{\alpha_{\ell}}.
#'
#' @section Details:
#' This function is used in the Metropolis-Hastings update step for the linkage
#' structure \eqn{\Lambda}.
#'
#' Given that the true value for the \eqn{\ell^{th}} field is \code{v_true}, the
#' empirical frequency of \code{v} in all of the records is
#' \eqn{\theta_{v,\ell}}, and the distortion probability for the \eqn{\ell^{th}}
#' field is \eqn{\alpha_{\ell}}, then the likelihood of \code{v} is given by
#' \deqn{(1-\alpha_{\ell})\delta(v,v_{true}) + \alpha_{\ell} \theta_{v,\ell},}
#' where \eqn{\delta(v,v_{true}) = 1} if \eqn{v = v_{true}}, and
#' \eqn{\delta(v,v_{true}) = 0} otherwise.
#'
#' @family Cluster Likelihood Functions
#'
#' @param v      A string. The potential value for the \eqn{\ell^{th}} field of
#'               the record that \code{hitmiss} calculates the likelihood for.
#' @param v_true A string. The true value of the \eqn{\ell^{th}} field of the
#'               record.
#' @param theta  A float. The empirical frequency of \code{v} in the set of
#'               records.
#' @param alpha  A float. The distortion probability for the \eqn{\ell^{th}}
#'               field.
#'
#' @return The likelihood of \code{v}.
#'
hitmiss <- function(v, v_true, theta, alpha){
  return((1-alpha) * I(v == v_true) + alpha * theta)
}

#' Cluster likelihood for the \eqn{\ell}th field.
#'
#' Given the \eqn{\ell^{th}} field of a cluster as a vector along with the
#' empirical frequencies \eqn{\theta_\ell} and distortion probability,
#' \code{recursive_cluster} uses a recursive formula to calculate the joint
#' likelihood for the \eqn{q^{th}} cluster \eqn{v_{C_q}}.
#'
#' @section Details:
#' This function is used in the Metropolis-Hastings update step for the linkage
#' structure \eqn{\Lambda}.
#'
#' The recursive formula for \eqn{P(v_{C_{q},\ell} \;|\; \alpha_\ell, \lambda)}
#' is given by \deqn{P(v_{C_{q},\ell} \;|\; \alpha_\ell, \lambda) =
#' \theta_{v_{ij\ell}},} if \eqn{v_{ij\ell} = v_{C_q,\ell},} and
#'
#' \eqn{P(v_{C_{q},\ell} \;|\; \alpha_\ell, \lambda) =}
#' \deqn{(1-\alpha_\ell)\left[\prod_{(i', j') \in C_{q\setminus(ij)}} (1 -
#' \alpha_\ell) \delta(v_{i' j'\ell}, v_{ij\ell}) + \alpha_\ell \theta_{v_{i'
#' j'},\ell}\right] \theta_{v_{ij\ell}} + \alpha_\ell \theta_{v_{ij\ell}}
#' P(V_{C_{q}\setminus(ij),\ell} \; | \; \alpha_\ell, \lambda)} otherwise.
#'
#' @family Cluster Likelihood Functions
#'
#' @param cluster  A vector. The vector of strings/factors for the
#'                 \eqn{\ell^{th}} field of the selected cluster.
#' @param thetas_l A matrix. The first column of \code{thetas_l} must be the
#'                 possible values that \code{cluster}, and the second column
#'                 should contain the corresponding frequencies of each value
#'                 in the whole set of records.
#' @param alpha_l  A float. The distortion probability of the \eqn{\ell^{th}}
#'                 field.
#'
#' @return The joint likelihood of the \eqn{\ell^{th}} field of a cluster.
#'
recursive_cluster <- function(cluster, thetas_l, alpha_l){
  if(length(cluster) == 0) return(1)
  v_true <- cluster[1]
  cluster <- cluster[-1]
  theta <- thetas_l[thetas_l[,1] == v_true, 2]
  if(length(cluster) == 0){
    return(theta)
  } else {
    likelihoods <- mapply(hitmiss,
                          v = cluster,
                          theta = thetas_l[match(cluster, thetas_l[,1]), 2],
                          MoreArgs = list(v_true = v_true, alpha = alpha_l))
    prob <- (1 - alpha_l) * prod(likelihoods) * theta
    return(prob + alpha_l * theta * recursive_cluster(cluster, thetas_l, alpha_l))
  }
}

#' Cluster log-likelihood.
#'
#' Given a data.frame of records for a specific cluster,
#' \code{log_recursive_cluster} loops over each key variable (each column of the
#' cluster), calls \code{recursive_cluster}, and returns the sums of the
#' log-likelihoods for each field.
#'
#' @section Details:
#' This function is used in the Metropolis-Hastings update step for the linkage
#' structure \eqn{\Lambda}.
#'
#' @seealso \code{\link{recursive_cluster}} for the actual recursion formula
#'   implemented for this function.
#'
#' @family Cluster Likelihood Functions
#'
#' @param cluster A data.frame. All of the records and their values for each key
#'                variable for a specific cluster.
#' @param thetas  A list of matrices. Each matrix corresponds to a key variable
#'                and contains two columns: 1) all of the values each key variable
#'                can take in the data, and 2) the corresponding empirical
#'                frequency of each value.
#' @param alpha   A data.frame. Contains a column \code{key_vars} with each of the
#'                key variables, and a column \code{alpha} with the corresponding
#'                distortion probabilities.
#'
#' @return The joint log-likelihood of the provided cluster.
#'
log_linkage_cluster <- function(cluster, thetas, alphas, key_vars){
  log_lik <- 0
  for (key_var in key_vars){
    alpha_l <- alphas[alphas$key_vars == key_var, 'alpha']
    cluster_l <- cluster[[key_var]]
    thetas_l <- thetas[[key_var]][thetas[[key_var]][,1] %in% cluster_l, ]
    log_lik <- log_lik + log(recursive_cluster(cluster_l, thetas_l, alpha_l))
  }
  return(log_lik)
}

#' Joint cluster likelihood for the \eqn{\ell}th field.
#'
#' Given all the values for the \eqn{\ell^{th}} field of the observed records in
#' \code{field}, as well as the cluster membership for each record in
#' \code{linkage}, \code{log_alpha_cluster} returns the joint log-likelihood of
#' all the clusters for the \eqn{\ell^{th}} field.
#'
#' @section Details:
#' This function is used in the Metropolis-Hastings update step for the
#' distortion probabilities \eqn{\alpha}.
#'
#' @seealso \code{\link{recursive_cluster}} for the actual recursion formula
#'   implemented for this function.
#'
#' @family Cluster Likelihood Functions
#'
#' @param linkage  An integer vector. The cluster labels for each record.
#' @param field    A character vector. The values for the \eqn{\ell^{th}}
#'                 field of the records.
#' @param alpha_l  A float. The distortion probability for the \eqn{\ell^{th}}
#'                 field.
#' @param thetas_l A matrix. The first column of \code{thetas_l} must contain
#'                 all the possible values that the \eqn{\ell^{th}} field can
#'                 take (i.e. the unique values of \code{field}), and the second
#'                 column must contain the corresponding frequencies of each value
#'                 in the whole set of records.
#'
#' @return The joint cluster log-likelihood for the \eqn{\ell^{th}} field.
#'
log_distortion_cluster <- function(linkage, field, alpha_l, thetas_l){
  log_lik <- 0
  for (link in linkage){
    cluster_l <- field[linkage == link]
    thetas_field <- thetas_l[thetas_l[,1] %in% field, ]
    log_lik <- log_lik + log(recursive_cluster(cluster_l, thetas_field, alpha_l))
  }
  return(log_lik)
}

#' Joint conditional regression log-likelihood.
#'
#' Conditional on the linkage structure \eqn{\Lambda} (\code{linkage}) and the
#' regression parameters \eqn{\beta} (\code{beta}), \eqn{\Sigma_{\widetilde{x}}
#' (\code{sigma_data}), \eqn{\Sigma_{x|\widetilde{x}} (\code{sigma_x}), and
#' \eqn{\sigma^2_{y|\widetilde{x}}} (\code{sigma_y}),
#' \code{log_joint_lm_cluster} calculates the joint log-likelihood of the
#' covariates \eqn{X} and the outcome variable \eqn{y} over all of the clusters.
#'
#' @section Details:
#' Let \eqn{C_j} be the cluster of records corresponding to the \eqn{j^{th}} true
#' entity, then assuming that there are \eqn{N_{pop}} true entities represented
#' in the records (i.e. there are \eqn{N_{pop}} unique elements of
#' \code{linkage}), \code{log_joint_lm_cluster} calculates
#' \deqn{\sum_{j=1}^{N_{pop}} \log P\left([y,X]_{C_j} \; | \; \Lambda, \beta,
#' \Sigma_{\widetilde{x}}, \Sigma_{x|\widetilde{x}},
#' \sigma^2_{y|\widetilde{x}}\right).}
#'
#' @family Cluster Likelihood Functions
#'
#' @param linkage    An integer vector. The cluster labels for each record.
#' @param model_data A data.frame. The model.frame of the regression. The first
#'                   column must contain the outcome variable \eqn{y} and the
#'                   remaining columns must contain the model covariates \eqn{X}.
#' @param beta       A numeric vector. The regression coefficients. Can also be
#'                   a single float in the case of univariate regression.
#' @param sigma_data A matrix. The variance-covariance matrix of the covariates.
#'                   Must be symetric positive-definite. Can also be a strictly
#'                   positive float in the case of univariate regression.
#' @param sigma_y    A float. The variance of the outcome variable \eqn{y}. Must
#'                   be strictly positive.
#' @param sigma_x    A matrix. The current variance-covariance matrix of the
#'                   distortion of the covariates. Must be symetric
#'                   positive-definite. Can also be a strictly positive float
#'                   in the case of univariate regression.
#'
#' @return The joint log-likelihood of \eqn{(y,X)}.
#'
log_joint_lm_cluster <- function(linkage, model_data, beta, sigma_data, sigma_y, sigma_x){
  bSb <- t(beta) %*% sigma_data %*% beta
  bS <- t(beta) %*% sigma_data
  B <- rbind(c(bSb, bS), cbind(t(bS), sigma_data))
  Sigma <- Matrix::as.matrix(Matrix::bdiag(sigma_y, sigma_x))

  log_lik <- 0

  for (link in unique(linkage)){
    # extract cluster and stack rows into a single vector
    cluster <- model_data[link == linkage, ]
    log_lik <- log_lik + log_lm_cluster(cluster, B, Sigma)
  }
  return(log_lik)
}

#' Conditional regression log-likelihood for single cluster.
#'
#' Conditional on the  the regression parameters \eqn{\beta} (\code{beta}),
#' \eqn{\Sigma_{\widetilde{x}} (\code{sigma_data}),
#' \eqn{\Sigma_{x|\widetilde{x}} (\code{sigma_x}), and
#' \eqn{\sigma^2_{y|\widetilde{x}}} (\code{sigma_y}),
#' \code{log_lm_cluster} calculates the log-likelihood of the
#' covariates \eqn{X} and the outcome variable \eqn{y} for the \eqn{j^{th}}
#' cluster \eqn{C_j} in the records.
#'
#' @section Details:
#' Let \eqn{C_j} be the cluster of records corresponding to the \eqn{j^{th}}
#' true entity, then assuming that there are \eqn{n} records in \eqn{C_j}, and
#' \eqn{p} covariates in \eqn{X}, \code{log_lm_cluster} calculates
#' \deqn{\log P\left([y,X]_{C_j} \; | \; \beta, \Sigma_{\widetilde{x}},
#' \Sigma_{x|\widetilde{x}}, \sigma^2_{y|\widetilde{x}}\right) =
#' \log N_{n(p+1)}\left(\mathbf{0}, \; (\mathbf{1}_n \mathbf{1}_n') \otimes
#' \mathbf{B} + \mathbf{I}_{n} \otimes \Sigma \right).}
#' where
#' \eqn{\mathbf{B}=}
#' \tabular{rccl}{
#' | \tab \eqn{\beta' \Sigma_{\widetilde{x}} \beta} \tab \eqn{\beta'
#' \Sigma_{\widetilde{x}}} \tab | \cr
#' | \tab \eqn{\Sigma_{\widetilde{x}} \beta} \tab \eqn{\Sigma_{\widetilde{x}}}
#' \tab |
#' }
#' and
#' \eqn{\Sigma=}
#' \tabular{lccr}{
#' | \tab \eqn{\sigma^2_{y|\widetilde{x}}} \tab \eqn{\mathbf{0}} \tab | \cr
#' | \tab \eqn{\mathbf{0}} \tab \eqn{\Sigma_{x|\widetilde{x}}} \tab |.
#' }
#'
#' @family Cluster Likelihood Functions
#'
#' @param cluster A data.fram. The model.frame of the regression corresponding
#'                to the records of a particular cluster. The first column must
#'                contain the outcome variable \eqn{y} and the remaining columns
#'                must contain the model covariates \eqn{X}.
#' @param B       A matrix. The first block component of the variance-covariance
#'                matrix.
#' @param Sigma   A matrix. The first block component of the variance-covariance
#'                matrix.
#'
#' @return The conditional log-likelihood of the given cluster.
#'
log_lm_cluster <- function(cluster, B, Sigma){

  if (nrow(cluster) == 0) return(0)

  # number of records in cluster
  n <- nrow(cluster)

  # flatten cluster into vector
  cluster <- c(t(cluster))

  # creating variance covariance matrix
  if (n == 1) {
    varcov <- B + Sigma
  } else {
    diagonal <- diag(n) %x% Sigma
    varcov <- diagonal + rep(1,n) %*% t(rep(1,n)) %x% B
  }

  # marginalizing out missing values
  if (any(is.na(cluster))){
    varcov <- varcov[!is.na(cluster), !is.na(cluster)]
    cluster <- cluster[!is.na(cluster)]
  }
  mn <- rep(0, length(cluster))
  log_lik <- mvnfast::dmvn(cluster, mu = mn, sigma = varcov, log = TRUE)
  return(log_lik)
}

#' Log prior probability of linkage struture
#'
#' \code{log_p_lambda} calculates \eqn{\log P(\lambda_{ij} \;|\;
#' \lambda_{-(ij)})} where \eqn{\lambda_{ij}} is the cluster label of the
#' \eqn{j^{th}} record in the \eqn{i^{th}} database.
#'
#' @param lambda       An integer. The cluster label corresponding to the record
#'                     in the \code{record} row of the data.
#' @param record       An integer. The row corresponding to the current record
#'                     being updated in the data.
#' @param linkage      An integer vector. The current cluster labels for each
#'                     record.
#' @param lambda_prior A list. A list containing the type of prior to use for
#'                     the linkage structure and the corresponding parameters
#'                     for the prior.
#'
#' @return \eqn{\log P(\lambda_{ij} \;|\; \lambda_{-(ij)})}
#'
log_p_lambda <- function(lambda, record, linkage, lambda_prior){

  N <- length(linkage)
  k_ij <- dplyr::n_distinct(linkage[-record])
  linkage_ij <- linkage[-record]

  if (lambda_prior$prior == "uniform") {
    if (lambda %in% linkage_ij) {
      return(-log(N))
    } else {
      return(log(N - k_ij) - log(N))
    }
  } else if (lambda_prior$prior == "PYP") {
    sigma <- lambda_prior$sigma
    nu <- lambda_prior$nu
    if (lambda %in% linkage_ij) {
      n_q <- sum(linkage_ij == lambda)
      return(log(n_q - sigma) - log(N - 1 + nu))
    } else {
      return(log(k_ij * sigma + nu) - log(N - 1 + nu))
    }
  } else {
    stop(paste(lambda_prior$prior, "is not an implemented prior for the linkage structure."))
  }
}
