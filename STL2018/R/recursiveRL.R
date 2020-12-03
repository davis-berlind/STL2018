#' Title
#'
#' @importFrom magrittr `%>%`
#' @importFrom foreach `%dopar%`
#'
#' @param data
#' @param key_vars
#' @param sample_size
#' @param burnin
#' @param thin
#' @param n_chains
#' @param progressbar
#' @param alpha_init
#' @param alpha_prior
#' @param alpha_proposal
#' @param alpha_drift
#' @param lambda_init
#' @param lambda_prior
#'
#' @return
#' @export
#'
#' @examples
recursiveRL <- function(data, key_vars, sample_size, burnin = 0, thin = 1,
                        n_chains = 1, progressbar = TRUE,
                        alpha_init = NULL, alpha_prior, alpha_proposal = NULL,
                        lambda_init = NULL, lambda_prior) {

  #### Argument Tests ####

  # creating index of records and databases
  if (is.data.frame(data)) {
    index <- list(recordID = 1:nrow(data))
  } else if (is.list(data) & all(sapply(data, is.data.frame))) {
    index <- list(dfID     = rep(1:length(data), times = sapply(data, nrow)),
                  recordID = Reduce(c, sapply(sapply(data, nrow), seq)))
  } else {
    stop("Argument data must be a data.frame or a list of data.frame objects.")
  }

  data <- data_tests(data, key_vars)
  sampler_tests(key_vars, burnin, sample_size, thin,
                n_chains, progressbar,
                alpha_init, lambda_init)

  #### prior tests ####
  lambda_prior_tests(lambda_prior)
  alpha_prior_tests(alpha_prior, key_vars)
  if (is.null(alpha_proposal)) {
    # default to reflected random walk with 0.1 drift if no proposal given
    alpha_proposal <- list(proposal = "RRW", drift = 0.1)
  } else {
    alpha_proposal_tests(alpha_proposal, key_vars)
  }

  #### Initializing Sampler ####

  iter <- burnin + sample_size * thin
  N <- nrow(data)       # number of records
  L <- length(key_vars) # number of fields

  # calculate and store thetas
  thetas <- emp_dist(data, key_vars)

  # creating IDs for each record
  data$RLid <- seq(N)

  # initialize linkage structure
  lambda_chain <- matrix(nrow = sample_size, ncol = N)

  # creating column to track current linkage structure
  if (is.null(lambda_init)) {
    data$linkage <- 1:N
  } else {
    data$linkage <- lambda_init
  }

  # initialize distortion probabilities
  alpha_chain <- matrix(nrow = sample_size, ncol = L)

  # creating column to track current distortion probabilities
  if (is.null(alpha_init)) {
    alphas <- data.frame(key_vars = key_vars,
                         alpha = rep(0.1, L))
  } else {
    alphas <- data.frame(key_vars = key_vars,
                         alpha = alpha_init)
  }

  #### sampler ####

  # setting up parallels
  numcores <- min(n_chains, parallel::detectCores())
  if (numcores > 1) {
    # only register parallel backend if n_chains > 1
    cl <- parallel::makeCluster(numcores)
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }

  ret <- foreach::foreach(chain = 1:n_chains, .combine = list) %dopar% {
    if(!exists("pb") & progressbar){
      pb <- tcltk::tkProgressBar(paste("Chain", chain), min = 1, max = iter)
    }
    for (i in 1:iter) {

      if(progressbar) tcltk::setTkProgressBar(pb, i)

      #### Lambda update step ####

      for (record in 1:N) {

        lambda_old <- data[record, "linkage"]  # get current cluster assignment
        lambda_new <- sample.int(N, size = 1)  # sample new cluster assignment

        old_cluster <- data[data$linkage == lambda_old, ]                       # get cluster of ij^th record under current assignment
        new_cluster <- data[data$linkage == lambda_new | data$RLid == record, ] # get cluster of ij^th record under new assignment
        old_cluster_ij <- old_cluster[old_cluster$RLid != record, ]             # old cluster without ij^th record
        new_cluster_ij <- new_cluster[new_cluster$RLid != record, ]             # new cluster without ij^th record

        # calculate update probability
        log.r <- sum(log_linkage_cluster(new_cluster, thetas, alphas, key_vars),
                     log_linkage_cluster(old_cluster_ij, thetas, alphas, key_vars),
                     log_p_lambda(lambda_new, record, data$linkage, lambda_prior),
                     -log_linkage_cluster(old_cluster, thetas, alphas, key_vars),
                     -log_linkage_cluster(new_cluster_ij, thetas, alphas, key_vars),
                     -log_p_lambda(lambda_new, record, data$linkage, lambda_prior))

        # metropolis update step
        if (log(runif (1)) < log.r) {
          data[record, "linkage"] <- lambda_new
        }
      }

      #### alpha update step ####

      for (key_var in key_vars) {
        if (key_var %in% names(alpha_proposal)) {
          alpha_proposal <- alpha_proposal[[key_var]][["proposal"]]
          alpha_drift <- alpha_proposal[[key_var]][["drift"]]
        } else {
          alpha_proposal <- alpha_proposal[["proposal"]]
          alpha_drift <- alpha_proposal[["drift"]]
        }


        alpha_new <- alpha_metropolis(current  = alphas[alphas$key_vars == key_var, "alpha"],
                                      linkage  = data[["linkage"]],
                                      field    = data[[key_var]],
                                      thetas_l = thetas[[key_var]],
                                      a        = alpha_prior[[key_var]][["alpha"]],
                                      b        = alpha_prior[[key_var]][["beta"]],
                                      proposal = alpha_proposal,
                                      drift    = alpha_drift)

        alphas[alphas$key_vars == key_var, "alpha"] <- alpha_new
      }

      # store update if outside burnin and thinning interval
      if ((i > burnin) & ((i - burnin) %% thin == 0)) {
        dex <- (i - burnin) / thin
        lambda_chain[dex, ] <- data$linkage
        alpha_chain[dex, ] <- alphas$alpha
      }
    }
    ret <- list(lambda = lambda_chain,
                alpha  = alpha_chain)
    ret
  }
  on.exit(if(exists("cl")) parallel::stopCluster(cl))
  if(n_chains > 1) names(ret) <- paste("Chain",seq(n_chains))
  ret[["index"]] <- index
  return(ret)
}



