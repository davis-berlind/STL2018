#' Title
#'
#' @param formula
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
#' @param lambda_init
#' @param lambda_prior
#' @param beta_init
#' @param beta_prior
#' @param beta_proposal
#' @param sigma_y_init
#' @param sigma_y_prior
#' @param sigma_y_proposal
#' @param sigma_x_diag
#' @param sigma_x_init
#' @param sigma_x_prior
#' @param sigma_x_proposal
#' @param sigma_data_fixed
#' @param sigma_data_diag
#' @param sigma_data_init
#' @param sigma_data_prior
#' @param sigma_data_proposal
#'
#' @return
#' @export
#'
#' @examples
regressionRL <- function(formula, data, key_vars,
                         sample_size, burnin = 0, thin = 1,
                         n_chains = 1, progressbar = TRUE,
                         alpha_init = NULL, alpha_prior, alpha_proposal = NULL,
                         lambda_init = NULL, lambda_prior,
                         beta_init = NULL, beta_prior = NULL, beta_proposal = NULL,
                         sigma_y_init = NULL, sigma_y_prior = NULL, sigma_y_proposal = NULL,
                         sigma_x_diag = FALSE,
                         sigma_x_init = NULL, sigma_x_prior = NULL, sigma_x_proposal = NULL,
                         sigma_data_fixed = TRUE, sigma_data_diag = FALSE,
                         sigma_data_init = NULL, sigma_data_prior = NULL, sigma_data_proposal = NULL){

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

  if (is.null(beta_prior)) {
    # default to flat prior
    beta_prior <- list(prior = "flat")
  } else {
    # TODO: implement beta prior tests
  }

  if (is.null(sigma_y_prior)) {
    # default to flat prior
    sigma_y_prior <- list(prior = "flat")
  } else {
    # TODO: implement sigma_y prior tests
  }

  if (is.null(sigma_x_prior)) {
    # default to flat prior
    sigma_x_prior <- list(prior = "flat")
  } else {
    # TODO: implement sigma_x prior tests
  }

  if (!sigma_data_fixed) {
    if (is.null(sigma_data_prior)) {
      # default to flat prior
      sigma_data_prior <- list(prior = "flat")
    } else {
      # TODO: implement sigma_data prior tests
    }
  }

  #### Initializing Sampler ####

  # creating model matrix
  rdata <- model.frame(formula, data = data, na.action = NULL)
  if (!all(apply(rdata, 2, is.numeric))) {
    stop("regressionRL does not currently support regressison with factor covariates.")
  }

  iter <- burnin + sample_size * thin
  N <- nrow(data)       # number of records
  L <- length(key_vars) # number of fields
  P <- ncol(rdata) - 1  # number of covariates

  # calculate and store thetas
  thetas <- emp_dist(data, key_vars)

  # estimate and store variance of covariates
  sigma_data <- var(rdata[,-1], na.rm = TRUE)

  if (all(is.na(sigma_data))) {
    warning("Broken Regression: No overlapping covariates. Assuming covariates are jointly independent.")
    sigma_data <- diag(apply(rdata[,-1], 2, var, na.rm = TRUE))
  } else if (sigma_x_diag) {
    sigma_data <- diag(apply(rdata[,-1], 2, var, na.rm = TRUE))
  } else if (any(is.na(sigma_data))) {
    warning("Broken Regression: Some of the covariates do not overlap. Assuming independence for non-overlaping covariates.")
    for (i in 1:P){
      for (j in 1:P){
        if (i == j & is.na(sigma_data[i,j])){
          sigma_data[i,j] <- var(rdata[,i+1], na.rm = TRUE)
        }
        else if (is.na(sigma_data[i, j])){
          sigma_data[i,j] <- 0
        }
      }
    }
  }

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

  if (is.null(alpha_proposal)) {
    # default to reflected random walk with 0.1 drift if no proposal given
    alpha_proposal <- list(proposal = "RRW", drift = 0.1)
  } else {
    alpha_proposal_tests(alpha_proposal, key_vars)
  }

  # creating column to track current distortion probabilities
  if (is.null(alpha_init)) {
    alphas <- data.frame(key_vars = key_vars,
                         alpha = rep(0.1, L))
  } else {
    alphas <- data.frame(key_vars = key_vars,
                         alpha = alpha_init)
  }

  # initial OLS w/o intercept
  ols <- rdata %>%
    dplyr::mutate_all(~ifelse(is.na(.), mean(., na.rm = TRUE), .)) %>%
    lm(update(formula, ~ . -1), data = .)

  resvar <- var(resid(ols))    # get variance of ols residuals
  se <- sqrt(diag(vcov(ols)))  # get standard errors of ols coefficients

  # intialize beta
  beta_chain <- matrix(nrow = sample_size, ncol = P)

  if (is.null(beta_proposal)) {
    # default to normal proposal w/ 95% CI drift param
    beta_proposal <- list(proposal = "normal", drift = 2*se)
  } else {
    # TODO: implement beta proposal tests
  }

  if(is.null(beta_init)) {
    beta <- ols %>% coef %>% as.vector
  } else {
    beta <- beta_init
  }

  # intialize sigma_x
  if (sigma_x_diag) {
    # only store diagonal elements
    sigma_x_chain <- matrix(nrow = sample_size, ncol = P)
  } else {
    # store lower trianular elements
    sigma_x_chain <- matrix(nrow = sample_size, ncol = P*(P + 1) / 2)
  }

  if (is.null(sigma_x_proposal) & sigma_x_diag) {
    # default to folded-normal proposal
    sigma_x_proposal <- list(proposal = "folded-normal", drift = 0.1)
  } else if (is.null(sigma_x_proposal) & !sigma_x_diag) {
    sigma_x_proposal <- list(proposal = "IW", drift = 1000)
  } else {
    # TODO: implement sigma_x proposal tests
  }

  if(is.null(sigma_x_init)) {
    sigma_x <- diag(rep(0.1, P))
  } else {
    sigma_x <- sigma_x_init
  }

  # intialize sigma_y
  sigma_y_chain <- rep(NA, sample_size)

  if (is.null(sigma_y_proposal)) {
    # default to folded-normal proposal
    sigma_y_proposal <- list(proposal = "folded-normal",
                             drift = ifelse(resvar <= 1, sqrt(resvar), sqrt(log(resvar))))
  } else {
    # TODO: implement sigma_y proposal tests
  }

  if(is.null(sigma_y_init)) {
    sigma_y <- var(resid(ols))
  } else {
    sigma_y <- sigma_y_init
  }

  # intialize sigma_data
  if (sigma_data_diag & !sigma_data_fixed) {
    # only store diagonal elements
    sigma_data_chain <- matrix(nrow = sample_size, ncol = P)
  } else if (!sigma_data_fixed) {
    # store lower trianular elements
    sigma_data_chain <- matrix(nrow = sample_size, ncol = P*(P + 1) / 2)
  }

  if (is.null(sigma_data_proposal) & sigma_data_diag) {
    # default to folded-normal proposal
    sigma_data_proposal <- list(proposal = "folded-normal",
                                drift = ifelse(diag(sigma_data) <= 1,
                                               sqrt(diag(sigma_data)),
                                               sqrt(log(diag(sigma_data))))
                                )
  } else if (is.null(sigma_data_proposal) & !sigma_data_diag) {
    sigma_data_proposal <- list(proposal = "IW", drift = 1000)
  } else {
    # TODO: implement sigma_x proposal tests
  }

  if(!is.null(sigma_data_init)) {
    sigma_data <- sigma_data_init
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

      bSb <- t(beta) %*% sigma_data %*% beta
      bS <- t(beta) %*% sigma_data
      B <- rbind(c(bSb, bS), cbind(t(bS), sigma_data))
      Sigma <- Matrix::bdiag(sigma_y, sigma_x)

      for (record in 1:N) {

        lambda_old <- data[record, "linkage"]  # get current cluster assignment
        lambda_new <- sample.int(N, size = 1)  # sample new cluster assignment

        old_cluster <- data[data$linkage == lambda_old, ]                       # get cluster of ij^th record under current assignment
        new_cluster <- data[data$linkage == lambda_new | data$RLid == record, ] # get cluster of ij^th record under new assignment
        old_cluster_ij <- old_cluster[old_cluster$RLid != record, ]             # old cluster without ij^th record
        new_cluster_ij <- new_cluster[new_cluster$RLid != record, ]             # new cluster without ij^th record

        old_rcluster <- rdata[data$linkage == lambda_old, ]                       # get (y,X) for cluster of ij^th record under current assignment
        new_rcluster <- rdata[data$linkage == lambda_new | data$RLid == record, ] # get (y,X) for cluster of ij^th record under new assignment
        old_rcluster_ij <- old_rcluster[old_cluster$RLid != record, ]             # (y,X) of old cluster without ij^th record
        new_rcluster_ij <- new_rcluster[new_cluster$RLid != record, ]             # (y,X) of new cluster without ij^th record

        # calculate update probability
        log.r <- sum(log_linkage_cluster(new_cluster, thetas, alphas, key_vars),
                     log_linkage_cluster(old_cluster_ij, thetas, alphas, key_vars),
                     log_lm_cluster(new_rcluster, B, Sigma),
                     log_lm_cluster(old_rcluster_ij, B, Sigma),
                     log_p_lambda(lambda_new, record, data$linkage, lambda_prior),
                     -log_linkage_cluster(old_cluster, thetas, alphas, key_vars),
                     -log_linkage_cluster(new_cluster_ij, thetas, alphas, key_vars),
                     -log_lm_cluster(old_rcluster, B, Sigma),
                     -log_lm_cluster(new_rcluster_ij, B, Sigma),
                     -log_p_lambda(lambda_old, record, data$linkage, lambda_prior))

        # metropolis update step
        if (log(runif (1)) < log.r) {
          data[record, "linkage"] <- lambda_new
        }
      }

      #### alpha update step ####

      for (key_var in key_vars) {
        if (key_var %in% names(alpha_proposal)) {
          proposal <- alpha_proposal[[key_var]][["proposal"]]
          drift <- alpha_proposal[[key_var]][["drift"]]
        } else {
          proposal <- alpha_proposal[["proposal"]]
          drift <- alpha_proposal[["drift"]]
        }

        alpha_new <- alpha_metropolis(current  = alphas[alphas$key_vars == key_var, "alpha"],
                                      linkage  = data[["linkage"]],
                                      field    = data[[key_var]],
                                      thetas_l = thetas[[key_var]],
                                      a        = alpha_prior[[key_var]][["alpha"]],
                                      b        = alpha_prior[[key_var]][["beta"]],
                                      proposal = proposal,
                                      drift    = drift)

        alphas[alphas$key_vars == key_var, "alpha"] <- alpha_new
      }

      ##### beta update step #####
      beta <- beta_metropolis(data$linkage, rdata, beta, sigma_data, sigma_y, sigma_x,
                              beta_prior, beta_proposal)

      ##### sigma_y update step #####
      sigma_y <- sigma_y_metropolis(data$linkage, rdata, beta, sigma_data, sigma_y, sigma_x,
                                    sigma_y_prior, sigma_y_proposal)

      ##### sigma_x update step #####
      sigma_x <- sigma_x_metropolis(data$linkage, rdata, beta, sigma_data, sigma_y, sigma_x,
                                    sigma_x_prior, sigma_x_proposal)

      ##### sigma_data update step #####
      if (!sigma_data_fixed) {
        sigma_data <- sigma_data_metropolis(data$linkage, rdata, beta, sigma_data, sigma_y, sigma_x,
                                            sigma_data_prior, sigma_data_proposal)
      }

      # store update if outside burnin and thinning interval
      if ((i > burnin) & ((i - burnin) %% thin == 0)) {
        dex <- (i - burnin) / thin
        lambda_chain[dex, ] <- data$linkage
        alpha_chain[dex, ] <- alphas$alpha
        beta_chain[dex, ] <- beta
        sigma_y_chain[dex] <- sigma_y
        if (sigma_x_diag) {
          sigma_x_chain[dex, ] <- diag(sigma_x)
        } else {
          sigma_x_chain[dex, ] <- MCMCpack::vech(sigma_x)
        }
        if (!sigma_data_fixed & sigma_data_diag) {
          sigma_data_chain[dex, ] <- diag(sigma_data)
        } else if (!sigma_data_fixed) {
          sigma_data_chain[dex, ] <- MCMCpack::vech(sigma_data)
        }
      }


    }
    ret <- list(lambda  = lambda_chain,
                alpha   = alpha_chain,
                beta    = beta_chain,
                sigma_y = sigma_y_chain,
                sigma_x = sigma_x_chain)
    if (!sigma_data_fixed) {
      ret[["sigma_data"]] <- sigma_data_chain
    }
    ret
  }
  on.exit(if(exists("cl")) parallel::stopCluster(cl))
  if(n_chains > 1) names(ret) <- paste("Chain",seq(n_chains))
  ret[["index"]] <- index
  return(ret)
}
