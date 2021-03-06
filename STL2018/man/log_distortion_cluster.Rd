% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihood_functions.R
\name{log_distortion_cluster}
\alias{log_distortion_cluster}
\title{Joint cluster likelihood for the \eqn{\ell}th field.}
\usage{
log_distortion_cluster(linkage, field, alpha_l, thetas_l)
}
\arguments{
\item{linkage}{An integer vector. The cluster labels for each record.}

\item{field}{A character vector. The values for the \eqn{\ell^{th}}
field of the records.}

\item{alpha_l}{A float. The distortion probability for the \eqn{\ell^{th}}
field.}

\item{thetas_l}{A matrix. The first column of \code{thetas_l} must contain
all the possible values that the \eqn{\ell^{th}} field can
take (i.e. the unique values of \code{field}), and the second
column must contain the corresponding frequencies of each value
in the whole set of records.}
}
\value{
The joint cluster log-likelihood for the \eqn{\ell^{th}} field.
}
\description{
Given all the values for the \eqn{\ell^{th}} field of the observed records in
\code{field}, as well as the cluster membership for each record in
\code{linkage}, \code{log_alpha_cluster} returns the joint log-likelihood of
all the clusters for the \eqn{\ell^{th}} field.
}
\section{Details}{

This function is used in the Metropolis-Hastings update step for the
distortion probabilities \eqn{\alpha}.
}

\seealso{
\code{\link{recursive_cluster}} for the actual recursion formula
  implemented for this function.

Other Cluster Likelihood Functions: 
\code{\link{hitmiss}()},
\code{\link{log_joint_lm_cluster}()},
\code{\link{log_linkage_cluster}()},
\code{\link{log_lm_cluster}()},
\code{\link{recursive_cluster}()}
}
\concept{Cluster Likelihood Functions}
