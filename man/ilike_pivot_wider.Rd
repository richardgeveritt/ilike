% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_output.R
\name{ilike_pivot_wider}
\alias{ilike_pivot_wider}
\title{Take output from the functions load_mcmc_output, load_smc_output or load_enk_output and convert it to a more standard format (non-tidy format with one column per dimension).}
\usage{
ilike_pivot_wider(output, variables = NULL)
}
\arguments{
\item{output}{Output from the functions load_mcmc_output, load_smc_output or load_enk_output.}

\item{variables}{(optional) Variables to include in the output (default is all).}
}
\description{
Take output from the functions load_mcmc_output, load_smc_output or load_enk_output and convert it to a more standard format (non-tidy format with one column per dimension).
}
