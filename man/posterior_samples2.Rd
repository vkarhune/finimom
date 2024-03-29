% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior_samples2.R
\name{posterior_samples2}
\alias{posterior_samples2}
\title{Title}
\usage{
posterior_samples2(
  beta,
  se,
  eaf,
  R,
  maxsize,
  tau0,
  r0,
  niter,
  burnin,
  p,
  seed = 456,
  excl.burnin = TRUE,
  n,
  a0 = 0.05,
  b0 = 0.95,
  inds0 = NULL,
  standardize = TRUE,
  msprior = NULL,
  verbose = TRUE,
  clump = TRUE,
  clump_r2 = 0.99^2,
  check_ld = FALSE,
  anchor = FALSE,
  anchor_cslevel = 0.5
)
}
\arguments{
\item{beta}{Vector of effect sizes.}

\item{se}{Vector of standard errors.}

\item{eaf}{Vector of effect allele frequencies.}

\item{R}{LD matrix.}

\item{maxsize}{The maximum number of causal variants.}

\item{tau0}{Parameter tau.}

\item{r0}{Parameter $r$.}

\item{niter}{Number of iterations.}

\item{burnin}{Number of burn-in samples.}

\item{p}{Number of variants.}

\item{seed}{Random seed.}

\item{excl.burnin}{Should the burn-in be excluded?}

\item{n}{Sample size (not used).}

\item{a0}{Hyperparameter a for the model size prior.}

\item{b0}{Hyperparameter b for the model size prior.}

\item{inds0}{Initial model indices (not used).}

\item{standardize}{Should the effect sizes be standardised? Defaults to TRUE.}

\item{msprior}{Model size prior (not used).}

\item{verbose}{Verbose output.}

\item{anchor_cslevel}{}
}
\value{
List.
}
\description{
Title
}
