# `latentnet`: Latent Position and Cluster Models for Statistical Networks

[![Build Status](https://travis-ci.org/statnet/latentnet.svg?branch=master)](https://travis-ci.org/statnet/latentnet)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/latentnet?color=2ED968)](https://cranlogs.r-pkg.org/)
[![cran version](https://www.r-pkg.org/badges/version/latentnet)](https://cran.r-project.org/package=latentnet)
[![R build status](https://github.com/statnet/latentnet/workflows/R-CMD-check/badge.svg)](https://github.com/statnet/latentnet/actions)

Fit and simulate latent position and cluster models for statistical networks. See Krivitsky and Handcock (2008) <10.18637/jss.v024.i05> and Krivitsky, Handcock, Raftery, and Hoff (2009) <10.1016/j.socnet.2009.04.001>.

## Public and Private repositories

To facilitate open development of the package while giving the core developers an opportunity to publish on their developments before opening them up for general use, this project comprises two repositories:
* A public repository `statnet/latentnet`
* A private repository `statnet/latentnet-private`

The intention is that all developments in `statnet/latentnet-private` will eventually make their way into `statnet/latentnet` and onto CRAN.

Developers and Contributing Users to the Statnet Project should read https://statnet.github.io/private/ for information about the relationship between the public and the private repository and the workflows involved.

## Latest Windows and MacOS binaries

A set of binaries is built after every commit to the repository. We strongly encourage testing against them before filing a bug report, as they may contain fixes that have not yet been sent to CRAN. They can be downloaded through the following links:

* [MacOS binary (a `.tgz` file in a `.zip` file)](https://nightly.link/statnet/latentnet/workflows/R-CMD-check.yaml/master/macOS-rrelease-binaries.zip)
* [Windows binary (a `.zip` file in a `.zip` file)](https://nightly.link/statnet/latentnet/workflows/R-CMD-check.yaml/master/Windows-rrelease-binaries.zip)

You will need to extract the MacOS `.tgz` or the Windows `.zip` file from the outer `.zip` file before installing. These binaries are usually built under the latest version of R and their operating system and may not work under other versions.
