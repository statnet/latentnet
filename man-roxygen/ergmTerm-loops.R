#  File man-roxygen/ergmTerm-loops.R in package latentnet, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#' @description It is recommended to only use this term if your network contains self-loops. Otherwise, use the \CRANpkg{ergm} term(s) \Sexpr[results=rd,stage=build]{statnet.common::paste.and(paste0("\\\\ergmTerm{ergm}{", strsplit("<%= ergm_analogue %>", " ")\[\[1\]\], "}{()}"), con = "or")}.
