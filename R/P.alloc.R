#  File R/P.alloc.R in package latentnet, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
## Lets the user call P.free.all from R. Only used for debugging.
P.free.all<-function() invisible(try(.C("P_free_all", PACKAGE="latentnet")))
