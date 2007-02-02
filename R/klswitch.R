klswitch.c <- function(qig,permute,Z.proc,Mu.proc,Z.var)
{
  nnodes <- dim(Z.proc)[1]
  p <- dim(Z.proc)[2]
  nsamp <- dim(Z.proc)[3]

  ngroups <- dim(Mu.proc)[1]
  npermute <- nrow(permute)

  output <- .C("call_klswitch",
               p = as.integer(p),
               nnodes = as.integer(nnodes),
               nsamp = as.integer(nsamp),
               ngroups = as.integer(ngroups),
               npermute = as.integer(npermute),
               Z = as.double(as.vector(Z.proc)),
               mu = as.double(as.vector(Mu.proc)),
               Z.var = as.double(as.vector(Z.var)),
               qig = as.double(as.vector(qig)),
               permute = as.integer(as.vector(permute)),
               minat = as.integer(rep(0,nsamp)),
               PACKAGE="latentnet")
  return(output)
}
