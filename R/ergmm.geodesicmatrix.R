ergmm.geodesicmatrix<-function(g, directed=is.directed(g)){
 ergmm.geodesicmatrix.edgelist(edgelist=as.matrix.network(g,matrix.type="edgelist"),
               n=g$gal$n, directed=directed)
}

ergmm.geodesicmatrix.edgelist <- function(edgelist, n=max(edgelist), directed=FALSE) {
# edgelist is an mx2 matrix of edges.  n=#nodes.
# This function returns an nxn matrix, where
#      M[i,j] : length of shortest path from vertex i to vertex j
  
# This function starts off just like ergmm.geodistn:
  edgelist<-edgelist[edgelist[,1]!=edgelist[,2],] # get rid of self-edges
  if (!directed) 
    edgelist<-rbind(edgelist,edgelist[,2:1])
  edgelist<-unique(edgelist)
  edgelist<-edgelist[order(edgelist[,1],edgelist[,2]),]
  nodelist<-match(1:n,edgelist[,1],nomatch=1)-1
  
# Now everything is ready.  Call the C code.
  ans<-.C("geodesic_matrix", as.integer(t(edgelist)), as.integer(n),
    as.integer(nodelist), as.integer(dim(edgelist)[1]), colors=integer(n),
    gmat=integer(n*n), queue=integer(n), PACKAGE="latentnet") $ gmat
  ans[ans==n]<-Inf # length n really means no path exists
  ans=matrix(ans,n,n,byrow=TRUE) # byrow=TRUE is only important when directed==TRUE
  ans
}
