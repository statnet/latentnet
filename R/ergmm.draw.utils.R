#  File R/ergmm.draw.utils.R in package latentnet, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
#######################################################################
ergmm.drawcircle <- function(center,radius,length=50,...)
{
  x0 <- seq(-radius,radius,length=length)
  x1 <- seq(radius,-radius,length=length)
  x <- c(x0,x1)
  y <- c(sqrt(radius^2 - x0^2),-sqrt(radius^2 - x1^2))
  #' @importFrom graphics lines
  lines(x+center[1],y+center[2],...)
}




#' Draw a pie chart at a specified location.
#' 
#' Used by \code{\link{plot.ergmm}} to draw pie charts to visualize soft
#' clusterings when \code{pie=TRUE}. Exported as a courtesy to dependent
#' packages.
#' 
#' 
#' @param center A numeric vector of length 2, specifying the horizontal and
#' the vertical coordinates of its center.
#' @param radius Radius of the pie chart.
#' @param probs A vector of probabilities/weights of each sector; they do not
#' have to sum to 1.
#' @param n Number of points to use to approximate the "circle".
#' @param cols A vector of colors to use for the sectors.
#' @param \dots Additional arguments, currently unused.
#' @author See COPYRIGHT.
#' @seealso plot.ergmm
#' @keywords graphs
#' @examples
#' 
#' plot(c(0,sum(1:11))*2,c(-10,10),type="n",asp=1)
#' for(i in 1:10) ergmm.drawpie(c(sum(1:i)*2,0), radius=i, probs=1:(i+1))
#' 
#' @export
ergmm.drawpie <- function(center,radius,probs,n=50,cols=1:length(probs),...)
{
  x <- c(0,cumsum(probs)/sum(probs))
  dx <- diff(x)
  np <- length(probs)
  for (i in 1:np)
  {
    t2p <- 2 * pi * seq(x[i], x[i + 1], length = n)
    xc <- center[1] + c(cos(t2p), 0) * radius
    yc <- center[2] + c(sin(t2p), 0) * radius
    #' @importFrom graphics polygon
    polygon(xc, yc, border = FALSE, col = cols[i])
  }
  ergmm.drawcircle(center=center,radius=radius,col=1)
}
