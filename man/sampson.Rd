\name{sampson}
\alias{sampson}
\alias{samplike}
\docType{data}
\title{Cumulative network of positive affection within a monastery as a ``network'' object}
\description{
Sampson (1969) recorded
the social interactions among a group of monks
while resident as an experimenter on vision, and collected
numerous sociometric rankings.  
During his stay, a political
``crisis in the cloister'' resulted in the expulsion of four monks
(Nos. 2, 3, 17, and 18) and the voluntary departure of several
others - most immediately, Nos. 1, 7, 14, 15, and 16. (In the
end, only 5, 6, 9, and 11 remained).
Of particular interest
is the data on positive affect relations (``liking''), in which each
monk was asked if they had positive relations to each of the
other monks.

The data were gathered at three times to capture changes in group sentiment
over time.
They were represent three time points in the period during which a
new cohort entered the monastery near the end of the study but
before the major conflict began. 

Each member ranked only his top three choices on ``liking''.
(Some subjects offered tied ranks for their top four choices).
A tie from monk A to monk B exists if A
nominated B as one of his three best friends at that that time point.

\code{samplike} is the time-aggregated network. 
It is the cumulative tie for ``liking'' over the three
periods.  For this, a tie from monk A to monk B exists if A
nominated B as one of his three best friends at any of the
three time points.

This data is standard in the social network analysis
literature, having been modeled by Holland and Leinhardt
(1981), Reitz (1982), Holland, Laskey and Leinhardt (1983),
and Fienberg, Meyer, and Wasserman (1981), Hoff, Raftery, and Handcock
(2002), etc.
This is only a small piece of the data collected by Sampson.

For additional examples of \code{network}s see the \code{netdata}
package which is part of the \code{statnet} package.
}
\usage{
 data(sampson)
}
\source{
 Sampson, S.~F. (1968), \emph{A novitiate in a period of change: An
 experimental and case study of relationships,} 
 Unpublished ph.d. dissertation, Department of Sociology, Cornell University.
}
\seealso{network, \code{\link[network]{plot.network}}, \code{\link{ergmm}}}
\references{
 White, H.C., Boorman, S.A. and Breiger, R.L. (1976).
 \emph{Social structure from multiple networks. I. Blockmodels of roles and
 positions.}
 American Journal of Sociology, 81(4), 730-780.
}
\examples{
\dontrun{
data(sampson)
plot(samplike)
}
}
\keyword{datasets}
\keyword{cluster}
\keyword{graphs}