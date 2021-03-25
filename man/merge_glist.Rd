\name{merge_glist}
\alias{merge_glist}

\title{
  Merging each GRanges objects in a GRangesList wihtin a windows.
}
\description{
Each GRanges object in a GRangesList are merged if the GRanges are wihtin the specified distance. Strand information are ignored.
}
\usage{merge_glist(x, tol = 10000)}

\arguments{
  \item{x}{A GRangesList.}
  \item{invert}{An integer specifying the maximum distance between adjacent windows.}
}
\details{}
\value{A GRangeslist containing the merged windows.}

\author{
Cheuk-Ting Law
}

\seealso{}
\examples{
grl <- GRangesList(GRanges(c(1,1), IRanges(c(10, 15), c(18, 20))))
merge_glist(grl)
}
