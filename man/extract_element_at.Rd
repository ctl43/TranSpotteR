\name{extract_element_at}
\alias{first_element}
\alias{last_element}
\alias{extract_element_at}

\title{
  Subsetting element from a list-like type object.
}
\description{
Extracting the element from each slot of a list. It behaves like \code{lapply}, e.g. \code{lapply(x, tail, n = 1} is similar to \code{last_element(x)}, but the latter one will be much faster.
For the speed, \code{extract_element_at(x, 1)} is slightly slower than \code{first_element}.
}
\usage{
first_element(x, invert = FALSE)
last_element(x, invert = FALSE)
extract_element_at(x, invert = FALSE)
}

\arguments{
  \item{x}{A list like object.}
  \item{invert}{A logical scalar specifying whether extract element except the selected one, \code{x}.}
}
\details{It also works with \linkS4class{GRangesList} and \linkS4class{List}.}
\value{A list containing the selected element from the each slot in the list.}

\author{Cheuk-Ting Law}

\examples{
x <- list(c(1:10), c(5:20))
first_element(x)
first_element(x, invert = TRUE)
last_element(x)
extract_element_at(x, 3)
}
