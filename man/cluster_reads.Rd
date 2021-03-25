\name{cluster_reads}
\alias{cluster_reads}

\title{Clustering anchor reads, sequencing reads passing a threshold of mapping quality (MAPQ).}
\description{
Anchors, reads with sufficiently high mapping quality (e.g. MAPQ > 10) are clustered together.
}
\usage{cluster_reads(total, anchor_min_mapq, min_cov = 5, size_tol = 0, max_reads = 100)}

\arguments{
  \item{total}{A \code{GRanges} object, imported from \code{import_files}, containing all the extracted sequence information.}
  \item{anchor_min_mapq}{A integer specifying the minimum mapping quality score for being a anchor reads.}
  \item{min_cov}{A integer specifying the minimum coverage that a cluster should have.}
  \item{max_reads}{A integer specifying the maxium reads that a cluster can have for reducing the computational cost at the later step, especially, read assembly.}
}
\details{
Here is an example of min_cov=3,
1112221112211111 (min_cov = 2)     1122333222111 (min_cov >= 3)
   --------                        -----------
------             <- filtered   -------        <- included
         -------                     ------
}
\value{
A GRange object with a \linkS4class{DataFrame} that contains the reads consisting the cluster and their partner reads.
}

\author{Cheuk-Ting Law}

\examples{}
