\name{construct_contigs}
\alias{construct_contigs}
\title{Construct contigs from read clusters}
\description{Perform Overlap-Layout-Consensus (OLC) to combine reads from cluster into a contigs.}

\usage{
construct_contigs(x, BPPARAM = MulticoreParam(workers = 10))
}

\arguments{
  \item{x}{A \link[GenomicRanges]{GRanges} object generated from \link{cluster_reads}.}
  \item{BPPARAM}{A \link[BiocParallel]{BiocParallelParam} object controlling parameters in parallelization.}
}
\details{
In the first round of assembly, reads consisting the clusters and their partner reads are assembled separately.
In the second round of assembly,  the cluster contigs and their partners contigs are assembled to create a longer contig.
In every assembly step, more than one contigs can be constructed.
If no longer contig is formed, the cluster contigs and partner contigs are remained. 
On the other hand, if a longer contig is formed, the cluster and partner contig are erased even the contigs do not involve in the composition of the longer contig.
The number of reads used in constructing each contig are also recorded.
}
\value{
A GRange object with a \link[S4class]{DataFrame} that contains reads consisting the cluster, their partner reads, contigs assembled from cluster reads (cluster_contigs), contigs assembled from partner reads (partner_contigs) and contig combining cluster and partner contig (long_contigs).
The number of reads used in constructing each contig are recorded in the \code{elementMetadata} of the contigs.
}

\author{Cheuk-Ting Law}

\seealso{
\code{\link{construct_contigs}}
}
\examples{}

