\name{assemble_reads}
\alias{assemble_reads}
\title{Assemble sequences into a contig}
\description{Perform Overlap-Layout-Consensus (OLC) to combine reads into a contig.}

\usage{
assemble_reads(vec, n_reads = NULL, msa_result = FALSE, 
    add_id = TRUE, min_len = 8L, min_pid = 85, consensus_min_len = 100)
}

\arguments{
  \item{vec}{A vector containing the sequences to be assembled.}
  \item{n_reads}{A integer vector to record the number of reads that are used to construct the sequences in \code{vec}. It can be ignored for general usage.}
  \item{msa_result}{A logical scalar specifying whether the multiple sequence alignment result is returned.}
  \item{add_id}{A logical scalar specifying whether Universally unique identifiers (UUID) are used to name the returned contigs. If it is \code{FALSE}, the names of the returned contigs will be the sequence indexs that are used in the assembly.}
  \item{min_len}{An integer specifying the minimum overlapping length in the step of overlapping reads.}
  \item{min_pid}{A numeric object specifying the minimum percentage of identity (PID) between the overlapping region.}
  \item{consensus_min_len}{An integer specifying the minimum length of the returned contigs.}
}
\details{
By default, the function does not return the reads that have no overlapping with others.
}
\value{
A list that contains
\item{consensus}{A string vector containing assembled contigs}
\item{n_reads}{An integer vector indicating the number of reads used to construct the contigs}
If \code{msa_result} is TRUE,
the list will have an extra slot.
\item{msa}{A character vector containing a view of MSA.}
}

\author{
Cheuk-Ting Law
}

\seealso{
\code{\link{construct_contigs}}
}
\examples{
seq <- c("ACTACGATCGTAGCTGCTGTGTATGCAT", "GCTGCTGTGTATGCATGCTGACGTCGT")
# ACTACGATCGTAGCTGCTGTGTATGCAT-----------
# ------------GCTGCTGTGTATGCATGCTGACGTCGT
assemble_reads(seq, consensus_min_len = 10, add_id = FALSE)
}

