\name{import_files}
\alias{import_files}
\title{Extract informative reads for insertion inference}
\description{Extract discordant reads and reads mapped to sequence in interest.}

\usage{
import_files <- function(extracted, seq_info)
}

\arguments{
  \item{extracted}{Paths to all reads extracted from \code{extract_reads}. The number of paths can be more than 1.}
  \item{seqinfo}{A \link[GenomeInfoDb]{Seqinfo} object containing the chromosome information. It can be converted from a chrom_sizes.txt file, retrieved from sam/bam header or retrieved from R packages directly, e.g. BSgenome.}
}
\details{
All extract data are imported into the R system.
}
\value{
A text file containing extracted reads.
The text file is basically simplified sam format, which contains the fields of QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, SEQUENCE and QNAME_id.
}

\author{
Cheuk-Ting Law
}

\examples{
# For hs37d5,
library(BSgenome.Hsapiens.1000genomes.hs37d5)
seq_info <- seqinfo(hs37d5)
import_files("my_extract.txt", seq_info)
}

