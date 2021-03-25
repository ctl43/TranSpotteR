\name{extract_reads}
\alias{extract_reads}
\title{Extract informative reads for insertion inference}
\description{Extract discordant reads and reads mapped to sequence in interest.}

\usage{
extract_reads(bam, tmp_dir = NULL,
              out_dir = getwd(), chromosome = c(1:22, "X", "Y", "KJ173426"),
              samtools = "samtools",
              interested_region = NULL,
              interested_seq = NULL,
              keep_intermediate = FALSE,
              threads = 5,
              high_mapq = 5)
}

\arguments{
  \item{bam}{A file path pointing to a bam file.}
  \item{tmp_dir}{A file path pointing to a folder for storage temporary files.
  If it is \code{NULL}, '/tmp' is used.}
  \item{out_dir}{A directory to store the output. By default, it is the working directory.}
  \item{chromosome}{Chromosome/RNAME of interest. Therefore, those supplementary chromosomes are filtered.}
  \item{samtools}{A path to call samtools if it cannot be called simply by samtools in the system}
  \item{interested_region}{A file path to a bed file containing the regions of interest.}
  \item{interested_seq}{A file path to a fasta file containing the sequence of interest, e.g. transposon sequence or virus genome.}
  \item{keep_intermediate}{A logical scalar specifying whether it keep the intermediate files.}
  \item{threads}{An integer specifying the number of CPU threads used in bwa alignment. If \code{threads} => 2, the reads in regions of interest and disordant reads are extracted in parallel by samtools.}
  \item{high_mapq}{An integer specifying the mapping quality (MAPQ) threshold. At least one reads should pass this threshold in a read pair.}
}
\details{
  The discordant read pairs and reads in regions of interest are extracted by samtools.
  Here are some filters to reduce the size of extracted reads.
  1. At least on reads mapped to the chromosomes/RNAME of interest.
  2. only reads pairs that are at least 5000bp apart are considered as discordant reads.
  3. To further narrow down the number of reads, if both reads of a read pair perfectly map to the sequence of interests, the read pairs will not be used because they cannot provide useful information to locate the potential insertion site.
  However, if both reads are mapped to the sequence of interests but with long clipped sequence (>20bp), they are included for the further analysis.
}
\value{
A text file containing extracted reads.
The text file is basically simplified sam format, which contains the fields of QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, SEQUENCE and QNAME_id.
}

\author{
Cheuk-Ting Law
}

\examples{
extract_reads("my.bam")
}

