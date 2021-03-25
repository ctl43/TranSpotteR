\name{annotate_seq}
\alias{annotate_seq}

\title{
Annotate sequences by aligning to different reference sequences sequentially.
}
\description{
Annotate sequences to the aligned regions.
}
\usage{
annotate_seq(seq, insert, genome, 
             customised_annotation = customised_annotation,
             BPPARAM = MulticoreParam(workers = 3L))

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seq}{A vector of strings with names.}
  \item{insert}{A file path pointing to a BWA indexed fasta sequence. The sequence are the targeted sequence intended to bed identified in the genome,e.g. LINE1, virus genome and plasmid sequence, etc.}
  \item{genome}{A file path pointing to BWA indexed genome sequence.}
  \item{customised_annotation}{A list object containing functions that further annotate the sequence, e.g. identifying the polyA sequence, etc.}
  \item{BPPARAM}{A \link[BiocParallel]{BiocParallelParam} object controlling parameters in parallelization.
  When the number of workers is larger than 3, it will not improve the performance.}
}
\details{
The annotation process includes the following steps,
1. Contigs are aligned to the insert sequence first and annotate the primary alignment (regardless the mapping quality) to the sequence.
2. The unaligned parts are extracted and aligned to the genome.
Only aligned regions with mapping quality greater than 10 are annotated to the contigs.
If supplementary alignments are found, the primary aligmnet will be assigned to the sequence first, then the supplementary alignment. The assignmnet of supplementary alignment do not follow a specific rule, it follows the order of appearance.
3. The unannotated reads are subjected to customised annotations to identify customised structures, for example, polyA sequence.
4. The unannotated parts in reads are remained to be the sequence itself.
}
\value{
A \link[data.table]{data.table} object with fields, start, end, width, QNAME, annotation, cigar and seq.
\item{start}{Starting location of the annotation in the read.}
\item{end}{Ending location of the annotation in the read.}
\item{width}{Length of the annotated part.}
\item{QNAME}{Names of the annotated sequence.}
\item{annotation}{Annotation of the read part, e.g. aligned genomic regions, aligned insertion or other customised annotation, etc.}
\item{cigar}{Concise Idiosyncratic Gapped Alignment Report (CIGAR)}.
\item{seq}{Original sequence of the annotated part.}
}

\author{Cheuk-Ting Law}

\seealso{
\link{construct_contigs}
}
\examples{}