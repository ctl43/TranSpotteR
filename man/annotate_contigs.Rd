\name{annotate_contigs}
\alias{annotate_contigs}

\title{
Contig annotation
}
\description{
Contigs generated from \link{construct_contigs} are annotated to aligned regions.
}
\usage{
annotate_contigs(x,
                 insert,
                 genome,
                 customised_annotation = customised_annotation,
                 BPPARAM = MulticoreParam(workers = 3L))

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A GRanges object generated from construct_contigs.}
  \item{insert}{A file path pointing to a BWA indexed fasta sequence. The sequence are the targeted sequence intended to bed identified in the genome,e.g. LINE1, virus genome and plasmid sequence, etc.}
  \item{genome}{A file path pointing to BWA indexed genome sequence.}
  \item{customised_annotation}{A list object containing functions that further annotate the sequence, e.g. identifying the polyA sequence, etc.}
  \item{BPPARAM}{A \link[BiocParallel]{BiocParallelParam} object controlling parameters in parallelization.
  When the number of workers is larger than 3, it will not improve the performance because it only parallelizes the annotation of cluster contigs, partner contigs and long contigs.}
}
\details{
The annotation process includes the following steps,
1. Contigs are aligned to the insert sequence first and annotate the primary alignment (regardless the mapping quality) to the sequence.
2. The unaligned parts are extracted and aligned to the genome.
Only aligned regions with mapping quality greater than 10 are annotated to the contigs.
If supplementary alignments are found, the primary aligmnet will be assigned to the sequence first, then the supplementary alignment. The assignmnet of supplementary alignment do not follow a specific rule, it follows the order of appearance.
3. The unannotated reads are subjected to customised annotations to identify customised structures, for example, polyA sequence.
4. The annotation from cluster contigs and partners contigs are joined together, while the annotation from long contigs will be remained the same.
}
\value{
A \link[data.table]{data.table} object with fields, contig_detail, nreads, cluster_origin and cluster_region.
\item{contig_detail}{A list of data.table contains annotations of contigs, see \link{annotate_seq}.}
\item{nreads}{The number of reads consisting the contig.}
\item{cluster_origin}{Index of cluster that the annotation originates.}
\item{cluster_region}{Genomic region of the cluster.}
}

\author{Cheuk-Ting Law}

\seealso{
\link{annotate_seq}
}
\examples{}
