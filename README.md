# TranSpotteR
This package aims to identify annotated and unannotated LINE1 insertion in the genome.

## Overall workflow
### 1. Extracting useful reads (extract_read)
Reads for inferring insertion are retrieved in this step.
Firstly, the discordant reads are extracted by samtools and the distance between the read pairs should be at least 5000 apart.
Secondly, the reads at the regions of interest can also be supplied.
For the LINE1 identification case, regions of interest are the three most recently evolved LINE1 (L1HS, L1PA2 and L1PA3) positions in hg19 that are annotated in RepeatMasker and those regions are flanked by 1000bp to increase the coverage.
To further reduce the number of reads for the downstream step, only read pairs that contatin  at least one read aligning to the sequences of interest (In the LINE1 identification case, it is a consensus sequence of Hot LINE1) or reads pair has long clipped regions(In this case, it is potentially polyA) are retrieved.
Combining these steps, the reads for finding reference and non-reference LINE1 insertion can be included.

### 2. Clustering reads (clustering_reads)
The uniquely mapped (MAPQ >= 10) reads are clustered to serve as an anchor region.
The cluster should have coverage of at least n reads (default:3 reads) at a point.

### 3. Locally assembling reads clusters (sequence_construction)
For the assembly step, an sequence assembly function was written (in R and C++ via Rcpp) due to the lack of a related function in R.
The sequence assembly employed the Overlaps-Layout-Consensus(OLC) method.
The clustered reads and their partner reads are assembled respectively, then the function will attempt to combine the contigs generated from the clustered reads and the contig from the partner reads to form a long contig.

### 4. Annotating the constructed reads (annotate_constructed_reads)
The contigs are then annotated by aligning to the sequence of interest and the genomic.
Contigs are aligned to sequences of interest (e.g a consensus sequence of Hot LINE1) and the unaligned parts of the read will be collected and proceeded to the next alignment to the genome.
The first alignment is served as a 'bait' to collect the sequence of interest, and the second aligment is to locate the seqeunce in the genome.
The annotations from both alignments are combined and be labelled to the contigs.
For example, ACTCGTGCTTTTCGCTATCGTAGATCGACTAGCA will be annotated as 1:1-10:+ TTCGCTAT 22:90-104:-
In principle, this read annotation function can be entended to many situations rather than just the identification of transposon insertion but also other sequences of intereset,for example, virus genome, integrated plasmid and other ERV.
Apart from annotating the sequence of interest, if a contig maps to two different genomic regions, the function can also annotate them.
Therefore, it can be extended to chromosomal translocation and this will be done in the future.

### 5. Inferring the LINE1 integration (line1_inference)
Under development and will be out soon.

## Example usage
```
extract_read(bam = "tesing.bam", out_dir = getwd())
reads <- import_files(extracted = "tesing_extracted.txt", anchor_min_mapq = 10)
clusters <- clustering_reads(reads)
clusters <- sequence_construction(clusters)
annotation <- annotate_constructed_reads(clusters)
result <- line1_inference(annotation) # under development
```

## In the future
The functions of this package will be extended to detect chromosomal translocation and virus insertion soon.

