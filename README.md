# TranSpotteR
This package aims to identify annotated and unannotated LINE1 insertion in the genome.

## Overall workflow
### 1. Extracting useful reads (extract_read)
In this step,  informative reads are retrieved for subsequent analysis.
Firstly, the discordant reads are extracted by samtools and the distance between the read pairs should be at least 5000 apart.
Secondly, the reads at the regions of interest can also be included.
For the LINE1 identification case, regions of interest are the three most recently eveolved LINE1 (L1HS, L1PA2 and L1PA3) positions in hg19 that are annotated in RepeatMasker and flanked by 1000bp.
To further reduce the number of reads for the downstream step, only the read pair that contains at least one reads can align to the sequences of interest (in this case, it is Hot_L1_polyA) or reads pair has long clipped reads (potentially polyA) are retrieved.
Combining these steps, the reads for finding reference and non-reference LINE1 insertion can be included.

### 2. Clustering reads (clustering_reads)
The uniquely mapped (MAPQ >= 10) reads are clustered and the cluster should have coverage of at least n reads (default:3 reads) at a point.

### 3. Locally assembling reads clusters (sequence_construction)
For the assembly step, an own assembly function was writtern in R combining C++ via Rcpp due to the lack of a sequence assembly function writtern in R.
The sequence assembly method is Overlaps-Layout-Consensus(OLC).

4. Annotating the constructed reads (annotate_constructed_reads)
5. Inferring the LINE1 integration (line1_inference)
