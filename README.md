# TranSpotteR
This package aims to identify annotated and unannotated LINE1 insertion in the genome.

## Overall workflow
### 1. Extracting useful reads (extract_read)
In this step,  informative reads are retrieved for subsequent analysis.
Firstly, the discordant reads are extracted by samtools and the distance between the read pairs should be at least 5000 apart.
Secondly, the reads at the regions of interest can also be included.
For the LINE1 identification case, regions of interest are the three most recently evolved LINE1 (L1HS, L1PA2 and L1PA3) positions in hg19 that are annotated in RepeatMasker and flanked by 1000bp.
To further reduce the number of reads for the downstream step, only the read pair that contains at least one reads can align to the sequences of interest (in this case, it is Hot_L1_polyA) or reads pair has long clipped regions (potentially polyA) are retrieved.
Combining these steps, the reads for finding reference and non-reference LINE1 insertion can be included.

### 2. Clustering reads (clustering_reads)
The uniquely mapped (MAPQ >= 10) reads are clustered and the cluster should have coverage of at least n reads (default:3 reads) at a point.

### 3. Locally assembling reads clusters (sequence_construction)
For the assembly step, an own assembly function was written (in R and C++ via Rcpp) due to the lack of a sequence assembly function written in R.
The sequence assembly employed the Overlaps-Layout-Consensus(OLC) method.
It assembles the clustered reads and their partner reads, then it attempts to assemble the contig generated from the clustered reads and the contig from the partner reads.

### 4. Annotating the constructed reads (annotate_constructed_reads)
The contigs are then annotated by aligning to the sequence of interest and the genomic.
Contigs are aligned to the sequence of interest (Hot_L1_polyA, in this case) and parts of the read that are not aligned will be collected and proceed to the next alignment to the genome.
The annotations from both alignments are combined and be labelled to the contigs.
For example, ACTCGTGCTTTTCGCTATCGTAGATCGACTAGCA will be annotated as 1:1-10:+ TTCGCTAT 22:90-104:-
In principle, this read annotation function can apply to many situations.
It does not only limit to annotate LINE1, it can also annotate other sequences of interest, for example, virus genome, integrated plasmid and other ERV.
Apart from  annotating the sequence of interest, if parts of a read map to two different genomic regions, the function can also annotate them.
Therefore, it can be extended to chromosomal translocation, which will be done in the future.

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
The functions of this package will be extended to detecting chromosomal translocation and virus insertion soon.

