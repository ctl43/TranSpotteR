# TranSpotteR
This package aims to identify annotated and unannotated LINE1 insertion in the genome.

## Overall workflow
### 1. Extracting useful reads (extract_reads)
Reads for inferring insertion are retrieved in this step.
Firstly, the discordant reads are extracted by samtools and the distance between the read pairs should be at least 5000 apart.
Secondly, the reads at the regions of interest can also be supplied.
For the LINE1 identification case, regions of interest are the three most recently evolved LINE1 (L1HS, L1PA2 and L1PA3) positions in hg19 that are annotated in RepeatMasker and those regions are flanked by 1000bp to increase the coverage.
To further reduce the number of reads for the downstream step, only read pairs that contatin  at least one read aligning to the sequences of interest (In the LINE1 identification case, it is a consensus sequence of Hot LINE1) or reads pair has long clipped regions(In this case, it is potentially polyA) are retrieved.
Combining these steps, the reads for finding reference and non-reference LINE1 insertion can be included.

### 2. Clustering reads (cluster_reads)
Reads with MAPQ >= 10 are considered as uniquely mapped and clustered together.

### 3. Assembling reads clusters (construct_contigs)
For the assembly step, an sequence assembly function was written (in R and C++ via Rcpp) due to the lack of a related function in R.
The sequence assembly employed the Overlaps-Layout-Consensus(OLC) method.
The clustered reads and their partner reads are assembled respectively, then the function will attempt to combine the cluster contig and the partner contig to form a long contig.

### 4. Annotating the constructed reads (annotate_contigs)
By aligning to the sequence of interest and the genomic, reads are annotated to the corresponding mapped regions.
Firetly, Contigs are aligned to sequences of interest (e.g a consensus sequence of Hot LINE1) and the unaligned parts of the read will be subjected to the next alignment aginst the genome.
The first alignment is served as a 'bait' to collect the sequence of interest, and the second aligment is to locate the insertion site in the genome.
Besides, the polyA sequence will also be identified and annotated in this step..

### 5. Inferring the LINE1 integration (infer_tranposon)
Under development and will be out soon.

## Example usage
```
extract_reads(bam = "tesing.bam", out_dir = getwd())
reads <- import_files(extracted = "tesing_extracted.txt", anchor_min_mapq = 10)
clusters <- cluster_reads(reads)
clusters <- construct_contigs(clusters)
annotation <- annotate_contigs(clusters, insert = "LINE1.fa", genome="hg19.fa")
result <- infer_tranposon(annotation) # under development
```

## In the future
The applications of this package will be extended to detect chromosomal translocation and virus insertion soon.
