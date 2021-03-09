# TranSpotteR
This package aims to identify annotated and unannotated LINE1 insertion in the genome.

## Overall workflow
### 1. Extracting useful reads (extract_read)
In this step, the informative reads are retrieved for subsequent analysis.
Firstly, the discordant reads are extracted and the distance between the read pairs should be at least 5000 apart.
Secondly, the reads at the interested regions are also retrieved.
For the LINE1 identification case, the interested regions are the three most recently eveolved LINE1 (L1HS, L1PA2 and L1PA3) positions in hg19, that are annotated in RepeatMasker and flanked by 1000bp.
Combining two steps, the reference LINE1 insertion and the non-reference LINE1 insertion can be included.
To further reduce the number of reads for the downstream step, only the read pair that contains at least one reads can align to the interested sequence (in this case, it is Hot_L1_polyA) or reads pair has long clipped reads (potentially polyA) are retrieved.
The collected reads are combined and exported for the next step.

2. Clustering reads (clustering_reads)
3. De novo assembly of the reads clusters (sequence_construction)
4. Annotating the constructed reads (annotate_constructed_reads)
5. Inferring the LINE1 integration (line1_inference)
