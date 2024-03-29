# TranSpotteR
This package aims to identify annotated and unannotated LINE1 insertion in the genome.

## Workflow
### 1. Extracting useful reads (`extract_reads`)
Reads that are potentially informative for inferring insertion are retrieved, including discordant reads, reads aligning to the sequence of intereset and reads at the region of interest.

### 2. Clustering reads (`cluster_reads`)
Reads with MAPQ >= 10 are considered as uniquely mapped and clustered together.

### 3. Assembling reads clusters (`construct_contigs`)
For the assembly step, a self-writtern assembly function employing the Overlaps-Layout-Consensus(OLC) method is applied to assemble reads in read cluster into longer contigs.

### 4. Annotating the constructed reads (`annotate_contigs`)
Contigs are annotated to the regions that they align to.
First, contigs are aligned to sequences of interest (e.g a consensus sequence of Hot LINE1).
Then, the unaligned parts of the read will be subjected to the next alignment aginst the genome.
The first alignment is served as a 'bait' to collect all the sequence of interest, and the second aligment is to locate the insertion site in the genome.
In addition, some useful feature, for example, polyA sequence will also be identified and annotated in this step.

### 5. Inferring the LINE1 integration (`infer_transposon`)
Under development and will be out soon.

### 6. Inferring viral/plasmid integration (`infer_simple_insertion`)

## Example usage (LINE1 insertion)
```r
extract_reads(bam = "tesing.bam", out_dir = getwd())
reads <- import_files(extracted = "tesing_extracted.txt", anchor_min_mapq = 10)
clusters <- cluster_reads(reads)
clusters <- construct_contigs(clusters) # The most time-consuming sttep;saving the result is recommended
saveRDS(clusters, "test_contigs.rds")
annotation <- annotate_contigs(clusters, insert = "LINE1.fa", genome="hg19.fa")
result <- infer_transposon(annotation) # under development (for LINE1)
```

## Example usage (viral insertion)
```r
extract_reads(bam = "tesing.bam", out_dir = getwd())
reads <- import_files(extracted = "tesing_extracted.txt", anchor_min_mapq = 10)
clusters <- cluster_reads(reads)
clusters <- construct_contigs(clusters)
# Or simply import the result generated for LINE1 insertion
# clusters <- readRDS("test_contigs.rds")
annotation <- annotate_contigs(clusters, insert = "HBV.fa", genome="hg19.fa")
result <- infer_simple_insertion(annotation)
```

## In the future/To-do list
1. The applications of this package will be extended to detect chromosomal translocation and transduction of viral genome.
2. Long-read sequencing reads will be accepted later.
3. Trace back the LINE1 origin from the genome based on sequence similarity (although it is likely that many of them are untraceable).
