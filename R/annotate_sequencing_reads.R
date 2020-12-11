#' @export
#' @importFrom S4Vectors 'elementMetadata<-' elementMetadata split
#' @importFrom BiocGenerics start end strand
#' @importFrom BiocParallel bplapply MulticoreParam

# A wrapper for read annotation
annotate_constructed_reads <- function(x, BPPARAM = MulticoreParam(workers = 3L))
  # A wrapper function process both plus and minus stranded read clusters.
  # Written by Cheuk-Ting Law
{
  p <- x[strand(x) == "+"]
  m <- x[strand(x) == "-"]
  out <- bplapply(list(p, m), .internal_annotation, BPPARAM = BPPARAM)
  out <- c(out[[1]], out[[2]])
  names(out) <- seq_along(out)
  return(out)
}

##################################
##################################
#' @export
#' @importFrom IRanges setdiff Views
add_back_seq <- function(x, ir = NULL){
  if(!is.null(ir)){
    start <- start(ir)
    end <- end(ir)
    to_keep <- elementMetadata(ir)$annotation
  }
  n_start <- length(start)
  n_end <- length(end)
  if(n_start == 0 | n_end == 0 ){
    return(list(start = 1, end = nchar(x), info = x, seq = x))
  }

  if(n_start != n_end){
    stop("The length of start and end are not the same.")
  }
  idx <- order(start)
  start <- start[idx]
  end <- end[idx]
  to_keep <- to_keep[idx]

  # To deal with overlapping range
  # The one on the left hand side will occupy the place first, then the second one.
  diff <- start[-1] - head(end, -1)
  diff[diff <= 0] <- 1
  start <- c(start[1], diff + head(end, -1))
  need_start <- c(1, end + 1)
  need_end <- c(start - 1, nchar(x))
  selected_seq <- stringr::str_sub(x, start = need_start, end = need_end)
  need_start <- need_start[selected_seq!=""]
  need_end <- need_end[selected_seq!=""]
  selected_seq <- selected_seq[selected_seq!=""]

  start <- c(start, need_start)
  end <- c(end, need_end)
  to_keep <- c(to_keep, selected_seq)
  idx <- order(start)
  start <- start[idx]
  end <- end[idx]
  to_keep <- to_keep[idx]
  seq <- stringr::str_sub(x, start = start, end = end)
  list(start = start, end = end, info = to_keep, seq = seq)
}

#' @export
#' @importFrom IRanges CharacterList
.internal_annotation <-  function(clusters,  BPPARAM = MulticoreParam(workers = 3))
  # A wrapper function to annotate clustered reads, the partner reads from the read cluster ,
  # and long contigs that consist of clustered reads and their partner reads.
  # Written by Cheuk-Ting Law
{
  strand <- as.character(unique(strand(clusters)))
  if(length(strand) > 1){
    stop("")
  }

  if(length(strand) == 0){
    return(NULL)
  }

  # Can be sped up by combining all the reads together and only doing once.
  cluster_anno <- .annotate_reads(unlist(clusters$cluster_contigs), BPPARAM = BPPARAM)
  partner_anno <- .annotate_reads(unlist(clusters$partner_contigs), BPPARAM = BPPARAM)
  long_anno <- .annotate_reads(unlist(clusters$long_contigs), BPPARAM = BPPARAM)
  cluster_grp <- factor(rep(seq_along(clusters$cluster_contigs), lengths(clusters$cluster_contigs)), levels = seq_along(clusters))
  partner_grp <- factor(rep(seq_along(clusters$partner_contigs), lengths(clusters$partner_contigs)), levels = seq_along(clusters))
  long_grp <- factor(rep(seq_along(clusters$long_contigs), lengths(clusters$long_contigs)), levels = seq_along(clusters))

  # Combining annotation
  cluster_anno <- lapply(cluster_anno, split, f = cluster_anno$grp)
  cluster_anno <- lapply(cluster_anno, split, f = cluster_grp)
  partner_anno <- lapply(partner_anno, split, f = partner_anno$grp)
  partner_anno <- lapply(partner_anno, split, f = partner_grp)
  long_anno <- lapply(long_anno, split, f = long_anno$grp)
  long_anno <- lapply(long_anno, split, f = long_grp)

  tmp_fun <- function(x, y){mapply(c, rep(x, each = length(y)), rep(y, length(x)), SIMPLIFY = FALSE)}
  cluster_anno$grp <- lapply(cluster_anno$grp, function(z)lapply(z, as.character))
  partner_anno$grp <- lapply(partner_anno$grp, function(z)lapply(z, as.character))
  long_anno$grp <- lapply(long_anno$grp, function(z)lapply(z, as.character))

  if(strand=="+"){ # Creating all possible combination of read structures
    combined_start <- mapply(tmp_fun, x = cluster_anno$start, y = partner_anno$start, SIMPLIFY = FALSE)
    combined_end <- mapply(tmp_fun, x = cluster_anno$end, y = partner_anno$end, SIMPLIFY = FALSE)
    combined_anno <- mapply(tmp_fun, x = cluster_anno$annotation, y = partner_anno$annotation, SIMPLIFY = FALSE)
    combined_seq <- mapply(tmp_fun, x = cluster_anno$seq, y = partner_anno$seq, SIMPLIFY = FALSE)
    combined_grp <- mapply(tmp_fun, x = cluster_anno$grp, y = partner_anno$grp, SIMPLIFY = FALSE)
    combined_nreads <- IntegerList(mapply(function(x, y)rep(x, each = length(y)) + y,
                                           x = elementMetadata(clusters$cluster_contigs)$n_reads,
                                           y = elementMetadata(clusters$partner_contigs)$n_reads,
                                           SIMPLIFY = FALSE))
  }else{
    combined_start <- mapply(tmp_fun, x = partner_anno$start, y = cluster_anno$start, SIMPLIFY = FALSE)
    combined_end <- mapply(tmp_fun, x = partner_anno$end, y = cluster_anno$end, SIMPLIFY = FALSE)
    combined_anno <- mapply(tmp_fun, x = partner_anno$annotation, y = cluster_anno$annotation, SIMPLIFY = FALSE)
    combined_seq <- mapply(tmp_fun, x = partner_anno$seq, y = cluster_anno$seq, SIMPLIFY = FALSE)
    combined_grp <- mapply(tmp_fun, x = partner_anno$grp, y = cluster_anno$grp, SIMPLIFY = FALSE)
    combined_nreads <- IntegerList(mapply(function(x, y)rep(x, each = length(y)) + y,
                                          x = elementMetadata(clusters$partner_contigs)$n_reads,
                                          y = elementMetadata(clusters$cluster_contigs)$n_reads,
                                          SIMPLIFY = FALSE))
  }
  tmp_fun_2 <- function(i, j, k, d, g){
    mapply(data.table, start = i, end = j, anno = k, QNAME = d, seq = g, SIMPLIFY = FALSE)
  }
  combined_info <- mapply(tmp_fun_2, i = combined_start, j = combined_end, k = combined_anno, d = combined_grp, g = combined_seq, SIMPLIFY = FALSE)
  long_info <- mapply(tmp_fun_2, i = long_anno$start, j = long_anno$end, k = long_anno$annotation, d = long_anno$grp, g = long_anno$seq, SIMPLIFY = FALSE)
  long_nreads <- elementMetadata(clusters$long_contigs)$n_reads
  combined_nreads[lengths(long_nreads) > 0] <- long_nreads[lengths(long_nreads) > 0]
  combined_info[lengths(long_nreads) > 0] <- long_info[lengths(long_nreads) > 0]
  elementMetadata(clusters)$read_annotation <- combined_info
  elementMetadata(clusters)$nreads <- combined_nreads
  table(lengths(clusters$read_annotation) == lengths(combined_nreads))
  clusters
}

#' @importFrom  GenomicAlignments cigarToRleList
cigar_convert <- function(cigar_string, from, to)
{
  rle_list <-  cigarToRleList(cigar_string)
  grp <- rep(seq_along(rle_list), lengths(rle_list))
  rle_list <- unlist(rle_list)
  rle_list[rle_list == from] <- to
  rle_list <- split(rle_list, grp)

  # Converting to cigar strings
  rl <- runLength(rle_list)
  rv <- runValue(rle_list)
  grp <- rep(seq_along(rl), lengths(rl))
  as.character(lapply(split(paste0(unlist(rl),unlist(rv)), grp), paste, collapse = ""))
}

#' @export
#' @importFrom GenomicAlignments cigarWidthAlongQuerySpace cigarRangesAlongQuerySpace
#' @importFrom S4Vectors split runValue runLength nchar
#' @importFrom IRanges RleList CharacterList

.annotate_reads <- function(seq, ref = "~/dicky/reference/fasta/line1_reference/hot_L1_polyA.fa",
                            BPPARAM = MulticoreParam(workers = 3))
  # This function annotates chimeric sequences that consist of any multi-mapped sequences/non-reference sequence and any genomic regions.
  # It first mapped the bait (repeated regions/non-reference sequence), then the genomic regions by bwa
  # Written by Cheuk-Ting Law
{
  if(length(seq) == 0){
    return(character())
  }

  if(any(duplicated(names(seq)))){
    stop("Sequence names are duplicated.")
  }

  if(is.null(names(seq))){
    names(seq) <- paste0("read_", seq_along(seq))
  }

  # Mapping to LINE1 sequence
  aln_1 <- bwa_alignment(seq, ref = ref, samtools_param = "")
  aln_1 <- convertingHtoS(aln_1, unique_id = "QNAME")

  # Annotating the reads
  aln_1$unified_cigar <- unify_cigar_strand(aln_1$CIGAR, flag = aln_1$FLAG, to = "+", along_query = TRUE)

  # Getting clipped sequences
  clipped_seq_1 <- get_unmapped_clipped_read(aln_1, include_middle_unmapped = TRUE) # need debug

  # Mapping to genome
  middle <- clipped_seq_1$middle
  clipped_seq_1 <- clipped_seq_1$clipped

  # Can be sped up by combining all the reads together and only doing once.
  aln_2 <- bplapply(clipped_seq_1, bwa_alignment, call_bwa = "bwa mem ",
                    samtools_param = "-F 128 -F 4", BPPARAM = BPPARAM)
  aln_2 <- lapply(aln_2, convertingHtoS, unique_id = "QNAME")
  aln_2 <- lapply(aln_2, function(x)x[x$MAPQ > 10,])
  aln_2 <- lapply(aln_2, function(x){x$unified_cigar <- unify_cigar_strand(x$CIGAR, flag = x$FLAG, to = "+", along_query = TRUE); x})

  mapping_2 <- lapply(aln_2, sam2gr)
  middle_aln_2 <- bwa_alignment(unlist(middle$seq), call_bwa = "bwa mem ", samtools_param = "-F 128 -F 4")
  middle_aln_2 <- convertingHtoS(middle_aln_2, unique_id = "QNAME")
  middle_aln_2$unified_cigar <- unify_cigar_strand(middle_aln_2$CIGAR, flag = middle_aln_2$FLAG, to = "+", along_query = TRUE)

  # Annotating the sequencing reads
  # Getting the location of the mapped location in the reads
  mapped_1 <- aln_1[!aln_1$CIGAR=="*",]
  mapping1_read_loc <- unlist(cigarRangesAlongQuerySpace(cigar_convert(mapped_1$unified_cigar, from = "I", to = "M"), ops = "M"))
  elementMetadata(mapping1_read_loc)$annotation <- as.character(sam2gr(mapped_1))
  elementMetadata(mapping1_read_loc)$QNAME <- mapped_1$QNAME

  # For left clipped reads
  left_read_loc <- unlist(cigarRangesAlongQuerySpace(cigar_convert(aln_2$left$unified_cigar, from = "I", to = "M"), ops = "M"))
  elementMetadata(left_read_loc)$QNAME <- aln_2$left$QNAME
  elementMetadata(left_read_loc)$annotation <- as.character(sam2gr(aln_2$left))

  # For right clipped reads
  right_clipped_read_loc <- unlist(cigarRangesAlongQuerySpace(cigar_convert(aln_2$right$unified_cigar, from = "I", to = "M"), ops = "M"))
  right_clipped_len <- cigarWidthAlongQuerySpace(aln_2$right$unified_cigar)
  right_clipped_start_after <- nchar(seq[aln_2$right$QNAME]) - right_clipped_len
  right_read_loc <- IRanges(start = right_clipped_start_after + start(right_clipped_read_loc),
                               end = right_clipped_start_after + end(right_clipped_read_loc))
  elementMetadata(right_read_loc)$QNAME <- aln_2$right$QNAME
  elementMetadata(right_read_loc)$annotation <- as.character(sam2gr(aln_2$right))

  # For middle clipped regions
  middle_grp <- sub("\\.[0-9]+$", "\\1", middle_aln_2$QNAME)
  mid_start <- unlist(middle$start[middle_grp])[middle_aln_2$QNAME]
  mid_end <- unlist(middle$end[middle_grp])[middle_aln_2$QNAME]
  mid_clipped_read_loc <- unlist(cigarRangesAlongQuerySpace(cigar_convert(middle_aln_2$unified_cigar, from = "I", to = "M"), ops = "M"))
  mid_read_loc <- IRanges(start = mid_start + start(mid_clipped_read_loc) - 1,
                          end = mid_start + end(mid_clipped_read_loc) - 1 )
  elementMetadata(mid_read_loc)$annotation <- as.character(sam2gr(middle_aln_2))
  elementMetadata(mid_read_loc)$QNAME <- gsub("\\..*","",middle_aln_2$QNAME)

  # For unmapped reads
  unmapped_read_loc <- unlist(cigarRangesAlongQuerySpace(cigar_convert(aln_2$unmapped$unified_cigar, from = "I", to = "M"), ops = "M"))
  elementMetadata(unmapped_read_loc)$QNAME <- aln_2$unmapped$QNAME
  elementMetadata(unmapped_read_loc)$annotation <- as.character(sam2gr(aln_2$unmapped))

  # Combining all together
  combined_read_loc <- c(mapping1_read_loc, left_read_loc, right_read_loc, mid_read_loc, unmapped_read_loc)
  combined_read_loc <- split(combined_read_loc, factor(elementMetadata(combined_read_loc)$QNAME, levels = names(seq)))
  system.time(anno_out <- mapply(add_back_seq, x = seq, ir = combined_read_loc, SIMPLIFY = FALSE))
  start <- lapply(anno_out, "[[", i = 1)
  end <- lapply(anno_out, "[[", i = 2)
  annotation <- lapply(anno_out, "[[", i = 3)
  seq <- lapply(anno_out, "[[", i = 4)
  grp <- factor(rep(names(start), lengths(start)), levels = names(start))

  return(list(start = unlist(start, use.names = FALSE),
              end = unlist(end, use.names = FALSE),
              annotation = unlist(annotation, use.names = FALSE),
              seq = unlist(seq, use.names = FALSE),
              grp = grp))
}
