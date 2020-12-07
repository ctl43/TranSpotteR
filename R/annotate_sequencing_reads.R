#' @export
#' @importFrom S4Vectors 'elementMetadata<-' elementMetadata split
#' @importFrom BiocGenerics start end strand
#' @importFrom BiocParallel bplapply MulticoreParam

# A wrapper for read annotation
annotate_constructed_reads <- function(x, partner_is_anchor = TRUE, BPPARAM = MulticoreParam(workers = 3L)){
  p <- x[strand(x) == "+"]
  m <- x[strand(x) == "-"]
  out <- bplapply(list(p, m), .internal_annotation, BPPARAM = BPPARAM, partner_is_anchor = partner_is_anchor)
  out <- c(out[[1]], out[[2]])
  names(out) <- seq_along(out)
  return(out)
}
##################################
##################################

#' @export
multiple_replacement <- function(x, ir = NULL, start = NULL, end = NULL,  to_replace = NULL, to_replace_element = NULL){
  if(!is.null(ir)){
    start <- start(ir)
    end <- end(ir)
    to_replace <- elementMetadata(ir)[[to_replace_element]]
  }
  n_start <- length(start)
  n_end <- length(end)
  if(n_start == 0 | n_end == 0 ){
    return(as.character(x))
  }

  if(n_start != n_end){
    stop("The length of start and end are not the same.")
  }


  idx <- order(start)
  start <- start[idx]
  end <- end[idx]
  to_replace <- to_replace[idx]
  # To deal with overlapping range
  # The one on the left hand side will occupy the place first, then the second one.
  diff <- start[-1] - head(end, -1)
  diff[diff <= 0] <- 1
  start <- c(start[1], diff + head(end, -1))
  to_replace[start > end] <- ""
  need_start <- c(1, end + 1)
  need_end <- c(start - 1, nchar(x))
  n <- length(need_start)
  collected <- c()
  for(i in seq_len(n - 1)){
    # Can be mapply
    collected <- c(collected, substr(x, need_start[i], need_end[i]))
    collected <- c(collected, to_replace[i])
  }
  out <- paste(collected[collected!=""], collapse = " ")
  sub("^ ", "", out)
}


#' @export
#' @importFrom IRanges CharacterList
.internal_annotation <-  function(clusters, partner_is_anchor = FALSE,  BPPARAM = MulticoreParam(workers = 3)){
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

  cluster_grp <- rep(seq_along(clusters$cluster_contigs), lengths(clusters$cluster_contigs))
  partner_grp <- rep(seq_along(clusters$partner_contigs), lengths(clusters$partner_contigs))
  long_grp <- rep(seq_along(clusters$long_contigs), lengths(clusters$long_contigs))

  # Annotating read cluster with anchor regions
  clusters$has_long <- lengths(clusters$long_contigs) != 0

  # Combining annotation
  cluster_anno <- split(cluster_anno, factor(cluster_grp, levels = seq_along(clusters)))
  partner_anno <- split(partner_anno, factor(partner_grp, levels = seq_along(clusters)))
  long_anno <- CharacterList(split(long_anno, factor(long_grp, levels = seq_along(clusters))))
  both_mapped <- lengths(cluster_anno) > 0 & lengths(partner_anno) > 0

  if(strand=="+"){ # Creating all possible combination of read structures
    combined_anno <- CharacterList(mapply(function(x, y){paste0(rep(x, each = length(y)), " NNNNN ",y)},
                                          x = cluster_anno,
                                          y = partner_anno,
                                          SIMPLIFY = FALSE))
    combined_n_reads <- IntegerList(mapply(function(x, y){rep(x, each = length(y)) + y},
                                           x = elementMetadata(clusters$cluster_contigs)$n_reads,
                                           y = elementMetadata(clusters$partner_contigs)$n_reads,
                                           SIMPLIFY = FALSE))
  }else{
    combined_anno <- CharacterList(mapply(function(x, y){paste0(rep(x, each = length(y)), " NNNNN ",y)},
                                          x = partner_anno,
                                          y = cluster_anno,
                                          SIMPLIFY = FALSE))
    combined_n_reads <- IntegerList(mapply(function(x, y){rep(x, each = length(y)) + y},
                                                                 x = elementMetadata(clusters$partner_contigs)$n_reads,
                                                                 y = elementMetadata(clusters$cluster_contigs)$n_reads,
                                                                 SIMPLIFY = FALSE))
  }
  n_reads <- elementMetadata(clusters$long_contigs)$n_reads
  n_reads[lengths(n_reads) == 0] <- combined_n_reads[lengths(n_reads) == 0]
  combined_anno[any(combined_anno == " NNNNN ")] <- CharacterList(character(0))
  combined_anno[!both_mapped] <- CharacterList(character(0))
  combined_storage <- CharacterList(mapply(c, long_anno, combined_anno, SIMPLIFY = FALSE))
  elementMetadata(combined_storage)$n_reads <- n_reads
  clusters$read_annotation <- combined_storage
  clusters
}

#' @importFrom  GenomicAlignments cigarToRleList
cigar_convert <- function(cigar_string, from, to){
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
                            BPPARAM = MulticoreParam(workers = 3)){
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
  elementMetadata(mapping1_read_loc)$mapped_regions <- as.character(sam2gr(mapped_1))
  elementMetadata(mapping1_read_loc)$QNAME <- mapped_1$QNAME

  # For left clipped reads
  left_read_loc <- unlist(cigarRangesAlongQuerySpace(cigar_convert(aln_2$left$unified_cigar, from = "I", to = "M"), ops = "M"))
  elementMetadata(left_read_loc)$QNAME <- aln_2$left$QNAME
  elementMetadata(left_read_loc)$mapped_regions <- as.character(sam2gr(aln_2$left))

  # For right clipped reads
  right_clipped_read_loc <- unlist(cigarRangesAlongQuerySpace(cigar_convert(aln_2$right$unified_cigar, from = "I", to = "M"), ops = "M"))
  right_clipped_len <- cigarWidthAlongQuerySpace(aln_2$right$unified_cigar)
  right_clipped_start_after <- nchar(seq[aln_2$right$QNAME]) - right_clipped_len
  right_read_loc <- IRanges(start = right_clipped_start_after + start(right_clipped_read_loc),
                               end = right_clipped_start_after + end(right_clipped_read_loc))
  elementMetadata(right_read_loc)$QNAME <- aln_2$right$QNAME
  elementMetadata(right_read_loc)$mapped_regions <- as.character(sam2gr(aln_2$right))

  # For middle clipped regions
  middle_grp <- sub("\\.[0-9]+$", "\\1", middle_aln_2$QNAME)
  mid_start <- unlist(middle$start[middle_grp])[middle_aln_2$QNAME]
  mid_end <- unlist(middle$end[middle_grp])[middle_aln_2$QNAME]
  mid_clipped_read_loc <- unlist(cigarRangesAlongQuerySpace(cigar_convert(middle_aln_2$unified_cigar, from = "I", to = "M"), ops = "M"))
  mid_read_loc <- IRanges(start = mid_start + start(mid_clipped_read_loc) - 1,
                          end = mid_start + end(mid_clipped_read_loc) - 1 )
  elementMetadata(mid_read_loc)$mapped_regions <- as.character(sam2gr(middle_aln_2))
  elementMetadata(mid_read_loc)$QNAME <- gsub("\\..*","",middle_aln_2$QNAME)

  # For unmapped reads
  unmapped_read_loc <- unlist(cigarRangesAlongQuerySpace(cigar_convert(aln_2$unmapped$unified_cigar, from = "I", to = "M"), ops = "M"))
  elementMetadata(unmapped_read_loc)$QNAME <- aln_2$unmapped$QNAME
  elementMetadata(unmapped_read_loc)$mapped_regions <- as.character(sam2gr(aln_2$unmapped))

  # Combining all together
  combined_read_loc <- c(mapping1_read_loc, left_read_loc, right_read_loc, mid_read_loc, unmapped_read_loc)
  combined_read_loc <- split(combined_read_loc, factor(elementMetadata(combined_read_loc)$QNAME, levels = names(seq)))
  anno_out <- mapply(multiple_replacement, x = seq, ir = combined_read_loc,  to_replace_element = "mapped_regions")
  return(anno_out)
}
