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
  out <- bplapply(list(p, m), .internal_annotation, BPPARAM = BPPARAM, customised_annotation = list(is_polyA = is_polyA))
  out <- c(out[[1]], out[[2]])
  names(out) <- seq_along(out)
  return(out)
}

##################################
##################################
#' @export
#' @importFrom Biostrings BStringSet letterFrequency
is_polyA <- function(x){
  j <- letterFrequency(BStringSet(x), letters = c("A", "T", "G", 'C'))
  tot <- rowSums(j)
  a_prop <- j[, 1] / tot
  t_prop <- j[, 2] / tot
  out <- (a_prop > 0.8|t_prop > 0.8) & (nchar(x) > 5)
  out[grepl(":", x)] <- FALSE
  out
}

#' @export
.resolve_overlap <- function(start, end){
  # To deal with overlapping range
  # The one on the left hand side will occupy the place first, then the second one.
  idx <- order(start)
  start <- start[idx]
  end <- end[idx]
  diff <- start[-1] - head(end, -1)
  diff[diff <= 0] <- 1
  start <- c(start[1], diff + head(end, -1))
  start <- start[order(idx)]
  end <- end[order(idx)]
  list(start = start, end = end)
}

#' @export
#' @importFrom stringr 'str_sub'
.add_back_seq <- function(x, occupied_info, QNAME){
  if(nrow(occupied_info) == 0){
    return(data.table(start = 1, end = nchar(x),
                      width = nchar(x),
                      annotation = x,
                      QNAME = QNAME,
                      cigar = "*",
                      seq = x))
  }
  occupied_start <- occupied_info$start
  occupied_end <- occupied_info$end
  diff <- occupied_start[-1] - head(occupied_end, -1)
  diff[diff <= 0] <- 1
  occupied_start <- c(occupied_start[1], diff + head(occupied_end, -1))
  need_start <- c(1, occupied_end + 1)
  need_end <- c(occupied_start - 1, nchar(x))
  selected_seq <- str_sub(x, start = need_start, end = need_end)
  need_start <- need_start[selected_seq!=""]
  need_end <- need_end[selected_seq!=""]
  selected_seq <- selected_seq[selected_seq!=""]
  selected <- data.table(start = need_start, end = need_end, width = need_end - need_start + 1,
                         annotation = selected_seq,
                         QNAME = rep(QNAME, length(selected_seq)),
                         cigar = rep("*", length(selected_seq)))
  occupied_info <- rbind(occupied_info, selected)
  occupied_info$seq <- str_sub(x, start = occupied_info$start, end = occupied_info$end)
  return(occupied_info)
}

#' @export
#' @importFrom IRanges CharacterList
.internal_annotation <-  function(clusters,  BPPARAM = MulticoreParam(workers = 3), customised_annotation = NULL)
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
  cluster_anno <- .annotate_reads(unlist(clusters$cluster_contigs), BPPARAM = BPPARAM, customised_annotation = customised_annotation)
  partner_anno <- .annotate_reads(unlist(clusters$partner_contigs), BPPARAM = BPPARAM, customised_annotation = customised_annotation)
  long_anno <- .annotate_reads(unlist(clusters$long_contigs), BPPARAM = BPPARAM, customised_annotation = customised_annotation)
  cluster_grp <- factor(rep(seq_along(clusters$cluster_contigs), lengths(clusters$cluster_contigs)), levels = seq_along(clusters))
  partner_grp <- factor(rep(seq_along(clusters$partner_contigs), lengths(clusters$partner_contigs)), levels = seq_along(clusters))
  long_grp <- factor(rep(seq_along(clusters$long_contigs), lengths(clusters$long_contigs)), levels = seq_along(clusters))

  # Combining annotation
  cluster_anno <- lapply(cluster_anno, split, f = cluster_anno$QNAME)
  cluster_anno <- lapply(cluster_anno, split, f = cluster_grp)
  partner_anno <- lapply(partner_anno, split, f = partner_anno$QNAME)
  partner_anno <- lapply(partner_anno, split, f = partner_grp)
  long_anno <- lapply(long_anno, split, f = long_anno$QNAME)
  long_anno <- lapply(long_anno, split, f = long_grp)

  tmp_fun <- function(x, y){mapply(c, rep(x, each = length(y)), rep(y, length(x)), SIMPLIFY = FALSE)}
  cluster_anno$QNAME <- lapply(cluster_anno$QNAME, function(z)lapply(z, as.character))
  partner_anno$QNAME <- lapply(partner_anno$QNAME, function(z)lapply(z, as.character))
  long_anno$QNAME <- lapply(long_anno$QNAME, function(z)lapply(z, as.character))

  if(strand=="+"){ # Creating all possible combination of read structures
    combined_info <- mapply(function(x,y)mapply(tmp_fun, x = x, y = y, SIMPLIFY = FALSE),
                            x = cluster_anno, y = partner_anno, SIMPLIFY = FALSE)
    combined_nreads <- IntegerList(mapply(function(x, y)rep(x, each = length(y)) + y,
                                           x = elementMetadata(clusters$cluster_contigs)$n_reads,
                                           y = elementMetadata(clusters$partner_contigs)$n_reads,
                                           SIMPLIFY = FALSE))
  }else{
    combined_info <- mapply(function(x,y)mapply(tmp_fun, x = x, y = y, SIMPLIFY = FALSE),
                            x = partner_anno, y = cluster_anno, SIMPLIFY = FALSE)
    combined_nreads <- IntegerList(mapply(function(x, y)rep(x, each = length(y)) + y,
                                          x = elementMetadata(clusters$partner_contigs)$n_reads,
                                          y = elementMetadata(clusters$cluster_contigs)$n_reads,
                                          SIMPLIFY = FALSE))
  }

  tmp_fun_2 <- function(...){
    mapply(data.table, ..., SIMPLIFY = FALSE)
  }
  combined_info <- do.call(function(...)mapply(tmp_fun_2, ...), combined_info)
  long_info <- do.call(function(...)mapply(tmp_fun_2, ...), long_anno)
  long_nreads <- elementMetadata(clusters$long_contigs)$n_reads
  combined_nreads[lengths(long_nreads) > 0] <- long_nreads[lengths(long_nreads) > 0]
  combined_info[lengths(long_nreads) > 0] <- long_info[lengths(long_nreads) > 0]
  elementMetadata(clusters)$read_annotation <- combined_info
  elementMetadata(clusters)$nreads <- combined_nreads
  clusters
}

#' @export
#' @importFrom GenomicAlignments cigarWidthAlongQuerySpace cigarRangesAlongQuerySpace cigarWidthAlongReferenceSpace cigarQNarrow
#' @importFrom S4Vectors split runValue runLength nchar
#' @importFrom IRanges RleList CharacterList
#' @importFrom BiocGenerics start end width

.annotate_reads <- function(seq, ref = "~/dicky/reference/fasta/line1_reference/hot_L1_polyA.fa",
                            BPPARAM = MulticoreParam(workers = 3), customised_annotation = NULL)
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
  aln_1$unified_cigar <- unify_cigar_strand(aln_1$CIGAR, flag = aln_1$FLAG, to = "+")

  # Getting clipped sequences
  clipped_seq_1 <- get_unmapped_clipped_read(aln_1, include_middle_unmapped = TRUE)

  # Mapping to genome
  middle <- clipped_seq_1$middle
  clipped_seq_1 <- clipped_seq_1$clipped

  # Can be sped up by combining all the reads together and only doing once.
  aln_2 <- bplapply(clipped_seq_1, bwa_alignment, call_bwa = "bwa mem ",
                    samtools_param = "-F 128 -F 4", BPPARAM = BPPARAM)
  aln_2 <- lapply(aln_2, convertingHtoS, unique_id = "QNAME")
  aln_2 <- lapply(aln_2, function(x)x[x$MAPQ > 10,])
  aln_2 <- lapply(aln_2, function(x){x$unified_cigar <- unify_cigar_strand(x$CIGAR, flag = x$FLAG, to = "+"); x})

  mapping_2 <- lapply(aln_2, sam2gr)
  middle_aln_2 <- bwa_alignment(unlist(middle$seq), call_bwa = "bwa mem ", samtools_param = "-F 128 -F 4")
  middle_aln_2 <- convertingHtoS(middle_aln_2, unique_id = "QNAME")
  middle_aln_2$unified_cigar <- unify_cigar_strand(middle_aln_2$CIGAR, flag = middle_aln_2$FLAG, to = "+")

  # Annotating the sequencing reads
  # Getting the location of the mapped location in the reads
  mapped_1 <- aln_1[!aln_1$CIGAR=="*",]
  mapping1_read_loc <- unlist(cigarRangesAlongQuerySpace(mapped_1$unified_cigar, ops = c("M", "I"), reduce.ranges = TRUE))
  elementMetadata(mapping1_read_loc)$annotation <- as.character(sam2gr(mapped_1))
  elementMetadata(mapping1_read_loc)$QNAME <- mapped_1$QNAME
  elementMetadata(mapping1_read_loc)$cigar <- mapped_1$unified_cigar

  # For left clipped reads
  left_read_loc <- unlist(cigarRangesAlongQuerySpace(aln_2$left$unified_cigar, ops = c("M", "I"), reduce.ranges = TRUE))
  elementMetadata(left_read_loc)$QNAME <- aln_2$left$QNAME
  elementMetadata(left_read_loc)$annotation <- as.character(sam2gr(aln_2$left))
  elementMetadata(left_read_loc)$cigar <- aln_2$left$unified_cigar

  # For right clipped reads
  right_clipped_read_loc <- unlist(cigarRangesAlongQuerySpace(aln_2$right$unified_cigar, ops = c("M", "I"), reduce.ranges = TRUE))
  right_clipped_len <- cigarWidthAlongQuerySpace(aln_2$right$unified_cigar)
  right_clipped_start_after <- nchar(seq[aln_2$right$QNAME]) - right_clipped_len
  right_read_loc <- IRanges(start = right_clipped_start_after + start(right_clipped_read_loc),
                               end = right_clipped_start_after + end(right_clipped_read_loc))
  elementMetadata(right_read_loc)$QNAME <- aln_2$right$QNAME
  elementMetadata(right_read_loc)$annotation <- as.character(sam2gr(aln_2$right))
  elementMetadata(right_read_loc)$cigar <- aln_2$right$unified_cigar

  # For middle clipped regions
  middle_grp <- sub("\\.[0-9]+$", "\\1", middle_aln_2$QNAME)
  mid_start <- unlist(middle$start[middle_grp])[middle_aln_2$QNAME]
  mid_end <- unlist(middle$end[middle_grp])[middle_aln_2$QNAME]
  mid_clipped_read_loc <- unlist(cigarRangesAlongQuerySpace(middle_aln_2$unified_cigar, ops = c("M", "I"), reduce.ranges = TRUE))
  mid_read_loc <- IRanges(start = mid_start + start(mid_clipped_read_loc) - 1,
                          end = mid_start + end(mid_clipped_read_loc) - 1 )
  elementMetadata(mid_read_loc)$annotation <- as.character(sam2gr(middle_aln_2))
  elementMetadata(mid_read_loc)$QNAME <- gsub("\\..*","",middle_aln_2$QNAME)
  elementMetadata(mid_read_loc)$cigar <- middle_aln_2$unified_cigar

  # For unmapped reads
  unmapped_read_loc <- unlist(cigarRangesAlongQuerySpace(aln_2$unmapped$unified_cigar, ops = c("M", "I"), reduce.ranges = TRUE))
  elementMetadata(unmapped_read_loc)$QNAME <- aln_2$unmapped$QNAME
  elementMetadata(unmapped_read_loc)$annotation <- as.character(sam2gr(aln_2$unmapped))
  elementMetadata(unmapped_read_loc)$cigar <- aln_2$unmapped$unified_cigar

  # Re-calculating data after resolving overlapped annotation
  combined_read_loc <- c(mapping1_read_loc, left_read_loc, right_read_loc, mid_read_loc, unmapped_read_loc)
  combined_read_loc <- data.table(data.frame(combined_read_loc))
  read_origin <- factor(combined_read_loc$QNAME, levels = names(seq))
  combined_read_loc <- split(combined_read_loc, read_origin)
  resolved <- mapply(.resolve_overlap, start = lapply(combined_read_loc, "[[", i = "start"),
                     end = lapply(combined_read_loc, "[[", i = "end"), SIMPLIFY = FALSE)
  start <- lapply(resolved, "[[", i = "start")
  end <- lapply(resolved, "[[", i = "end")
  resolved <- data.table(start = unlist(start), end = unlist(end))
  combined_read_loc <- do.call(rbind, combined_read_loc)
  tmp_cigar <- sub("^[0-9]+S", "", combined_read_loc$cigar)
  combined_read_loc$cigar <- sub("[0-9]+S$", "", tmp_cigar)
  cigar_start <- resolved$start - combined_read_loc$start + 1
  cigar_end <- cigar_start + resolved$end - resolved$start
  combined_read_loc <- combined_read_loc[!resolved$start > resolved$end,] # in case, there is complete overlapping
  resolved <- resolved[!resolved$start > resolved$end,]
  combined_read_loc$cigar <- cigarQNarrow(combined_read_loc$cigar, start = cigar_start, end = cigar_end)
  combined_read_loc$start <- resolved$start
  combined_read_loc$end <- resolved$end

  # Correcting the annotation cos some read annotation are trimmed in add_back_seq
  is_plus <- grepl(":\\+", combined_read_loc$annotation)
  plus_gr <- convert_character2gr(combined_read_loc$annotation[is_plus])
  start(plus_gr) <- end(plus_gr) - cigarWidthAlongReferenceSpace(combined_read_loc$cigar[is_plus]) + 1
  combined_read_loc$annotation[is_plus] <- as.character(plus_gr)
  is_minus <- grepl(":-", combined_read_loc$annotation)
  minus_gr <- convert_character2gr(combined_read_loc$annotation[is_minus])
  end(minus_gr) <- start(minus_gr) + cigarWidthAlongReferenceSpace(combined_read_loc$cigar[is_minus]) - 1
  combined_read_loc$annotation[is_minus] <- as.character(minus_gr)

  # Adding back those unannotation parts
  combined_read_loc <- combined_read_loc[order(combined_read_loc$start),] #make sure reads are ordered within each read
  combined_read_loc <- combined_read_loc[order(combined_read_loc$QNAME),]
  combined_read_loc <- split(combined_read_loc, factor(combined_read_loc$QNAME, levels = names(seq)))
  combined_read_loc <- mapply(.add_back_seq, x = seq, occupied_info = combined_read_loc, QNAME = names(combined_read_loc),SIMPLIFY = FALSE)
  combined_read_loc <- data.table::rbindlist(combined_read_loc)
  combined_read_loc <- combined_read_loc[order(combined_read_loc$start),] #make sure reads are ordered within each read
  combined_read_loc <- combined_read_loc[order(factor(combined_read_loc$QNAME, levels = names(seq))),]

  # Customized annotation
  if(any(is.null(names(customised_annotation)))){
    names(customised_annotation) <- paste0("customised_anno_", seq_along(customised_annotation))
  }
  extra <- lapply(customised_annotation, function(p)p(combined_read_loc$annotation))
  setDT(extra)
  return(cbind(combined_read_loc, extra))
}

