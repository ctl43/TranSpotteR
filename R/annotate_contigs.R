#' @export
#' @importFrom S4Vectors 'elementMetadata<-' elementMetadata split
#' @importFrom BiocGenerics start end strand
#' @importFrom BiocParallel bplapply MulticoreParam

# A wrapper for read annotation
annotate_contigs <- function(x,
                             insert = NULL,
                             genome = NULL,
                             customised_annotation = list(anno_polyA),
                             BPPARAM = MulticoreParam(workers = 3L))
  # A wrapper function process both plus and minus stranded read clusters.
  # Written by Cheuk-Ting Law
{
  p <- x[strand(x) == "+"]
  m <- x[strand(x) == "-"]
  p_out <- .internal_annotation(p, BPPARAM = BPPARAM, insert = insert, genome = genome,
                                customised_annotation = customised_annotation,
                                unique_mapq = 30)
  m_out <- .internal_annotation(m, BPPARAM = BPPARAM, insert = insert, genome = genome,
                                customised_annotation = customised_annotation,
                                unique_mapq = 30)
  out <- rbind(p_out, m_out)
  names(out) <- c("contig_detail", "nreads", "cluster_origin")
  out$cluster_region <- as.character(x[out$cluster_origin])
  return(out)
}

##################################
##################################

#' @importFrom Biostrings BStringSet letterFrequency
anno_polyA <- function(x){
  # Annotating the polyA part of the constructed reads
  j <- letterFrequency(BStringSet(x), letters = c("A", "T", "G", 'C'))
  tot <- rowSums(j)
  a_prop <- j[, 1] / tot
  is_pa <- (a_prop > 0.85) & (nchar(x) > 5)
  t_prop <- j[, 2] / tot
  is_pt <- (t_prop > 0.85) & (nchar(x) > 5)
  is_pt <- !grepl(":|N", x) & is_pt # Ignoring sequence with N
  is_pa <- !grepl(":|N", x) & is_pa # Ignoring sequence with N
  x[is_pa] <- paste0("polyA", ":", 1, "-", nchar(x[is_pa]), ":", strand = "+")
  x[is_pt] <- paste0("polyA", ":", 1, "-", nchar(x[is_pt]), ":", strand = "-")
  x
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
  #Adding back unmapped part of the constructed read to the annotation
  if(nrow(occupied_info) == 0){
    return(data.table(start = 1,
                      end = nchar(x),
                      width = nchar(x),
                      QNAME = QNAME,
                      annotation = x,
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
  selected <- data.table(start = need_start,
                         end = need_end,
                         width = need_end - need_start + 1,
                         QNAME = rep(QNAME, length(selected_seq)),
                         annotation = selected_seq,
                         cigar = rep("*", length(selected_seq)))
  occupied_info <- rbind(occupied_info, selected)
  occupied_info$seq <- str_sub(x, start = occupied_info$start, end = occupied_info$end)
  return(occupied_info)
}

#' @export
#' @importFrom IRanges CharacterList
.internal_annotation <-  function(clusters,  BPPARAM = MulticoreParam(workers = 3),
                                  customised_annotation = list(anno_polyA),
                                  insert, genome, unique_mapq = 30)
  # Annotating the clustered reads, the partner reads from the read cluster ,
  # and long contigs that consist of clustered reads and their partner reads.
  # 1. Aligning the sequence to the insert, LINE1 in this case
  # 2. Aligning the remaining part to the genome in this case
  # 3. Annotating the whole constructed reads
  # 4. Identifying the polyA/polyT or customised annotation
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
  flat_cluster_contigs <- unlist(clusters$cluster_contigs)
  flat_partner_contigs <- unlist(clusters$partner_contigs)
  flat_long_contigs <- unlist(clusters$long_contigs)
  cluster_anno <- annotate_seq(flat_cluster_contigs, BPPARAM = BPPARAM, customised_annotation = customised_annotation,
                               insert = insert, genome = genome, unique_mapq = unique_mapq)
  partner_anno <- annotate_seq(flat_partner_contigs, BPPARAM = BPPARAM, customised_annotation = customised_annotation,
                               insert = insert, genome = genome, unique_mapq = unique_mapq)
  long_anno <- annotate_seq(flat_long_contigs, BPPARAM = BPPARAM, customised_annotation = customised_annotation,
                            insert = insert, genome = genome, unique_mapq = unique_mapq)
  cluster_grp <- factor(rep(names(clusters), lengths(clusters$cluster_contigs)), levels = names(clusters))
  partner_grp <- factor(rep(names(clusters), lengths(clusters$partner_contigs)), levels = names(clusters))
  long_grp <- factor(rep(names(clusters), lengths(clusters$long_contigs)), levels = names(clusters))

  # Combining annotation
  cluster_anno <- lapply(cluster_anno, split, f = factor(cluster_anno$QNAME, levels = names(flat_cluster_contigs)))
  cluster_anno <- lapply(cluster_anno, split, f = cluster_grp)
  partner_anno <- lapply(partner_anno, split, f = factor(partner_anno$QNAME, levels = names(flat_partner_contigs)))
  partner_anno <- lapply(partner_anno, split, f = partner_grp)
  long_anno <- lapply(long_anno, split, f = factor(long_anno$QNAME, levels = names(flat_long_contigs)))
  long_anno <- lapply(long_anno, split, f = long_grp)
  tmp_fun <- function(x, y){mapply(c, rep(x, each = length(y)), rep(y, length(x)), SIMPLIFY = FALSE)}

  if(strand == "+"){ # Creating all possible combination of read structures
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
  combined_info <- do.call(function(...)mapply(tmp_fun_2, ..., SIMPLIFY = FALSE), combined_info)
  long_info <- do.call(function(...)mapply(tmp_fun_2, ..., SIMPLIFY = FALSE), long_anno)
  long_nreads <- elementMetadata(clusters$long_contigs)$n_reads
  combined_nreads[lengths(long_nreads) > 0] <- long_nreads[lengths(long_nreads) > 0]
  combined_info[lengths(long_nreads) > 0] <- long_info[lengths(long_nreads) > 0]

  if(length(combined_info)==0){
    return(NULL)
  }

  names(combined_info) <- names(clusters)
  combined_info_grp <- rep(names(combined_info), lengths(combined_info))
  flat_combined_info <- unlist(combined_info, recursive = FALSE, use.names = FALSE)
  flat_combined_nreads <- unlist(combined_nreads)
  combined_anno_dt <- data.table(flat_combined_info, flat_combined_nreads, combined_info_grp)
  return(combined_anno_dt)
}

#' @export
#' @importFrom GenomicAlignments cigarWidthAlongQuerySpace cigarRangesAlongQuerySpace cigarWidthAlongReferenceSpace cigarQNarrow
#' @importFrom S4Vectors split runValue runLength nchar
#' @importFrom IRanges RleList CharacterList
#' @importFrom BiocGenerics start end width 'start<-' 'end<-'

annotate_seq <- function(seq, insert, genome,
                         BPPARAM = MulticoreParam(workers = 3),
                         customised_annotation = NULL,
                         unique_mapq = 20)
  # This function annotates chimeric sequences that consist of any multi-mapped sequences/non-reference sequence and any genomic regions.
  # It first mapped the bait (repeated regions/non-reference sequence), then the genomic regions by bwa
  # Written by Cheuk-Ting Law
{
  if(length(seq) == 0){
    return(data.table())
  }

  if(any(duplicated(names(seq)))){
    stop("Sequence names are duplicated.")
  }

  if(is.null(names(seq))){
    names(seq) <- paste0("read_", seq_along(seq))
  }

  # Mapping to LINE1 sequence
  aln_1 <- bwa_alignment(seq, ref = insert, samtools_param = "")
  aln_1 <- convertingHtoS(aln_1, unique_id = "QNAME")

  # Annotating the reads
  aln_1$unified_cigar <- unify_cigar_strand(aln_1$CIGAR, flag = aln_1$FLAG, to = "+")

  # Getting clipped sequences
  clipped_seq_1 <- get_unmapped_clipped_read(aln_1, include_middle_unmapped = TRUE) #Need to improve

  # Mapping to genome
  middle <- clipped_seq_1$middle
  clipped_seq_1 <- clipped_seq_1$clipped

  # Can be sped up by combining all the reads together and only doing once.
  aln_2 <- bplapply(clipped_seq_1, bwa_alignment, call_bwa = "bwa mem ",
                    samtools_param = "-F 128 -F 4", BPPARAM = BPPARAM, ref = genome)
  aln_2 <- lapply(aln_2, convertingHtoS, unique_id = "QNAME")
  aln_2 <- lapply(aln_2, function(x)x[x$MAPQ > unique_mapq,])
  aln_2 <- lapply(aln_2, function(x){x$unified_cigar <- unify_cigar_strand(x$CIGAR, flag = x$FLAG, to = "+"); x})

  mapping_2 <- lapply(aln_2, sam2gr)
  middle_aln_2 <- bwa_alignment(unlist(middle$seq), call_bwa = "bwa mem ", samtools_param = "-F 128 -F 4", ref = genome)
  middle_aln_2 <- convertingHtoS(middle_aln_2, unique_id = "QNAME")
  middle_aln_2$unified_cigar <- unify_cigar_strand(middle_aln_2$CIGAR, flag = middle_aln_2$FLAG, to = "+")

  # Annotating the sequencing reads
  # Getting the location of the mapped location in the reads
  mapped_1 <- aln_1[!aln_1$CIGAR=="*",]
  mapping1_read_loc <- unlist(cigarRangesAlongQuerySpace(mapped_1$unified_cigar, ops = c("M", "I"), reduce.ranges = TRUE))
  elementMetadata(mapping1_read_loc)$QNAME <- mapped_1$QNAME
  elementMetadata(mapping1_read_loc)$annotation <- as.character(sam2gr(mapped_1))
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
  elementMetadata(mid_read_loc)$QNAME <- gsub("\\..*", "", middle_aln_2$QNAME)
  elementMetadata(mid_read_loc)$annotation <- as.character(sam2gr(middle_aln_2))
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
  for(i in customised_annotation){
    combined_read_loc$annotation <- i(combined_read_loc$annotation)
  }
  return(combined_read_loc)
}
