#' @export
#' @importFrom GenomicRanges GRanges split
#' @importFrom S4Vectors 'elementMetadata<-' Rle runValue elementMetadata
#' @importFrom IRanges subsetByOverlaps IntegerList
#' @importFrom csaw mergeWindows
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom Biostrings DNAStringSet DNAStringSetList reverseComplement
#' @importFrom BiocParallel bplapply MulticoreParam

sequence_construction <- function(x, BPPARAM = MulticoreParam(workers = 20)){
  p <- x[strand(x) == "+"]
  m <- x[strand(x) == "-"]
  p <- .internal_construction(p, BPPARAM = BPPARAM)
  m <- .internal_construction(m, BPPARAM = BPPARAM)
  return(c(p, m))
}

.internal_construction <- function(x, BPPARAM){
  sign <- runValue(strand(x))
  if(length(sign)!=1){
    stop()
  }

  info <- S4Vectors::elementMetadata(unlist(x$partners))
  partner_seq <- info$SEQUENCE
  partner_seq <- DNAStringSet(partner_seq)
  if(sign == "+"){
    partner_seq[!info$is_rc] <- reverseComplement(partner_seq[!info$is_rc])
  }else{
    partner_seq[info$is_rc] <- reverseComplement(partner_seq[info$is_rc])
  }
  partner_seq <- as.character(partner_seq)
  grp <- rep(seq_along(x$partners), lengths(x$partners))
  grp <- factor(grp, levels = seq_along(x$partners))
  partner_seq <- split(partner_seq, grp)

  # Assembling partner seq
  partner_contigs <- bplapply(partner_seq, greedy_scs, BPPARAM = BPPARAM)

  # Assembling cluster sequences
  grp <- rep(seq_along(x$members), lengths(x$members))
  grp <- factor(grp, levels = seq_along(x$members))
  cluster_seq <- split(unlist(x$members)$SEQUENCE, grp)
  cluster_contigs <- bplapply(cluster_seq, greedy_scs, BPPARAM = BPPARAM)

  # Assembling both
  has_both <- lengths(partner_contigs) > 0 & lengths(cluster_contigs) > 0
  merged_contigs <- rep(list(c()), length(partner_contigs))
  merged_seq <- mapply(function(x,y)c(x,y),x = partner_contigs[has_both], y = cluster_contigs[has_both], SIMPLIFY = FALSE)
  merged_contigs[has_both] <- bplapply(merged_seq, greedy_scs, BPPARAM = BPPARAM)
  names(merged_contigs) <- names(partner_contigs)

  # Storing data
  names(cluster_contigs) <- names(partner_contigs) <- names(merged_contigs) <- seq_along(x)
  partner_contigs <- DNAStringSetList(lapply(partner_contigs, DNAStringSet)) # to preserve the sequence names
  cluster_contigs <- DNAStringSetList(lapply(cluster_contigs, DNAStringSet))
  merged_contigs <- DNAStringSetList(lapply(merged_contigs, DNAStringSet))
  has_long <- lengths(merged_contigs) > 0
  x$long_contigs <- merged_contigs
  partner_contigs[has_long] <- DNAStringSetList(DNAStringSet())
  x$partner_contigs <- partner_contigs
  cluster_contigs[has_long] <- DNAStringSetList(DNAStringSet())
  x$cluster_contigs <- cluster_contigs
  x
}
