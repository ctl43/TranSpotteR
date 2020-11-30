#' @export
#' @importFrom S4Vectors 'elementMetadata<-' runValue elementMetadata
#' @importFrom Biostrings DNAStringSet DNAStringSetList reverseComplement
#' @importFrom BiocParallel bplapply MulticoreParam

sequence_construction <- function(x, BPPARAM = MulticoreParam(workers = 10), prevent_overload = TRUE, max_n = 100){
  # prevent overloading (need to be moved to another place)
  # Down sample to an amount that can be handled properly
  if(prevent_overload){
    is_large <- lengths(x$members) >= max_n
    large <- x[is_large]
    selected <- lapply(large$members, function(x){set.seed(20190721); sample(seq_along(x), max_n)})
    large$members <- mapply(function(x, y)x[y], x = large$members, y = selected)
    large$partners <- mapply(function(x, y)x[y], x = large$partners, y = selected)
    x[is_large] <- large
  }
  p <- x[strand(x) == "+"]
  m <- x[strand(x) == "-"]
  p <- .internal_construction(p, BPPARAM = BPPARAM)
  m <- .internal_construction(m, BPPARAM = BPPARAM)
  return(c(p, m))
}

.internal_construction <- function(x, BPPARAM){
  sign <- runValue(strand(x))
  if(length(sign)>1){
    stop()
  }

  if(length(sign)==0){
    return(NULL)
  }

  info <- elementMetadata(unlist(x$partners))
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
  partner_grp <- as.integer(cut(seq_along(partner_seq), breaks = BPPARAM$workers)) # increase the efficiency of using multicore
  partner_seq <- split(partner_seq, partner_grp)
  # system.time(partner_contigs <- bplapply(partner_seq, greedy_scs, BPPARAM = BPPARAM))
  partner_contigs <- bplapply(partner_seq, function(x)lapply(x, greedy_scs), BPPARAM = BPPARAM)
  partner_contigs <- unlist(partner_contigs, recursive = FALSE, use.names = FALSE)

  # Assembling cluster sequences
  cluster_seq <- lapply(x$members, function(x)x$SEQUENCE)
  cluster_grp <- as.integer(cut(seq_along(cluster_seq), breaks = BPPARAM$workers)) # increase the efficiency of using multicore
  cluster_seq <- split(cluster_seq, cluster_grp)
  # cluster_contigs <- bplapply(cluster_seq, greedy_scs, BPPARAM = BPPARAM)
  cluster_contigs <- bplapply(cluster_seq, function(x)lapply(x, greedy_scs), BPPARAM = BPPARAM)
  cluster_contigs <- unlist(cluster_contigs, recursive = FALSE, use.names = FALSE)

  # Assembling both
  has_both <- lengths(partner_contigs) > 0 & lengths(cluster_contigs) > 0
  merged_contigs <- rep(list(c()), length(partner_contigs))
  merged_seq <- mapply(function(x,y)c(x,y),x = partner_contigs[has_both], y = cluster_contigs[has_both], SIMPLIFY = FALSE)
  merged_grp <- as.integer(cut(seq_along(merged_seq), breaks = BPPARAM$workers)) # increase the efficiency of using multicore
  merged_seq <- split(merged_seq, merged_grp)
  # merged_contigs[has_both] <- bplapply(merged_seq, greedy_scs, BPPARAM = BPPARAM)
  tmp <- bplapply(merged_seq, function(x)lapply(x, greedy_scs), BPPARAM = BPPARAM)
  merged_contigs[has_both] <- unlist(tmp, recursive = FALSE, use.names = FALSE)


  # Storing data
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
