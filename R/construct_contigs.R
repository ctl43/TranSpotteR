#' @export
#' @importFrom S4Vectors 'elementMetadata<-' runValue elementMetadata
#' @importFrom Biostrings DNAStringSet DNAStringSetList reverseComplement
#' @importFrom BiocParallel bplapply MulticoreParam

construct_contigs <- function(x, BPPARAM = MulticoreParam(workers = 10)){
  p <- x[strand(x) == "+"]
  m <- x[strand(x) == "-"]
  p <- .internal_construction(p, BPPARAM = BPPARAM)
  m <- .internal_construction(m, BPPARAM = BPPARAM)
  return(c(p, m))
}

#' @export
.internal_construction <- function(x, BPPARAM){
  # Constructing reads for clustered reads, partner reads and then long contig
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
  partner_contigs <- bplapply(partner_seq, function(x)lapply(x, assemble_reads), BPPARAM = BPPARAM)
  partner_contigs <- unlist(partner_contigs, recursive = FALSE, use.names = FALSE)

  # Assembling cluster sequences
  cluster_seq <- split(unlist(x$members)$SEQUENCE, rep(seq_along(x), lengths(x$members)))
  cluster_grp <- as.integer(cut(seq_along(cluster_seq), breaks = BPPARAM$workers)) # increase the efficiency of using multicore
  cluster_seq <- split(cluster_seq, cluster_grp)
  cluster_contigs <- bplapply(cluster_seq, function(x)lapply(x, assemble_reads), BPPARAM = BPPARAM)
  cluster_contigs <- unlist(cluster_contigs, recursive = FALSE, use.names = FALSE)

  ## Tidying up the number reads
  merged_contigs <- rep(CharacterList(character(0L)), length(partner_contigs))
  elementMetadata(merged_contigs)$n_reads <- rep(IntegerList(integer()),length(merged_contigs))
  partner_n_reads <- IntegerList(lapply(partner_contigs, "[[", i = 2))
  cluster_n_reads <- IntegerList(lapply(cluster_contigs, "[[", i = 2))
  partner_contigs <- CharacterList(lapply(partner_contigs, "[[", i = 1))
  cluster_contigs <- CharacterList(lapply(cluster_contigs, "[[", i = 1))
  elementMetadata(partner_contigs)$n_reads <- partner_n_reads
  elementMetadata(cluster_contigs)$n_reads <- cluster_n_reads

  # Assembling both clustered reads and their partner reads
  ## Requiring to have both cluster contig and partner contig
  has_both <- lengths(partner_contigs) > 0 & lengths(cluster_contigs) > 0
  partner_contigs[!has_both] <- CharacterList(character(0L))
  cluster_contigs[!has_both] <- CharacterList(character(0L))
  partner_n_reads[!has_both] <- IntegerList(integer(0L))
  cluster_n_reads[!has_both] <- IntegerList(integer(0L))
  elementMetadata(partner_contigs)$n_reads <- partner_n_reads
  elementMetadata(cluster_contigs)$n_reads <- cluster_n_reads

  if(sum(has_both) > 0){
    merged_seq <- mapply(function(x, y, p, q)list(seq = c(x, y), n_reads = c(p, q)),
                         x = partner_contigs[has_both],
                         y = cluster_contigs[has_both],
                         p = partner_n_reads[has_both],
                         q = cluster_n_reads[has_both],
                         SIMPLIFY = FALSE)
    merged_grp <- as.integer(cut(seq_along(merged_seq), breaks = BPPARAM$workers)) # increase the efficiency of using multicore
    merged_seq <- split(merged_seq, merged_grp)
    tmp <- bplapply(merged_seq, function(x)lapply(x, function(x)assemble_reads(x[[1]], x[[2]])), BPPARAM = BPPARAM)
    tmp <- unlist(tmp, recursive = FALSE, use.names = FALSE)
    tmp_merged_n_reads <- IntegerList(lapply(tmp, "[[", i = 2))
    merged_contigs[has_both] <- CharacterList(lapply(tmp, "[[", i = 1))
    elementMetadata(merged_contigs)$n_reads <- IntegerList(integer(0L))
    elementMetadata(merged_contigs)$n_reads[has_both] <- tmp_merged_n_reads
  }

  # Storing data
  has_long <- lengths(merged_contigs) > 0
  x$has_long <- has_long
  x$long_contigs <- merged_contigs
  partner_contigs[has_long] <- CharacterList(character(0L))
  elementMetadata(partner_contigs)[["n_reads"]][has_long] <- IntegerList(integer(0L))
  x$partner_contigs <- partner_contigs
  cluster_contigs[has_long] <- CharacterList(character(0L))
  elementMetadata(cluster_contigs)[["n_reads"]][has_long] <- IntegerList(integer(0L))
  x$cluster_contigs <- cluster_contigs
  x
}
