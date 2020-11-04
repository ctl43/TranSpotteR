sequence_construction <- function(x){
  sign <- runValue(strand(x)) 
  if(length(sign)!=1){
    stop()
  }
  info <- do.call(rbind, x$partners)
  partner_seq <- info$SEQUENCE
  partner_seq <- DNAStringSet(partner_seq)
  if(sign=="+"){
    partner_seq[!info$is_rc] <- reverseComplement(partner_seq[!info$is_rc])
  }else{
    partner_seq[info$is_rc] <- reverseComplement(partner_seq[info$is_rc])
  }
  
  grp <- rep(seq_along(x$partners), sapply(x$partners, nrow))
  grp <- factor(grp, levels = seq_along(x$partners))
  partner_seq <- split(partner_seq, grp)
  
  # Assembling partner seq
  partner_contigs <- DNAStringSetList(mclapply(partner_seq, function(y).assembly(y)[[2]], mc.cores = 60L))
  
  # Assembling cluster sequences
  grp <- rep(seq_along(x$members), lengths(x$members))
  grp <- factor(grp, levels = seq_along(x$members))
  cluster_seq <- split(DNAStringSet(unlist(x$members)$SEQUENCE), grp)
  cluster_contigs <- DNAStringSetList(mclapply(cluster_seq, function(y).assembly(y)[[2]], mc.cores = 60L))

  # Assembling both
  merged_seq <- mapply(function(x,y)c(x,y),x=partner_contigs, y=cluster_contigs)
  merged_contigs <- DNAStringSetList(mclapply(merged_seq, function(y).assembly(y)[[2]], mc.cores = 60L))
  
  # Attempting to constructing the break point between LINE1 and genome
  if(sign=="+"){
    combined_contigs <- DNAStringSetList(mcmapply(join_reads, p=cluster_contigs, q=partner_contigs, min_overlap=10, mc.cores = 3L)) # Be careful about the direction
  }else{
    combined_contigs <- DNAStringSetList(mcmapply(join_reads, p=partner_contigs, q=cluster_contigs, min_overlap=10, mc.cores = 3L)) # Be careful about the direction
  }

  # Finding the longer reads as final assembled reads
  pmax_len <- pmax(max(nchar(cluster_contigs)), max(nchar(partner_contigs)))
  is_improved <- max(nchar(merged_contigs))>(pmax_len+30)
  has_combined <- lengths(combined_contigs)!=0
  merged_contigs[!is_improved] <- DNAStringSetList(DNAStringSet())
  merged_better <- max(nchar(merged_contigs))>=max(nchar(combined_contigs))
  selected_contigs <- rep(DNAStringSetList(DNAStringSet()), length(merged_better))
  selected_contigs[merged_better] <- merged_contigs[merged_better]
  selected_contigs[!merged_better] <- combined_contigs[!merged_better]
  has_longer <- lengths(selected_contigs)>0

  # Storing data
  names(cluster_contigs) <- names(partner_contigs) <- names(selected_contigs) <- seq_along(x)
  x$long_contigs <- selected_contigs
  partner_contigs[has_longer] <- DNAStringSetList(DNAStringSet())
  x$partner_contigs <- partner_contigs
  cluster_contigs[has_longer] <- DNAStringSetList(DNAStringSet())
  x$cluster_contigs <- cluster_contigs
  x
}
