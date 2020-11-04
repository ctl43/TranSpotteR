# potential transduction cases
adding_filter2pair <- function(p){
  # p <- ol_db[4]
  elementMetadata(p) <- NULL
  not_only_polyA <- !mapply(function(x,y)all(c(start(x), start(y))>6000), x=p@first$combined_bait, y=p@second$combined_bait) # removing those clusters that seem making sense but actually only polyA
  make_sense <- !((end(p@first)-start(p@second))<(-1000)| # not far away
                    seqnames(p@first)!=seqnames(p@second)& # on the same chromosome
                    (p@first$has_long&p@second$has_long)) # have long contigs
  make_sense <- as.logical(make_sense)
  p_is_correct <- lengths(p@first$anchor_left_mapping2genome)>0|
    lengths(p@first$anchor_unmapped_mapping2genome)>0|
    lengths(p@first$long_left_mapping2genome)>0
  m_is_correct <- lengths(p@second$anchor_right_mapping2genome)>0|
    lengths(p@second$anchor_unmapped_mapping2genome)>0|
    lengths(p@second$long_right_mapping2genome)>0
  
  elementMetadata(p)$p_is_correct <- p_is_correct
  elementMetadata(p)$m_is_correct <- m_is_correct
  elementMetadata(p)$not_only_polyA <- not_only_polyA
  elementMetadata(p)$make_sense <- make_sense
  
  pass <- rowSums(as.matrix(elementMetadata(p)))==ncol(elementMetadata(p))
  
  anchor_left_insert <- end(p@first$anchor_left_mapping2genome)
  anchor_right_insert <- start(p@second$anchor_right_mapping2genome)
  anchor_left_insert <- max(anchor_left_insert)
  anchor_right_insert <- max(anchor_right_insert)
  anchor_left_insert[anchor_left_insert<0] <- NA
  anchor_right_insert[anchor_right_insert<0] <- NA
  
  long_left_insert <- end(p@first$long_left_mapping2genome)
  long_right_insert <- start(p@second$long_right_mapping2genome)
  long_left_insert <- max(long_left_insert)
  long_right_insert <- max(long_right_insert)
  long_left_insert[long_left_insert<0] <- NA
  long_right_insert[long_right_insert<0] <- NA
  
  left_exact <- long_left_insert
  right_exact <- long_right_insert
  left_exact[is.na(left_exact)] <- anchor_left_insert[is.na(left_exact)]
  right_exact[is.na(right_exact)] <- anchor_right_insert[is.na(right_exact)]
  
  left_approx <- max(end(p@first$anchor_unmapped_mapping2genome))
  left_approx[left_approx<0] <- NA
  right_approx <- max(end(p@second$anchor_unmapped_mapping2genome))
  right_approx[right_approx<0] <- NA
  right_approx[!is.na(right_exact)] <- NA
  left_approx[!is.na(left_exact)] <- NA
  
  elementMetadata(p)$chr <- as.character(seqnames(p@first))
  elementMetadata(p)$left_exact <- left_exact
  elementMetadata(p)$right_exact <- right_exact
  elementMetadata(p)$left_approx <- left_approx
  elementMetadata(p)$right_approx <- right_approx
  elementMetadata(p)$pass <- pass
  
  # selected <- p[pass&p@first$has_long&p@second$has_long]
  # combined <- mcmapply(function(x,y)single_join_reads(x,y), 
  #                         x=selected@first$long_contigs, 
  #                         y=selected@second$long_contigs,
  #                         SIMPLIFY = FALSE, mc.cores = 10)
  # combined <- unlist(DNAStringSetList(combined))
  # anno <- .annotate_reads(combined)
  
  # test <- p[unique_strand_only&pass]
  right_bait <- max(end(p@second$combined_bait))
  right_bait_strand <- sapply(strand(p@second$combined_bait), function(x)as.character(x)[1])
  right_bait[right_bait_strand=="-"&!is.na(right_bait_strand)] <- min(start(p@second$combined_bait))[right_bait_strand=="-"&!is.na(right_bait_strand)]
  elementMetadata(p)$right_bait <- right_bait
  
  left_bait <- min(start(p@first$combined_bait))
  left_bait_strand <- sapply(strand(p@first$combined_bait), function(x)as.character(x)[1])
  left_bait[left_bait_strand=="-"&!is.na(left_bait_strand)] <- max(end(p@first$combined_bait))[left_bait_strand=="-"&!is.na(left_bait_strand)]
  elementMetadata(p)$left_bait <- left_bait
  
  unique_strand_only <- mapply(function(x,y)length(unique(as.character(strand(c(x,y)))))==1, x=p@first$combined_bait, y=p@second$combined_bait)
  elementMetadata(p)$unique_strand_only <- unique_strand_only
  left_pos <- p@elementMetadata$left_exact
  left_pos[is.na(left_pos)] <- p@elementMetadata$left_approx[is.na(left_pos)]
  elementMetadata(p)$left_pos <- left_pos
  right_pos <- p@elementMetadata$right_exact
  right_pos[is.na(right_pos)] <- p@elementMetadata$right_approx[is.na(right_pos)]
  elementMetadata(p)$right_pos <- right_pos
  
  p_on_left <- start(p@first)<=start(p@second)
  elementMetadata(p)$p_on_left <- p_on_left

  p
  
}
