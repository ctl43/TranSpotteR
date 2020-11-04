filter_clusters <- function(x){
  # Flattening the data before filtering
  anno_grp <- rep(x$group, lengths(x$combined_read_anno))
  flat_anno <- unlist(x$combined_read_anno)
  
  ## Filtering group without information
  ok_a <- sum(grepl(":", flat_anno)|grepl("[A-Z]", flat_anno))>1  # To filter contig that fully mapped to a single region with no extra information
  ok_b <- !(any(grepl("Hot_L1_polyA", flat_anno))&sum(grepl(":", flat_anno))==1) # To filter contig only mapping to LINE1
  gr <- convert_character2gr(flat_anno)
  is_long <- unname(!any(grepl("NNNNN", flat_anno)))
  elementMetadata(gr)$is_long <- is_long
  elementMetadata(gr)$anno <- flat_anno
  elementMetadata(gr)$grp <- anno_grp
  elementMetadata(gr)$count <- rep(x$count, lengths(x$combined_read_anno))
  elementMetadata(gr)$peak_count <- rep(x$peak_count, lengths(x$combined_read_anno))
  
  gr_start <- range(start(range(gr)))
  gr_start[abs(gr_start[,1]) > 999999999,] <- 0
  ok_c <- !((gr_start[,2]-gr_start[,1]<=50000)&lengths(gr)>1) # To filter contig that map to close regions (unlikely to be transduction event)
  ok <- ok_a&ok_b&ok_c
  gr <- gr[ok]
  
  # Filtering groups that are clearly not insertion
  non_orphan_grp <- unique(elementMetadata(gr)$grp[duplicated(elementMetadata(gr)$grp)])
  not_only_polyA <- any(start(gr[seqnames(gr)=="Hot_L1_polyA"])<6000)
  not_only_polyA <- elementMetadata(gr)$grp[not_only_polyA]
  
  ok <- (elementMetadata(gr)$grp%in%non_orphan_grp)&(elementMetadata(gr)$grp%in%not_only_polyA)
  gr <- gr[ok]
}