relating_cluster <- function(cluster){

  anchors <- cluster$combined_anchor
  names(anchors) <- seq_along(anchors)
  origins <- rep(seq_along(anchors), lengths(anchors))
  
  pool_anchor <- unlist(anchors, use.names = FALSE)
  pool_anchor$origins <- origins
  
  unique_id <- paste0(origins, "_", names(pool_anchor))
  tmp_grp <- factor(unique_id, levels=unique(unique_id))
  
  pool_anchor <- split(pool_anchor, tmp_grp)
  ol <- findOverlaps(pool_anchor, pool_anchor+1000, ignore.strand=TRUE)
  origins_2 <- as.integer(gsub("_.*","", names(pool_anchor)))
  df <- cbind.data.frame(origins_2[ol@from], origins_2[ol@to])
  loc <- unique(unlist(df))
  graph <- graph_from_data_frame(df, directed = FALSE)

  cluster$group[loc] <- as.integer(components(graph)$membership)
  cluster
}