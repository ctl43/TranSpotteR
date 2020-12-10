#' @export
#' @importFrom GenomicRanges findOverlaps
#' @importFrom igraph graph_from_data_frame components
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges CharacterList
#' @importFrom S4Vectors grepl
#' @importFrom BiocGenerics match

relating_cluster <- function(cluster, insert_seqname = "Hot_L1_polyA", overlap_distance = 1000){
  # Getting the anchor information
  grp <- factor(rep(seq_along(cluster), lengths(cluster$read_annotation)), levels = seq_along(cluster))
  anno <- unname(unlist(cluster$read_annotation))
  # sub_grp <- factor(rep(seq_along(anno), lengths(anno)), levels = seq_along(anno))
  anno <- CharacterList(lapply(anno, function(x)elementMetadata(x)$anno))
  anno_ranges <- anno[S4Vectors::grepl(":",anno)]
  tmp_grp <- factor(rep(seq_along(anno_ranges), lengths(anno_ranges)), levels = seq_along(anno_ranges))
  gr <- unlist(convert_character2gr(unlist(anno_ranges)))
  grp <- rep(grp, table(tmp_grp)) # complex grp
  gr <- S4Vectors::split(gr, grp)
  gr <- unique(gr)

  # Groupping related anchors
  anchors <- gr[!(BiocGenerics::match(seqnames(gr), insert_seqname, nomatch = 0)>0)]
  names(anchors) <- seq_along(anchors)
  ol <- findOverlaps(anchors, anchors + overlap_distance, ignore.strand = TRUE)
  df <- cbind.data.frame(ol@from, ol@to)
  loc <- unique(unlist(df))
  graph <- graph_from_data_frame(df, directed = FALSE)
  elementMetadata(cluster)$group <- NA
  elementMetadata(cluster)$group[loc] <- as.integer(components(graph)$membership)
  is_solo <- is.na(elementMetadata(cluster)$group)
  elementMetadata(cluster)$group[is_solo] <- max(as.integer(components(graph)$membership)) + seq_len(sum(is_solo)) # filling up the solo group
  elementMetadata(cluster)$gr <- gr
  cluster
}
