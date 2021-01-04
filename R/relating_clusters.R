#' @export
#' @importFrom GenomicRanges findOverlaps
#' @importFrom igraph graph_from_data_frame components
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges CharacterList
#' @importFrom S4Vectors grepl
#' @importFrom BiocGenerics match

relating_cluster <- function(info, insert_seqname = "Hot_L1_polyA", overlap_distance = 1000){
  # Getting the anchor information
  anno <- info[[1]]
  anno <- CharacterList(lapply(anno, "[[", i = "annotation"))
  anno_ranges <- anno[S4Vectors::grepl(":",anno)]
  gr <- convert_character2gr(unlist(anno_ranges))
  grp <- factor(rep(info$cluster_origin, lengths(anno_ranges)), levels = unique(info$cluster_origin)) # complex grp
  gr <- S4Vectors::split(gr, grp)
  gr <- unique(gr)

  # Groupping related anchors
  anchors <- gr[!(BiocGenerics::match(seqnames(gr), insert_seqname, nomatch = 0) > 0)]
  ol <- findOverlaps(anchors, anchors + overlap_distance, ignore.strand = TRUE)
  df <- cbind.data.frame(ol@from, ol@to)
  loc <- unique(unlist(df))
  graph <- graph_from_data_frame(df, directed = FALSE)
  cluster_origin_with_grp <- unique(info$cluster_origin)[loc]
  origin_counts <- table(factor(info$cluster_origin, levels = unique(info$cluster_origin)))
  group <- rep(as.integer(components(graph)$membership), origin_counts[cluster_origin_with_grp])
  info$group <- NA
  info$group[info$cluster_origin %in% cluster_origin_with_grp] <- group
  is_solo <- is.na(info$group)
  solo_count <- table(factor(info$cluster_origin[is_solo], levels = unique(info$cluster_origin)))
  solo_grp <- max(as.integer(components(graph)$membership)) + seq_len(sum(solo_count != 0))
  solo_grp <- rep(solo_grp, solo_count[solo_count>0])
  info$group[is_solo] <-  solo_grp# filling up the solo group
  info
}
