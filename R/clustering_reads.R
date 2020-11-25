#' @export
#' @importFrom GenomicRanges findOverlaps GRangesList split
#' @importFrom S4Vectors 'elementMetadata<-'
#' @importFrom BiocGenerics strand start end 'strand<-'

clustering_reads <- function(total, min_peak = 5, size_tol = 0){
  gr <- total[total$is_anchor | total$is_disc]
  p <- gr[strand(gr)=="+"]
  m <- gr[strand(gr)=="-"]
  p_cluster <- filter_read_by_coverage(p, min_peak = min_peak)
  names(p_cluster) <- seq_along(p_cluster)
  m_cluster <- filter_read_by_coverage(m, min_peak = min_peak)
  names(m_cluster) <- seq_along(m_cluster)
  strand(p_cluster) <- "+"
  strand(m_cluster) <- "-"

  # Identifying the reads members of the anchor clusters
  ## For positive clusters
  p_ol <- findOverlaps(p_cluster, p)
  p_group <- p_ol@from
  names(p_group) <- p[p_ol@to]$QNAME
  p_cluster$members <- split(p[p_ol@to], p_group)

  ## For negative clusters
  m_ol <- findOverlaps(m_cluster, m)
  m_group <- m_ol@from
  names(m_group) <- m[m_ol@to]$QNAME
  m_cluster$members <- split(m[m_ol@to], m_group)
  combined <- c(p_cluster, m_cluster)

  # Getting partner reads
  partner_ids <- .get_partner_read_id(unlist(combined$members, use.names = FALSE)$QNAME_id)
  combined$partners <- split(total[partner_ids], rep(seq_along(combined), lengths(combined$members)))
  combined
}

.get_partner_read_id <- function(x){
  out <- rep("", length(x))
  is_first <- grepl("_1", x)
  out[is_first] <- sub("_1", "_2", x[is_first])
  out[!is_first] <- sub("_2", "_1", x[!is_first])
  out
}
