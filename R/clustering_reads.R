clustering_reads <- function(gr, min_peak=5, size_tol=0){
  p <- gr[strand(gr)=="+"]
  m <- gr[strand(gr)=="-"]
  p_cluster <- filter_read_by_coverage(p, min_peak=min_peak, size_tol = size_tol)
  names(p_cluster) <- seq_along(p_cluster)
  m_cluster <- filter_read_by_coverage(m, min_peak=min_peak, size_tol = size_tol)
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
  list(p=p_cluster, m=m_cluster)
}
