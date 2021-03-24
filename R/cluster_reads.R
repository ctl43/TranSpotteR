#' @export
#' @importFrom GenomicRanges findOverlaps GRangesList split
#' @importFrom S4Vectors 'elementMetadata<-'
#' @importFrom BiocGenerics strand start end 'strand<-'

clustering_reads <- function(total, min_peak = 5, size_tol = 0, max_reads = 100){
  # min_peak, the minimum coverage that the cluster should have
  # If min_peak = 3
  #1112221112211111 (peak <= 3)     1122333222111 (peak >= 3)
  #   --------                        -----------
  #------             <- filtered   -------        <- included
  #         -------                     ------
  #
  gr <- total[total$is_anchor]
  p <- gr[strand(gr) == "+"]
  m <- gr[strand(gr) == "-"]
  p_cluster <- .filter_read_by_coverage(p, min_peak = min_peak)
  names(p_cluster) <- seq_along(p_cluster)
  m_cluster <- .filter_read_by_coverage(m, min_peak = min_peak)
  names(m_cluster) <- seq_along(m_cluster)
  strand(p_cluster) <- "+"
  strand(m_cluster) <- "-"

  # Down sampling the number of reads to max_reads
  ## For positive clusters
  p_ol <- findOverlaps(p_cluster, p)
  p_group <- p_ol@from
  p <- p[p_ol@to]
  p_is_large <- table(p_group) > max_reads
  p_grp <- data.table(which = seq_along(p), grp = p_group)
  p_large_which <- p_grp[p_grp$grp %in% which(p_is_large)]
  set.seed(20190721)
  p_selected_large <- p_large_which[, list(which = sample(which, max_reads)), by = grp]
  p_ok_group <- which(!p_is_large)
  p_selected <- sort(c(p_selected_large$which, p_grp$which[p_grp$grp %in% p_ok_group]))
  p_selected <- p_grp[p_selected]
  p <- p[p_selected$which]
  p_partner <- total[.get_partner_read_id(p$QNAME_id)]
  p_cluster$members <- split(p, p_selected$grp)
  p_cluster$partners <- split(p_partner, p_selected$grp)

  ## For negative clusters
  m_ol <- findOverlaps(m_cluster, m)
  m_group <- m_ol@from
  m <- m[m_ol@to]
  m_is_large <- table(m_group) > max_reads
  m_grp <- data.table(which = seq_along(m), grp = m_group)
  m_large_which <- m_grp[m_grp$grp %in% which(m_is_large)]
  set.seed(20190721)
  m_selected_large <- m_large_which[, list(which = sample(which, max_reads)), by = grp]
  m_ok_group <- which(!m_is_large)
  m_selected <- sort(c(m_selected_large$which, m_grp$which[m_grp$grp %in% m_ok_group]))
  m_selected <- m_grp[m_selected]
  m <- m[m_selected$which]
  m_partner <- total[.get_partner_read_id(m$QNAME_id)]
  m_cluster$members <- split(m, m_selected$grp)
  m_cluster$partners <- split(m_partner, m_selected$grp)
  combined <- c(p_cluster, m_cluster)
  names(combined) <- seq_along(combined)
  combined
}

.get_partner_read_id <- function(x){
  out <- rep("", length(x))
  is_first <- grepl("_1", x)
  out[is_first] <- sub("_1", "_2", x[is_first])
  out[!is_first] <- sub("_2", "_1", x[!is_first])
  out
}

#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors 'elementMetadata<-' Rle
#' @importFrom IRanges coverage slice subsetByOverlaps IntegerList
#' @importFrom csaw mergeWindows
#' @importFrom GenomeInfoDb seqinfo

.filter_read_by_coverage <- function(x, min_read = 2,
                                    min_peak = 5,
                                    window_tol = 500,
                                    rescue_tol = 500,
                                    rescue = FALSE){
  cvg <- coverage(x)
  sliced <- slice(cvg, min_read)
  max_count <- IntegerList(lapply(sliced, max))
  chr <- Rle(names(sliced), lengths(sliced))
  collapsed <- GRanges(seqnames = chr,
                       IRanges(unlist(start(sliced)), unlist(end(sliced))),
                       seqinfo = seqinfo(x))
  collapsed$max_count <- unlist(max_count)
  peaks <- collapsed[collapsed$max_count >= min_peak]
  if(rescue){
    temp <- collapsed[!collapsed$max_count >= min_peak] # Recusing reads that are around the main clusters
    rescued <- subsetByOverlaps(temp, peaks + rescue_tol)
    collapsed <- c(peaks, rescued)
  }else{
    collapsed <- peaks
  }
  collapsed <- mergeWindows(collapsed, tol = window_tol)$regions # Merge cluster within the tolerated range
  collapsed
}
