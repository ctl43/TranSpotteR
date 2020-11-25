#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors 'elementMetadata<-' Rle
#' @importFrom IRanges coverage slice subsetByOverlaps IntegerList
#' @importFrom csaw mergeWindows
#' @importFrom GenomeInfoDb seqinfo

filter_read_by_coverage <- function(x, min_read = 2,
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
