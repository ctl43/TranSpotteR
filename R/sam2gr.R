#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomicAlignments cigarWidthAlongReferenceSpace


sam2gr <- function(mapping, seqinfo=NULL){
  if(is.null(mapping)){
    return(GRanges())
  }
  gr <- GRanges(seqnames = mapping$RNAME,
                IRanges(mapping$POS, width=cigarWidthAlongReferenceSpace(mapping$CIGAR)),
                strand = ifelse(bitwAnd(mapping$FLAG, 0x10), "-", "+"),
                seqinfo = seqinfo)
  names(gr) <- mapping$QNAME
  gr
}
