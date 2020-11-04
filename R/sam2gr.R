sam2gr <- function(mapping, seqinfo=NULL){
  if(is.null(mapping)){
    return(GRanges())
  }
  gr <- GenomicRanges::GRanges(seqnames=mapping$RNAME, 
                IRanges::IRanges(mapping$POS, width=GenomicAlignments::cigarWidthAlongReferenceSpace(mapping$CIGAR)),
                strand=ifelse(bitwAnd(mapping$FLAG, 0x10), "-", "+"), 
                seqinfo=seqinfo)
  names(gr) <- mapping$QNAME
  gr
}
