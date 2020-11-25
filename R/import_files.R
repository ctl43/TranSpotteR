#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom data.table fread
#' @importFrom S4Vectors 'elementMetadata<-'

import_files <- function(x, mapq_filter=NULL, include_unmapped=TRUE){
  what <- c("character","integer","character", "integer", "integer", "character","character")
  x <- do.call(rbind, lapply(x, fread, colClasses=what))
  x$is_first <- !!bitwAnd(x$FLAG, 0x40)
  x$is_rc <- !!bitwAnd(x$FLAG, 0x10)
  x$is_supp <- !!bitwAnd(x$FLAG, 0x800)
  x$QNAME_id <- paste0(x$QNAME,"_",as.integer(!x$is_first)+1)

  # Converting all hard clipped to soft clipped
  x <- convertingHtoS(x)
  hs37d5_seq_info <- readRDS("/home/ctlaw/dicky/reference/hs37d5_seq_info.rds")
  if(include_unmapped){
    hs37d5_seq_info@seqnames <- c(hs37d5_seq_info@seqnames, "UNMAPPED")
    hs37d5_seq_info@seqlengths <- c(hs37d5_seq_info@seqlengths, 99999L)
    hs37d5_seq_info@is_circular <- c(hs37d5_seq_info@is_circular, FALSE)
    hs37d5_seq_info@genome <- c(hs37d5_seq_info@genome, "hs37d5_KJ173426")
    x[x$CIGAR=="*", ]$RNAME <- "UNMAPPED"
    x[x$CIGAR=="*", ]$POS <- 1
    x[x$CIGAR=="*", ]$CIGAR <- paste0(nchar(x[x$CIGAR=="*", ]$SEQUENCE),"M")
  }
  gr <- sam2gr(x, seqinfo=hs37d5_seq_info)
  elementMetadata(gr) <- x

  if(!is.null(mapq_filter)){
    discard <- gr$QNAME[gr$MAPQ<mapq_filter]
    gr <- gr[!gr$QNAME%in%discard, ]
    gr
  }else{
    gr
  }
}
