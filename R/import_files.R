#' @export
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom data.table fread
#' @importFrom S4Vectors 'elementMetadata<-'
#' @importFrom IRanges start end CharacterList IntegerList
#' @importFrom GenomicAlignments cigarWidthAlongReferenceSpace
#' @importFrom data.table rbindlist fread

import_files <- function(extracted, seq_info){
  what <- c("character","integer","character", "integer", "integer",
            "character", "character", "character")
  x <- rbindlist(lapply(x, fread, colClasses = what))
  dt <- .internal_import(extracted)
  # Determine split reads (multi and unique mapped reads are on the same reads)
  dt <- convertingHtoS(dt) # Converting all hard clipped to soft clipped
  dt[dt$CIGAR=="*", ]$RNAME <- "UNMAPPED"
  dt[dt$CIGAR=="*", ]$POS <- 1
  dt[dt$CIGAR=="*", ]$CIGAR <- paste0(nchar(dt$SEQUENCE[dt$CIGAR=="*"]), "M")
  seq_info@seqnames <- c(seq_info@seqnames, "UNMAPPED")
  seq_info@seqlengths <- c(seq_info@seqlengths, 99999L)
  seq_info@is_circular <- c(seq_info@is_circular, FALSE)
  seq_info@genome <- c(seq_info@genome, seq_info@genome[1])
  gr <- GRanges(seqnames = dt$RNAME,
                IRanges(dt$POS, width = cigarWidthAlongReferenceSpace(dt$CIGAR)),
                strand = ifelse(bitwAnd(dt$FLAG, 16), "-", "+"),
                seqinfo = seq_info)
  names(gr) <- dt$QNAME_id
  elementMetadata(gr) <- dt
  gr
}
