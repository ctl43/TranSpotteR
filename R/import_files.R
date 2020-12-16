#' @export
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom data.table fread
#' @importFrom S4Vectors 'elementMetadata<-'
#' @importFrom IRanges start end CharacterList IntegerList
#' @importFrom  GenomicAlignments cigarWidthAlongReferenceSpace

import_files <- function(mm, disc, include_unmapped = TRUE,
                         disc_min_mapq = 20, anchor_min_mapq = 5){
  # Total multiple mapped reads
  mm_dt <- .internal_import(x = mm)
  # Determine split reads (multi and unique mapped reads are on the same reads)
  mm_dt$is_anchor <- mm_dt$MAPQ >= anchor_min_mapq & !mm_dt$is_supp
  mm_dt$is_disc <- FALSE

  # For discordant reads
  disc_dt <- .internal_import(x = disc)
  has_multiple_seqnames <- lengths(unique(CharacterList(split(disc_dt$RNAME, disc_dt$QNAME)))) > 1
  start_range <- range(IntegerList(split(disc_dt$POS, disc_dt$QNAME)))
  keep <- rownames(start_range)[(start_range[,2] - start_range[,1] > 5000) | has_multiple_seqnames] # Extracting read that are far apart enough
  disc_dt <- disc_dt[disc_dt$QNAME %in% keep,]
  disc_dt$is_anchor <- FALSE
  disc_dt$is_disc <- TRUE
  to_keep <- disc_dt$QNAME[disc_dt$MAPQ > disc_min_mapq]
  disc_dt <- disc_dt[disc_dt$QNAME %in% to_keep, ]
  dt <- rbind(mm_dt, disc_dt)
  dt <- convertingHtoS(dt) # Converting all hard clipped to soft clipped
  dt[dt$CIGAR=="*", ]$RNAME <- "UNMAPPED"
  dt[dt$CIGAR=="*", ]$POS <- 1
  dt[dt$CIGAR=="*", ]$CIGAR <- paste0(nchar(dt$SEQUENCE[dt$CIGAR=="*"]),"M")
  seq_info <- readRDS("/home/ctlaw/dicky/reference/hs37d5_seq_info.rds")
  seq_info@seqnames <- c(seq_info@seqnames, "UNMAPPED")
  seq_info@seqlengths <- c(seq_info@seqlengths, 99999L)
  seq_info@is_circular <- c(seq_info@is_circular, FALSE)
  seq_info@genome <- c(seq_info@genome, "hs37d5_KJ173426")
  gr <- GRanges(seqnames = dt$RNAME,
                IRanges(dt$POS, width = cigarWidthAlongReferenceSpace(dt$CIGAR)),
                strand = ifelse(bitwAnd(dt$FLAG, 16), "-", "+"),
                seqinfo = seq_info)
  names(gr) <- dt$QNAME_id
  elementMetadata(gr) <- dt
  gr
}

#' @export
#' @importFrom data.table rbindlist fread
.internal_import <- function(x){
  what <- c("character","integer","character", "integer", "integer", "character", "character")
  x <- rbindlist(lapply(x, fread, colClasses = what))
  x$is_first <- !!bitwAnd(x$FLAG, 0x40)
  x$is_rc <- !!bitwAnd(x$FLAG, 0x10)
  x$is_supp <- !!bitwAnd(x$FLAG, 0x800)
  x$QNAME_id <- paste0(x$QNAME,"_",as.integer(!x$is_first) + 1)
  x
}
