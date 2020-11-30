#' @export
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom data.table fread
#' @importFrom S4Vectors 'elementMetadata<-'
#' @importFrom IRanges start end CharacterList IntegerList
#' @importFrom BiocParallel bplapply MulticoreParam

import_files <- function(mm, disc, include_unmapped = TRUE,
                         disc_min_mapq = 20, anchor_min_mapq = 5,
                         BPPARAM = MulticoreParam(workers = 10L)){
  # Total multiple mapped reads
  mm_gr <- .internal_import(x = mm, BPPARAM = BPPARAM)

  # Determine split reads (multi and unique mapped reads are on the same reads)
  mm_gr$is_anchor <- mm_gr$MAPQ >= anchor_min_mapq & !mm_gr$is_supp
  mm_gr$is_disc <- FALSE

  # For discordant reads
  disc_gr <- .internal_import(x = disc, mapq_filter = disc_min_mapq, BPPARAM = BPPARAM)
  seqname <- CharacterList(split(as.character(seqnames(disc_gr)), disc_gr$QNAME))
  multiple_seqnames <- lengths(unique(seqname)) > 1
  start_range <- IntegerList(split(start(disc_gr), disc_gr$QNAME))
  start_range <- range(start_range)
  keep <- rownames(start_range)[(start_range[,2] - start_range[,1] > 5000)|multiple_seqnames] # Extracting read that are far apart enough
  disc_gr <- disc_gr[disc_gr$QNAME%in%keep]
  disc_gr$is_anchor <- FALSE
  disc_gr$is_disc <- TRUE
  out <- c(mm_gr, disc_gr)
  names(out) <- out$QNAME_id
  out
}

.internal_import <- function(x, mapq_filter = NULL,
                             include_unmapped = TRUE,
                             BPPARAM = BPPARAM){
  what <- c("character","integer","character", "integer", "integer", "character", "character")
  x <- do.call(rbind, bplapply(x, fread, colClasses = what, BPPARAM = BPPARAM))
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
