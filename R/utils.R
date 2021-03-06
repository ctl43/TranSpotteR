#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors split
#' @importFrom IRanges IRanges
#' @importFrom IRanges CharacterList

# Converting read annotation to genomic range
convert_character2gr <- function(y, extra_info = NULL){

  is_list <- class(y)%in%c("list", "CompressedCharacterList")
  if(is_list){
    list_grp <- factor(rep(seq_along(y), lengths(y)), levels = seq_along(y))
    y <- unlist(y, use.names = FALSE)
  }

  if(length(y) == 0 & is_list){
    return(GRangesList(GRanges()))
  }

  if(length(y) == 0 & !is_list){
    return(GRanges())
  }

  if(any(!grepl(":", y))){
    stop()
  }
  y <- unlist(y)
  y <- do.call(rbind, strsplit(y, ":|-"))
  y[, 4][y[, 4] == ""] <- "-" # adding back the negative strand
  y <- data.frame(y)
  y[, 2] <- as.integer(y[, 2])
  y[, 3] <- as.integer(y[, 3])
  gr <- GRanges(seqnames=y[, 1], IRanges(start = y[, 2], end = y[, 3]), strand = y[, 4])

  if(!is.null(extra_info)){
    extra_info <- lapply(extra_info, unlist)
    for(i in names(extra_info)){
      elementMetadata(gr)[[i]] <- extra_info[[i]]
    }
  }

  if(is_list){
    return(S4Vectors::split(gr, list_grp))
  }else{
    return(gr)
  }
}

#' @export
#' @importFrom S4Vectors runValue runLength
rle_to_cigar <- function(x){
  rl <- runLength(x)
  rv <- runValue(x)
  what <- class(x)
  if(what =="CompressedRleList"){
    grp <- rep(seq_along(rv), lengths(rv))
    converted <- aggregate(paste0(unlist(rl), unlist(rv)), list(grp), paste, collapse="")[,2]
  }

  if(what == "Rle"){
    converted <- paste(paste0(rl, rv), collapse = "")
  }
  return(converted)
}

#' @export
#' @importFrom S4Vectors runValue runLength aggregate
#' @importFrom GenomicAlignments cigarToRleList

unify_cigar_strand <- function(cigar, flag = NULL, from, to, along_query = FALSE){
  out <- rep("", length(cigar))
  is_unmapped <- cigar == "*"
  if(all(is_unmapped)){
    return(cigar)
  }
  cigar <- cigar[!is_unmapped]
  flag <- flag[!is_unmapped]
  if(!is.null(flag)){
    from <- ifelse(bitwAnd(flag, 0x10), "-", "+")
  }
  cigar <- cigarToRleList(cigar)
  is_reversed <- from != to
  cigar[is_reversed] <- revElements(cigar[is_reversed])
  rv <- runValue(cigar)
  rl <- runLength(cigar)
  out[!is_unmapped] <- rle_to_cigar(cigar)
  out[is_unmapped] <- "*"
  out
}


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

#' @export
#' @import Biostrings reverseComplement DNAStringSet
string_reverseComplement <- function(strings){
  as.character(reverseComplement(DNAStringSet(strings)))
}

#' @export
#' @importFrom Biostrings reverseComplement DNAStringSet
convertingHtoS <- function(x, unique_id="QNAME_id"){
  is_rc <- !!bitwAnd(x$FLAG, 0x10)
  is_primary <- !bitwAnd(x$FLAG, 0x100)
  non_supp <- !bitwAnd(x$FLAG, 0x800)
  x$SEQUENCE[x$SEQUENCE=="*"] <- ""
  seq <- x$SEQUENCE
  seq <- DNAStringSet(seq)
  seq[is_rc] <- reverseComplement(seq[is_rc]) # Converting all sequence to + strand
  non_hard_seq <- seq[non_supp&is_primary]
  names(non_hard_seq) <- x[non_supp&is_primary,][[unique_id]]
  tmp_seq <- non_hard_seq[as.character(x[[unique_id]])]
  tmp_seq[is_rc] <- reverseComplement(tmp_seq[is_rc])
  x$SEQUENCE <- as.character(tmp_seq)
  x$CIGAR <- gsub("H","S", x$CIGAR)
  x
}
