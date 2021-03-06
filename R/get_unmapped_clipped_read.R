#' @export
#' @importFrom S4Vectors split
#' @importFrom IRanges CharacterList IntegerList
#' @importFrom GenomicAlignments cigarToRleList explodeCigarOpLengths
#' @importFrom Biostrings DNAStringSet reverseComplement DNAStringSetList

get_unmapped_clipped_read <- function(x, include_middle_unmapped=TRUE){
  # This function extracts the unmapped, clipped sequence
  # and the umapped sequence that is in the middle of a read.

  x$QNAME <- as.character(x$QNAME)

  # Initialisation
  unmapped_result <- r_result <- l_result <- rep("", sum(!duplicated(x$QNAME)))
  m_result <- rep(CharacterList(character()), sum(!duplicated(x$QNAME)))
  m_end <-  m_start <- rep(IntegerList(0), sum(!duplicated(x$QNAME)))
  names(m_end) <-  names(m_start) <- names(m_result) <- names(unmapped_result) <- names(r_result) <- names(l_result) <- unique(x$QNAME)

  # Getting unmapped reads
  is_unmapped <- x$CIGAR=="*"
  unmapped <- x[is_unmapped,]
  unmapped_result[unmapped$QNAME] <- unmapped$SEQUENCE

  # Processing mapped reads and getting their clipped sequences
  x <- x[!is_unmapped,]
  group <- x$QNAME
  cigar <- x$CIGAR
  seq <- x$SEQUENCE
  flag <- x$FLAG
  names(seq) <- group

  if(length(seq)==0){
    result <- list(unmapped = DNAStringSet(unmapped_result), left = DNAStringSet(l_result), right = DNAStringSet(r_result))
    list(clipped = result, middle = list(start = m_start,end = m_end, seq = m_result))
  }

  # Making sure the sequence are the same length as the cigar string
  len <- cigarWidthAlongQuerySpace(cigar)
  if(!all(len == nchar(seq))){
    stop()
  }

  # Unifying the strand of sequences
  strand <- ifelse(bitwAnd(flag, 0x10), "-", "+")
  is_rc <- strand == "-"
  seq[is_rc] <- reverseComplement(DNAStringSet(seq[is_rc]))

  # Unifying the strand of CIGAR string
  unified_cigar <- unify_cigar_strand(cigar, from = strand, to = "+")
  unified_cigar <- gsub("[0-9]+D", "", unified_cigar)

  # Get clipped length
  left_len <- .get_clip_length(unified_cigar, start = TRUE)
  right_len <- .get_clip_length(unified_cigar, start = FALSE)
  left_len <- min(IntegerList(split(left_len, group)))
  right_len <- min(IntegerList(split(right_len, group)))
  seq <- unlist(unique(DNAStringSetList(split(seq, group))))
  len <- max(IntegerList(split(len, group)))
  l_clipped <- substr(seq, 1, left_len)
  r_clipped <- substr(seq, len-right_len + 1, len)
  l_result[names(l_clipped)] <- l_clipped
  r_result[names(r_clipped)] <- r_clipped
  result <- list(unmapped = DNAStringSet(unmapped_result), left = DNAStringSet(l_result), right = DNAStringSet(r_result))

  if(include_middle_unmapped){
    # The middle unmapped regions are defined after combining all mapping
    cigar_rle <- as.list(cigarToRleList(unified_cigar))
    cigar_rle <- base::split(cigar_rle, group)
    cigar_rle <- lapply(cigar_rle, .combined_rle)
    combined_cigar <- sapply(cigar_rle, rle_to_cigar)
    middle_clipped_pos <- .get_middle_unmapped_pos(combined_cigar)
    m_start_tmp <- middle_clipped_pos$start
    m_end_tmp <- middle_clipped_pos$end
    m_seq <- rep(seq, lengths(m_start_tmp))
    m_grp <- rep(names(combined_cigar), lengths(m_start_tmp))
    m_clipped <- substr(m_seq, unlist(m_start_tmp), unlist(m_end_tmp))
    m_clipped <- CharacterList(split(unname(m_clipped), m_grp))
    m_clipped[any(m_clipped == "")] <- CharacterList(character(0))
    m_clipped <- lapply(m_clipped, function(x){names(x) <- seq_along(x); x})
    m_result[names(m_clipped)] <- CharacterList(m_clipped)
    m_start[names(m_start_tmp)] <- IntegerList(m_start_tmp)
    m_end[names(m_end_tmp)] <- IntegerList(m_end_tmp)
    result <- list(clipped = result, middle = list(start=m_start,end=m_end, seq=m_result))
  }
  result
}

#' @export
.get_clip_length <- function(cigars, start = TRUE) {
  cliplen <- integer(length(cigars))
  for (op in c("H", "S")) { # hard clips before soft clips.
    if (start) {
      finder <- paste0("^[0-9]+", op)
      keeper <- paste0("^([0-9]+)", op, ".*")
    } else {
      finder <- paste0("[0-9]+", op, "$")
      keeper <- paste0(".*[^0-9]([0-9]+)", op, "$")
    }
    present <- grepl(finder, cigars)
    cliplen[present] <- cliplen[present] + as.integer(sub(keeper, "\\1", cigars[present]))
    cigars <- sub(finder, "", cigars)
  }
  return(cliplen)
}

#' @export
#' @importFrom IRanges IntegerList
.get_middle_unmapped_pos <- function(cigars){
  end <- start <- rep(IntegerList(0), length(cigars))
  present <- grepl("([0-9]+[A-Z][0-9]+S[0-9]+[A-Z])", cigars)

  .internal <- function(cigar){
    storage <- c()
    current <- sub("[0-9]+[A-Z]$", "", cigar)
    last_one <- substr(current, nchar(current), nchar(current))
    continue <- grepl("^([0-9]+[A-Z][0-9]+[A-Z])", current)
    while(continue){
      if(last_one == "S"){
        storage <- c(storage, current)
      }
      current <- sub("[0-9]+[A-Z]$", "", current)
      last_one <- substr(current, nchar(current), nchar(current))
      continue <- grepl("^([0-9]+[A-Z][0-9]+[A-Z])", current)
    }
    out <- rev(IntegerList(strsplit(storage, "[A-Z]")))
    out <- t(sapply(cumsum(out), tail, n = 2))
    list(start = out[, 1] + 1, end = out[,2])
  }
  result <- lapply(cigars[present], .internal)
  start[present] <- IntegerList(lapply(result, "[[", i = "start"))
  end[present] <- IntegerList(lapply(result, "[[", i = "end"))
  start <- lapply(start, function(x){ names(x) <- seq_along(x);x })
  end <- lapply(end, function(x){ names(x) <- seq_along(x);x })
  names(end) <- names(start) <- names(cigars)
  list(start = start, end = end)
}

#' @importFrom BiocGenerics cbind do.call
.combined_rle <- function(x){
  if(length(x) == 1){
    return(x[[1]])
  }
  x <- do.call(cbind, x)
  result <- x[, 2]
  for(i in 3 : ncol(x)){
    result[result == "S"] <- x[, i][result == "S"]
  }
  out <- Rle(result, x[, 1])
  out
}
