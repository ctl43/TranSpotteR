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