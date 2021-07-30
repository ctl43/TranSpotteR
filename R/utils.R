#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors split
#' @importFrom IRanges IRanges
#' @importFrom IRanges CharacterList

# Converting read annotation to genomic range
convert_character2gr <- function(y, extra_info = NULL){
  # Converting character to a GenomicRange object
  # X:1-100:+ = GRanges("X", IRanges(1, 100), strand = "+")
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
  y <- strsplit(y, ":")
  y <- do.call(rbind, y)
  chr <- y[,1]
  strand <- y[,3]
  y <- IntegerList(strsplit(y[, 2], "-"))
  start <- unlist(first_element(y))
  end <- unlist(last_element(y))
  gr <- GRanges(seqnames=chr, IRanges(start = start, end = end), strand = strand)

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
  # Converting the rle object into CIGAR string
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
#' @importFrom S4Vectors runValue runLength aggregate revElements
#' @importFrom GenomicAlignments cigarToRleList

unify_cigar_strand <- function(cigar, flag = NULL, from, to, along_query = FALSE){
  # Operating the strand of CIGAR strand
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
  # Converting a sam derived table to a GenomicRange object
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
#' @importFrom Biostrings reverseComplement DNAStringSet
string_reverseComplement <- function(strings){
  as.character(reverseComplement(DNAStringSet(strings)))
}

#' @export
#' @importFrom Biostrings reverseComplement DNAStringSet
convertingHtoS <- function(x, unique_id="QNAME_id"){
  # Converting the hard clipped reads into soft clipped reads and adding back the hard clipped sequence
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

#' @export
#' @importFrom Biostrings DNAStringSet
#' @importFrom ShortRead writeFasta

bwa_alignment <- function(seq, ref="/home/ctlaw/reference/Homo_sapiens/hs37d5/hs37d5_KJ173426.fa",
                          working_dir=NULL, samtools_param="-F 4", call_bwa = "bwa mem ", threads = 5){
  # A warpper for calling the bwa in R
  if(length(seq)==0){
    return(data.frame("QNAME"= character(), "FLAG"= integer(),
                      "RNAME"= character(), "POS"= integer(),
                      "MAPQ"= integer(), "CIGAR"= character(),
                      "SEQUENCE"= character()))
  }
  read_names <- names(seq)
  if(is.null(read_names)){
    names(seq) <- paste0("read_", seq_along(seq))
  }

  if(any(duplicated(read_names))){
    stop("Duplicated sequence name is detected.")
  }

  if (is.null(working_dir)) {
    tempdir()
    working_dir <- tempfile()
    dir.create(working_dir)
    on.exit(unlink(working_dir, recursive = TRUE))
  }
  if(samtools_param==""){
    call_samtools <- "|samtools view /dev/stdin"
  }else{
    call_samtools <- paste0("|samtools view ", samtools_param, " /dev/stdin")
  }
  temp_fa <- file.path(working_dir, basename(working_dir))
  writeFasta(DNAStringSet(seq), temp_fa)
  cmd <- paste0(call_bwa, "-t ",threads, " ",ref, " ", temp_fa, call_samtools,"|cut -f 1,2,3,4,5,6,10")
  result <- system(cmd, intern = TRUE)
  if(length(result)==0){
    return(data.frame(
      "QNAME"= character(), "FLAG" = integer(),
      "RNAME"= character(), "POS" = integer(),
      "MAPQ"= integer(), "CIGAR" = character(),
      "SEQUENCE"= character()))
  }else{
    bwa_aln <- read.delim(text = result, stringsAsFactors = FALSE, header = FALSE)
    colnames(bwa_aln) <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "SEQUENCE")
    bwa_aln$QNAME <- as.character(bwa_aln$QNAME)
    bwa_aln
  }
}

#' @export
#' @importFrom S4Vectors split
#' @importFrom IRanges CharacterList IntegerList
#' @importFrom GenomicAlignments cigarToRleList explodeCigarOpLengths
#' @importFrom Biostrings DNAStringSet reverseComplement DNAStringSetList

get_unmapped_clipped_read <- function(x, include_middle_unmapped = TRUE){
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

#' @export
first_element <- function(x, invert = FALSE){
  if(length(x) == 0){
    return(x)
  }
  grp <- factor(rep(seq_along(x), lengths(x)), levels = seq_along(x))
  flat <- unlist(x, use.names = FALSE)
  is_dup <- duplicated(grp)
  if(invert){
    return(split(flat[is_dup], grp[is_dup]))
  }else{
    return(split(flat[!is_dup], grp[!is_dup]))
  }
}

#' @export
last_element <- function(x, invert = FALSE){
  if(length(x) == 0){
    return(x)
  }
  grp <- factor(rep(seq_along(x), lengths(x)), levels = seq_along(x))
  flat <- unlist(x, use.names = FALSE)
  rev_grp <- rev(grp)
  rev_flat <- rev(flat)
  is_dup <- duplicated(rev_grp)
  if(invert){
    return(split(rev(rev_flat[is_dup]), rev(rev_grp[is_dup])))
  }else{
    return(split(rev(rev_flat[!is_dup]), rev(rev_grp[!is_dup])))
  }
}

#' @export
#' @importFrom BiocGenerics unstrand start end
#' @importFrom IRanges LogicalList IntegerList
how_many_regions_in_range <- function(x, tol = 10000, ignore_strand = TRUE){
  lengths(merge_granges(x, tol)[[2]], ignore_strand = ignore_strand)
}

#' @export
#' @importFrom S4Vectors split
merge_granges <- function(x, tol, ignore_strand = FALSE){
  is_grlist <- class(x)=="CompressedGRangesList"
  if(is_grlist){
    grp <- get_list_grp(x, as_factor = TRUE)
    gr <- unlist(x)
  }else{
    gr <- x
    grp <- rep(1, length(gr), as_factor = TRUE)
  }
  chrom <- as.integer(seqnames(gr))
  start <- as.integer(start(gr))
  end <- as.integer(end(gr))

  if(ignore_strand){
    strand <- rep(1, length(gr), as_factor = TRUE)
    gr <- unstrand(gr)
  }else{
    strand <- as.integer(factor(strand(gr), levels = c("+", "-", "*")))
  }

  o <- order(grp, strand, chrom, start, end) # MUST BE SORTED LIKE THIS
  idx <- rep(0, length(gr))
  idx[o] <- cxx_merge_ranges(chrom[o], start[o], end[o], grp = grp[o], strand = strand[o], tol = tol)
  y <- IntegerList(split(idx, grp))
  members <- split(gr, unlist(idx))
  megred <- unlist(range(members))
  out_grp <- y - unname(cumsum(c(0, head(lengths(unique(y)), -1)))) # converting the idx to each merged
  out_grp <- IntegerList(out_grp)

  if(is_grlist){
    merged_grp <- unique(y)
    out_merged <- split(megred, get_list_grp(merged_grp, as_factor = TRUE))
    names(out_grp) <- names(out_merged) <- names(x)
    return(list(idx = out_grp, regions = out_merged))
  }else{
    return(list(idx = unlist(out_grp, use.names = FALSE), regions = megred))
  }
}

#' @export
#' @importFrom BiocGenerics paste
#' @importFrom stringr str_count
grl_to_character <- function(x){
  grp <- factor(rep(seq_along(x), lengths(x)), levels = seq_along(x))
  x <- unlist(x)
  x <- as.character(x)
  x[str_count(x, ":") == 1] <- paste0(x[str_count(x, ":") == 1], ":*")
  out <- paste(CharacterList(split(x, grp)), collapse = ",")
  out[out == ""] <- NA
  return(out)
}

#' @export
extract_element_at <- function (x, at, invert = FALSE){
  grp <- factor(rep(seq_along(x), lengths(x)), levels = seq_along(x))
  pos <- unlist(cumsum(IntegerList(split(duplicated(grp), grp))) + 1)
  flat <- unlist(x, use.names = FALSE)
  selected <- flat[pos == 2]
  selected_grp <- grp[pos == 2]
  if(invert == TRUE){
    split(flat[pos != at], grp[pos != at])
  }else{
    split(flat[pos == at], grp[pos == at])
  }
}
