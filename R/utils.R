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
#' @importFrom S4Vectors runValue runLength aggregate
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
first_element <- function(x, invert = FALSE){
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
how_many_regions_in_range <- function(x, tol = 10000){
  n <- lengths(x)
  grp <- factor(rep(seq_along(x), n), levels = seq_along(x))
  flat <- unlist(x)
  flat$grp <- grp
  flat <- sort(unstrand(flat))
  no_first_start <- start(flat)[duplicated(flat$grp)]
  no_first_start <- IntegerList(split(no_first_start, flat$grp[duplicated(flat$grp)]))
  no_last_end <- rev(rev(end(flat))[duplicated(rev(flat$grp))])
  no_last_end <- IntegerList(split(no_last_end, rev(rev(flat$grp)[duplicated(rev(flat$grp))])))
  dist_diff <- no_first_start - no_last_end
  seq_is_dup <- duplicated(split(seqnames(flat), flat$grp))
  tmp_grp <- factor(rep(seq_along(seq_is_dup), lengths(seq_is_dup)), levels = seq_along(seq_is_dup))
  tmp_1 <- unlist(seq_is_dup, use.names = FALSE)
  seq_is_dup <- LogicalList(split(tmp_1[duplicated(tmp_grp)], tmp_grp[duplicated(tmp_grp)]))
  dist_is_diff <- abs(dist_diff) > tol
  n_diff_regions <- sum(!((!dist_is_diff) & seq_is_dup)) + 1
  n_diff_regions[lengths(x)==0] <- 0
  n_diff_regions
}


merge_glist <- function(x, tol = 10000){
  n <- lengths(x)
  grp <- factor(rep(seq_along(x), n), levels = seq_along(x))
  flat <- unlist(x)
  strand(flat) <- "*"
  flat$grp <- grp
  flat <- sort(unstrand(flat))
  no_first_start <- start(flat)[duplicated(flat$grp)]
  no_first_start <- IntegerList(split(no_first_start, flat$grp[duplicated(flat$grp)]))
  no_last_end <- rev(rev(end(flat))[duplicated(rev(flat$grp))])
  no_last_end <- IntegerList(split(no_last_end, rev(rev(flat$grp)[duplicated(rev(flat$grp))])))
  dist_diff <- no_first_start - no_last_end
  seq_is_dup <- duplicated(split(seqnames(flat), flat$grp))
  tmp_grp <- factor(rep(seq_along(seq_is_dup), lengths(seq_is_dup)), levels = seq_along(seq_is_dup))
  tmp_1 <- unlist(seq_is_dup, use.names = FALSE)
  seq_is_dup <- LogicalList(split(tmp_1[duplicated(tmp_grp)], tmp_grp[duplicated(tmp_grp)]))
  dist_is_diff <- abs(dist_diff) > tol
  is_ol <- (!dist_is_diff) & seq_is_dup
  if(sum(sum(is_ol)) == 0){
    return(x)
  }
  ol_grp <- rep(seq_along(is_ol), lengths(is_ol))
  is_ol <- LogicalList(split(c(rep(FALSE, length(is_ol)), unlist(is_ol, use.names = FALSE)),
                             c(seq_along(is_ol), ol_grp))) # appending a TRUE at the first position
  is_ol[[which(lengths(x) == 0)]] <- logical(0L)
  merged_grp <- cumsum(!is_ol)
  ordered_flat <- unlist(split(flat, flat$grp))
  merged_grp <- paste0(unlist(merged_grp), "_", ordered_flat$grp)
  final_grp <- unlist(unique(IntegerList(split(as.integer(ordered_flat$grp), merged_grp))), use.names = FALSE)
  final_grp <- factor(final_grp, levels = levels(grp))
  out <- split(unlist(range(split(ordered_flat, merged_grp)), use.names = FALSE), final_grp)
  return(out)
}

grl_to_character <- function(x){
  grp <- factor(rep(seq_along(x), lengths(x)), levels = seq_along(x))
  x <- unlist(x)
  x <- as.character(x)
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


#
# test <- GRangesList(GRanges(c(1, 1, 2, 2), IRanges(c(1, 5, 2, 1), c(10, 20, 1000,50))),GRanges(),
#                     GRanges(c(1, 1, 3,3), IRanges(c(100000, 500000, 1, 1000), c(200000, 600000, 2, 50000))))
# # test <- rep(test, 10)
# merged_out <- merge_glist(test, tol = 1000)
# csaw_out <- lapply(test, function(x)csaw::mergeWindows(x, tol = 1000)$regions)
