#' @export
#' @importFrom BiocParallel MulticoreParam bplapply bpmapply
#' @importFrom IRanges CharacterList NumericList IntegerList LogicalList
#' @importFrom BiocGenerics width unlist
#' @importFrom S4Vectors List elementMetadata 'elementMetadata<-'
#' @importFrom data.table setDT

line1_inference <- function(clusters, BPPARAM = MulticoreParam(workers = 10L)){
  usable_clusters <- clusters[elementMetadata(clusters)$is_usable]
  if(length(usable_clusters) == 0){
    return(NULL)
  }

  # Tidying up the elementMetadata after groupping
  nreads <- usable_clusters$nreads
  grp <- rep(usable_clusters$group, lengths(usable_clusters$read_annotation))
  anno <- split(unlist(usable_clusters$read_annotation, recursive = FALSE, use.names = FALSE), grp)
  nreads <- split(unlist(usable_clusters$nreads, recursive = FALSE, use.names = FALSE), grp)
  # selected <- 1:100
  # for_parallel_y <- split(anno[selected], as.integer(cut(seq_along(anno)[selected], breaks = BPPARAM$workers)))
  # for_parallel_n_reads <- split(nreads[selected], as.integer(cut(seq_along(nreads)[selected], breaks = BPPARAM$workers)))
  for_parallel_y <- split(anno, as.integer(cut(seq_along(anno), breaks = BPPARAM$workers)))
  for_parallel_n_reads <- split(nreads, as.integer(cut(seq_along(nreads), breaks = BPPARAM$workers)))
  out <- bpmapply(function(a ,b){
    mapply(.internal_inference, y = a,
           n_reads = b,
           SIMPLIFY = FALSE)
    }, a = for_parallel_y, b = for_parallel_n_reads,
    BPPARAM = BPPARAM) # Need a better way to store the number of consisted reads
  out <- unlist(out, recursive = FALSE, use.names = FALSE)
  anno_out <- lapply(out, "[[", i = 1)
  df_out <- lapply(out, "[[", i = 2)
  col_names <- names(df_out[[1]])
  df_out <- data.frame(rbindlist(df_out))
  colnames(df_out) <- col_names
  anno_out <- List(anno_out)
  rownames(df_out) <- names(anno_out) <- names(anno)
  df_out$idx <- rownames(df_out)
  list(anno_out[lengths(anno_out) > 0], df_out[lengths(anno_out) > 0, ])
}


#' @export
#' @importFrom BiocGenerics start end
#' @importFrom GenomicRanges GRangesList
#' @importFrom IRanges subsetByOverlaps CharacterList
#' @importFrom S4Vectors elementMetadata 'elementMetadata<-' revElements
#' @importFrom Biostrings DNAStringSetList letterFrequency complement BString
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicAlignments GAlignments mapToAlignments
#' @importFrom stringr 'str_sub'

.internal_inference <- function(y, n_reads = NULL, search_range = 1000){
  # Initialise a list to collect information
  storage <- list("5p_chr" = NA,
                  "5p_transducted_genomic_region" = NA,
                  "5p_genomic_regions_start" = NA,
                  "5p_genomic_regions_end" = NA,
                  "5p_genomic_break" = NA,
                  "5p_insert_start" = NA,
                  "5p_insert_end" = NA,
                  "5p_insert_Orientation" = NA,
                  "5p_insert_break" = NA,
                  "5p_duplicated_seq" = NA,
                  "5p_is_exact" = NA,
                  "3p_chr" = NA,
                  "3p_genomic_regions_start" = NA,
                  "3p_genomic_regions_end" = NA,
                  "3p_genomic_break" = NA,
                  "3p_insert_start" = NA,
                  "3p_insert_end" = NA,
                  "3p_insert_break" = NA,
                  "3p_insert_Orientation" = NA,
                  "3p_transducted_genomic_region" = NA,
                  "3p_duplicated_seq" = NA,
                  "3p_is_exact" = NA,
                  "has_polyA" = NA,
                  "overlapping_genomic_region" = NA)
  names(y) <- seq_along(y)
  y_copy <- y
  y <- CharacterList(lapply(y, "[[", i = "annotation"))
  is_anno <- BiocGenerics::grepl(":", y)
  y <- y[is_anno]
  x <- convert_character2gr(y)
  elementMetadata(x)$read_annotation <- y_copy

  if(!is.null(n_reads)){
    elementMetadata(x)$n_reads <- n_reads
  }
  names(x) <- seq_along(x)

  # Finding the Granges having insertion
  is_insert <- seqnames(x) == "Hot_L1_polyA"
  ins_gr <- x[is_insert]
  has_insert <- any(is_insert)

  if(sum(has_insert) == 0){
    message("No insert regions is found")
    return(list(list(data.table()), storage))
  }

  # Choosing the one with the earliest start site
  insert_gr <- x[is_insert]
  l1_start <- min(start(insert_gr))
  min_l1 <- min(l1_start)
  is_min_l1 <- l1_start == min_l1

  # If no meaning anchor is found
  if(sum(is_min_l1) == 0 | min_l1 > 6000){
    return(list(list(data.table()), storage))
  }

  if(length(which(is_min_l1)) > 1){
    # Chossing the one with the longest insert as the representative one
    is_min_l1 <- seq_along(x) %in% which(is_min_l1)[which.max(max(width(insert_gr[is_min_l1])))]
    # Or we can choose the one has exact breakpoint later
  }

  # Determining which starter to use
  starter_gr <- x[which(is_min_l1)] # Getting the starting GRange
  non_insert <- seqnames(starter_gr[[1]]) != "Hot_L1_polyA"
  non_insert <- as.logical(non_insert)
  left_match <- non_insert[1]
  right_match <- tail(non_insert, n = 1)

  # Determining the direction of the anchor cluster
  if(any(c(left_match, right_match))){
    start_direction <- ifelse(left_match, "left", "right")
  }else{
    return(list(list(data.table()), storage))
  }

  # Selecting genomic ranges for searching next GRanges
  starter <- starter_gr[[1]][non_insert]

  # For checking whether breakpoints are exact
  start_anno <- elementMetadata(starter_gr)$read_annotation[[1]]

  if(start_direction == "left"){
    starter <- starter[1] # in case it has more than 1 starter
    insert_start <- starter_gr[[1]][!non_insert]
    insert_start <- insert_start[which.min(start(insert_start))] # in case it has two regions mapping to the insert sequence
    storage[["5p_chr"]] <- as.character(seqnames(starter))
    storage[["5p_genomic_regions_start"]] <- start(starter)
    storage[["5p_genomic_regions_end"]] <- end(starter)
    storage[["5p_genomic_break"]] <- end(starter)
    storage[["5p_insert_start"]] <- start(insert_start)
    storage[["5p_insert_end"]] <- end(insert_start)
    storage[["5p_insert_break"]] <- start(insert_start)
    storage[["5p_insert_Orientation"]] <- as.character(strand(insert_start))

    # Checking whether breakpoints are exact
    is_start <- which(start_anno$annotation==as.character(starter))
    is_insert_in_start <- grep("Hot_L1_polyA", start_anno$annotation)
    is_insert_in_start <- max(is_insert_in_start)
    start_is_exact <- start_anno[is_start]$QNAME == start_anno[is_insert_in_start]$QNAME
    storage[["5p_is_exact"]] <- start_is_exact
  }else{
    starter <- tail(starter, n = 1)
    insert_start <- starter_gr[[1]][!non_insert]
    insert_start <- insert_start[which.min(start(insert_start))]
    storage[["3p_chr"]] <- as.character(seqnames(starter))
    storage[["3p_genomic_regions_start"]] <- start(starter)
    storage[["3p_genomic_regions_end"]] <- end(starter)
    storage[["3p_genomic_break"]] <- end(starter)
    storage[["3p_insert_start"]] <- start(insert_start)
    storage[["3p_insert_end"]] <- end(insert_start)
    storage[["3p_insert_break"]] <- start(insert_start)
    storage[["3p_insert_Orientation"]] <- as.character(strand(insert_start))

    # Checking whether breakpoints are exact
    is_start <- which(start_anno$annotation==as.character(starter))
    is_insert_in_start <- grep("Hot_L1_polyA", start_anno$annotation)
    is_insert_in_start <- min(is_insert_in_start)
    start_is_exact <- start_anno[is_start]$QNAME == start_anno[is_insert_in_start]$QNAME
    storage[["3p_is_exact"]] <- start_is_exact
  }

  next_gr <- subsetByOverlaps(x[!is_min_l1], starter + search_range)#, ignore.strand=TRUE)

  # Checking whether they are in the right orientation
  .determine_direction <- function(q, z){
    tot <- seq_along(q)
    ol <- findOverlaps(q, z + search_range, ignore.strand = TRUE)@from
    if(length(ol) > 1){
      return("undetermined")
    }
    is_match <- tot == ol
    right_match <- tail(is_match, 1)
    left_match <- is_match[1]
    end_direction <- "undetermined"
    if(right_match){
      end_direction <- "right"
    }
    if(left_match){
      end_direction <- "left"
    }
    end_direction
  }

  next_direction <- sapply(next_gr, .determine_direction, z = starter)
  is_sensible <- next_direction != start_direction & next_direction != "undetermined"

  if(sum(is_sensible)==0){
    message("no sensible end")
    return(list(y_copy[is_min_l1], storage))
  }

  sensible <- next_gr[is_sensible]
  has_polyA <- LogicalList(lapply(elementMetadata(sensible)$read_annotation, "[[", i = "has_polyA"))
  has_polyA <- any(any(has_polyA))
  storage[["has_polyA"]] <- has_polyA
  sensible <- sensible[has_polyA]

  if(!any(has_polyA)){
    message("no sensible end (no end with polyA)")
    return(list(y_copy[is_min_l1], storage))
  }

  end_direction <- ifelse(start_direction == "left", "right", "left")

  get_the_best_guess <- function(y, end_direction){
    # Gathering information to judge the best guess
    # Priority: Most information > long > closest to the breakpoint
    n_info <- lengths(y) # the number of provided information
    most_informative <- n_info == max(n_info)

    # Determining the one with the closest breakpoint (if they are not long contig)
    if(end_direction == "left"){
      p <- sapply(y, function(x)end(x[1]))
      close_p <- max(p)
    }else{
      p <- sapply(y, function(x)start(x[1]))
      close_p <- min(p)
    }
    is_best_guess <- p == close_p

    bi_score <- as.integer(paste0(as.integer(most_informative), as.integer(elementMetadata(y)$is_long), as.integer(is_best_guess)))
    bi_score
  }

  # Searching regions with correct orientation
  if(length(sensible) > 1){
    score <- get_the_best_guess(sensible, end_direction = end_direction)
    is_ok <- score == max(score)
    end_partner <- sensible[is_ok][1]
  }else{
    end_partner <- sensible
  }


  end_anno <- elementMetadata(end_partner)$read_annotation[[1]]
  is_insert_in_end <- which(end_anno$has_polyA)

  if(end_direction == "left"){
    next_target_gr <- end_partner[[1]][-1]
    ending <- end_partner[[1]][1]
    storage[["5p_chr"]] <- as.character(seqnames(ending))
    storage[["5p_genomic_regions_start"]] <- start(ending)
    storage[["5p_genomic_regions_end"]] <- end(ending)
    storage[["5p_genomic_break"]] <- start(ending)

    # Checking whether the ending is exact
    end_genomic <- as.character(ending)
    is_end <- min(which(end_anno$annotation == end_genomic))
    is_insert_in_end <- min(is_insert_in_end)
    end_is_exact <- end_anno[is_end]$QNAME == end_anno[is_insert_in_end]$QNAME
    storage[["5p_is_exact"]] <- end_is_exact

  }else{

    next_target_gr <- head(end_partner[[1]], -1)
    ending <- tail(end_partner[[1]], n = 1)
    storage[["3p_chr"]] <- as.character(seqnames(ending))
    storage[["3p_genomic_regions_start"]] <- start(ending)
    storage[["3p_genomic_regions_end"]] <- end(ending)
    storage[["3p_genomic_break"]] <- end(ending)

    # Checking whether the ending is exact
    end_genomic <- as.character(ending)
    is_end <- which(end_anno$annotation == end_genomic)
    is_insert_in_end <- max(is_insert_in_end)
    end_is_exact <- end_anno[is_end]$QNAME == end_anno[is_insert_in_end]$QNAME
    storage[["3p_is_exact"]] <- end_is_exact
  }


  # Finding 3' transduction
  # initiation
  next_is_not_insert <- seqnames(next_target_gr) != "Hot_L1_polyA"
  next_target <- next_target_gr[next_is_not_insert]
  middle_regions <- GRangesList()
  transduced_regions <- next_target
  last_time <- rep(FALSE, length(x))
  is_thing4search <- !(names(x) %in% names(next_gr) | is_min_l1)

  while(sum(is_thing4search) > 0){
    thing4search <- x[is_thing4search]
    possible_middle <- subsetByOverlaps(thing4search, next_target + search_range, ignore.strand = TRUE)
    if(length(possible_middle) == 0){
      break
    }
    ol_region <- GRangesList(lapply(possible_middle, function(x)subsetByOverlaps(x, next_target + search_range, ignore.strand = TRUE)))
    strands <- lapply(strand(ol_region), function(x)unique(as.character(x)))
    ol_region <- ol_region[lengths(strands) == 1]
    possible_middle <- possible_middle[lengths(strands) == 1]

    # Flipping the strand if the strand does not match but overlap
    strand_not_match <- as.character(unique(strand(ol_region))) != as.character(strand(next_target))
    possible_middle[strand_not_match] <- .rc_granges(possible_middle[strand_not_match])
    middle_direction <- sapply(possible_middle, .determine_direction, z = next_target)
    sensible_middle <- possible_middle[middle_direction != start_direction & middle_direction != "undetermined"]

    if(length(sensible_middle) == 0){
      break
    }

    # In case the number of sensible middle is more than 1
    best_score <- get_the_best_guess(sensible_middle, end_direction = end_direction)
    middle <- sensible_middle[which.max(best_score)]
    if(end_direction == "left"){
      middle_regions <- c(middle_regions, middle)
    }else{
      middle_regions <- c(middle, middle_regions)
    }
    transduced_regions <- c(transduced_regions, middle[[1]][seqnames(middle[[1]]) != "Hot_L1_polyA"])
    next_target <- subsetByOverlaps(middle[[1]], next_target + search_range, invert = TRUE)
    next_target <- next_target[seqnames(next_target) != "Hot_L1_polyA"]

    if(length(next_target) == 0){
      break
    }

    is_thing4search <- !(elementMetadata(x)$id %in% elementMetadata(sensible_middle)$id | !is_thing4search)

    if(all(last_time == is_thing4search)|all(is_thing4search == FALSE)){
      break
    }
    last_time <- is_thing4search
  }

  # Getting information of the ending TE position
  insert_end <- c(middle_regions, end_partner)
  insert_end <- unlist(insert_end[seqnames(insert_end)=="Hot_L1_polyA"])
  if(length(insert_end) == 0){
    polyA_loc <- tail(which(elementMetadata(end_partner)$read_annotation[[1]]$has_polyA), n = 1)
    polyA_seq <- elementMetadata(end_partner)$read_annotation[[1]]$seq[polyA_loc]
    at_freq <- Biostrings::letterFrequency(BString(polyA_seq), letters = c("A", "T"))
    inferred_strand <- ifelse(at_freq[1] > at_freq[2], "+", "-")
    insert_end <- GRanges("Hot_L1_polyA", IRanges(6022, 6050), strand = inferred_strand)
  }
  insert_end <- insert_end[which.min(start(insert_end))]

  if(start_direction == "left"){
    storage[["3p_insert_start"]] <- start(insert_end)
    storage[["3p_insert_end"]] <- end(insert_end)
    storage[["3p_insert_break"]] <- end(insert_end)
    storage[["3p_insert_Orientation"]] <- as.character(strand(insert_end))
    if(length(transduced_regions)!=0){
      storage[["3p_transducted_genomic_region"]] <- paste(as.character(range(transduced_regions)), collapse = ",")
    }
    result <- c(starter_gr, middle_regions, end_partner)
  }else{
    storage[["5p_insert_start"]] <- start(insert_end)
    storage[["5p_insert_end"]] <- end(insert_end)
    storage[["5p_insert_break"]] <- end(insert_end)
    storage[["5p_insert_Orientation"]] <- as.character(strand(insert_end))
    if(length(transduced_regions)!=0){
      storage[["5p_transducted_genomic_region"]] <- paste(as.character(range(transduced_regions)), collapse = ",")
    }
    result <- c(end_partner, middle_regions, starter_gr)
  }

  # Getting TSD
  read_anno <- elementMetadata(result)$read_annotation
  left_anno <- read_anno[[1]][read_anno[[1]]$cigar!="*"][1]
  right_anno <- tail(read_anno, n = 1)[[1]]
  right_anno <- tail(right_anno[right_anno$cigar != "*", ],n = 1)
  left_gr <- result[[1]][1]
  right_gr <- tail(tail(result, n = 1)[[1]], n = 1)
  lr_intersect <- GenomicRanges::intersect(left_gr, right_gr)
  if(length(lr_intersect) != 0){
    lr_gr <- c(left_gr, right_gr)
    # Translating the overlapping regions to the read position
    alignments <- GAlignments(seqnames(lr_gr),
                              pos = start(lr_gr),
                              cigar = c(left_anno$cigar, right_anno$cigar),
                              strand = strand(lr_gr),
                              names=c("left", "right"))
    mapped_to_reads <- mapToAlignments(lr_intersect, alignments)
    left_ok <- end(mapped_to_reads)[1] == cigarWidthAlongQuerySpace(left_anno$cigar)
    right_ok <- start(mapped_to_reads)[2] == 1
    if(all(c(left_ok, right_ok))){
      genomic_overlap_seq <- str_sub(c(left_anno$seq, right_anno$seq), start = start(mapped_to_reads), end = end(mapped_to_reads))
      storage[["overlapping_genomic_region"]] <- as.character(lr_intersect)
      storage[["5p_duplicated_seq"]] <- genomic_overlap_seq[1]
      storage[["3p_duplicated_seq"]] <- genomic_overlap_seq[2]
    }
  }

  return(list(elementMetadata(result)[[1]], storage))
}

# Require testthat
.rc_granges <- function(x){
  # subsetByOverlaps(usable_clusters, GRanges(8, IRanges(110581038, 110581082)))
  # x <- readRDS("/home/ctlaw/dicky/analysis/Enhancer_hijack/temp/cb1a53da-3b96-11eb-80ad-00d861a997b7.rds")
  if(length(x) == 0){
    return(GRangesList())
  }
  # reverse complement the GRange first
  x <- revElements(x)
  original <- CharacterList(strand(x))
  strand(x)[original=="+"] <- "-"
  strand(x)[original=="-"] <- "+"

  # Reverse complement the information in elementMetadata
  anno <- elementMetadata(x)$read_annotation
  anno <- lapply(anno, .rc_dt)
  elementMetadata(x)$read_annotation <- anno
  return(x)
}

.rc_dt <- function(z){
  z <- split(z, factor(z$QNAME, levels = unique(z$QNAME)))
  z <- lapply(z, function(p){
    # Inverting the read position
    q <- p[, c("end", "start")]
    q <- abs(q - max(q)) + 1
    p[[1]] <- q[[1]]
    p[[2]] <- q[[2]]
    p <- p[order(p[["start"]])]
    is_not_seq <- grepl(":", p$annotation)
    # RC information in read annotation
    p$annotation[!is_not_seq] <- as.character(reverseComplement(DNAStringSet(p$annotation[!is_not_seq])))
    to_be_replaced <- p$annotation[is_not_seq]
    is_minus <- grepl(":-", to_be_replaced)
    is_plus <- grepl(":\\+", to_be_replaced)
    to_be_replaced[is_minus] <- sub(":-", ":+", to_be_replaced[is_minus])
    to_be_replaced[is_plus] <- sub(":\\+", ":-", to_be_replaced[is_plus])
    p$annotation[is_not_seq] <- to_be_replaced
    # RC the sequence
    p$seq <- as.character(reverseComplement(DNAStringSet(p$seq)))
    return(p)
    })
  return(do.call(rbind, z))
}
