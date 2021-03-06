#' @export
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicAlignments GAlignments mapToAlignments
#' @importFrom stringr 'str_sub'
#' @importFrom IRanges subsetByOverlaps CharacterList
#' @importFrom BiocGenerics width unlist start end grepl which
#' @importFrom gtools permutations
#' @importFrom S4Vectors 'elementMetadata<-' elementMetadata split
#' @importFrom GenomicRanges GRangesList
#' @importFrom data.table rbindlist

line1_inference <- function(x, n_reads = NULL, search_range = 10000){
  # This functions infers the line1 insertion from the annotated constructed reads.
  storage <- list("genomic_qname_5p" = NA,
                  "genomic_rname_5p" = NA,
                  "genomic_start_5p" = NA,
                  "genomic_end_5p" = NA,
                  "genomic_orientation_5p" = NA,
                  "genomic_break_5p" = NA,
                  "insert_qname_5p" = NA,
                  "insert_rname_5p" = NA,
                  "insert_start_5p" = NA,
                  "insert_end_5p" = NA,
                  "insert_orientation_5p" = NA,
                  "transduced_genomic_region_5p" = NA,
                  "insert_break_5p" = NA,
                  "duplicated_seq_5p" = NA,
                  "is_exact_5p" = NA,
                  "n_reads_5p" = NA,
                  "genomic_qname_3p" = NA,
                  "genomic_rname_3p" = NA,
                  "genomic_start_3p" = NA,
                  "genomic_end_3p" = NA,
                  "genomic_orientation_3p" = NA,
                  "genomic_break_3p" = NA,
                  "insert_qname_3p" = NA,
                  "insert_rname_3p" = NA,
                  "insert_start_3p" = NA,
                  "insert_end_3p" = NA,
                  "insert_orientation_3p" = NA,
                  "insert_break_3p" = NA,
                  "transduced_genomic_region_3p" = NA,
                  "duplicated_seq_3p" = NA,
                  "is_exact_3p" = NA,
                  "n_reads_3p" = NA,
                  "has_polyA" = NA,
                  "overlapping_genomic_region" = NA)

  # Preprocessing and reorganizing the data
  n_anno <- sapply(x, nrow)
  if(!is.null(n_reads)){
    n_reads <- rep(n_reads, sum(n_anno))
  }else{
    n_reads <- rep(NA, sum(n_anno))
  }
  origins <- factor(rep(seq_along(x), n_anno), levels = seq_along(x))
  x <- rbindlist(x)
  x$origin <- origins
  x$n_reads <- n_reads
  y <- x$annotation
  is_anno <- BiocGenerics::grepl(":", y)
  is_a <- x$has_polyA
  is_t <- x$has_polyT
  is_a_seq <- is_a & !is_anno
  is_t_seq <- is_t & !is_anno
  y[is_a_seq] <- convert_polyAT_to_character(y[is_a_seq], "+")
  y[is_t_seq] <- convert_polyAT_to_character(y[is_t_seq], "-")
  x$processed <- y
  is_anno <- BiocGenerics::grepl(":", y)
  y <- x[is_anno]
  extra_info <- list(QNAME = y$QNAME,origin = y$origin, is_polyA = y$has_polyA|y$has_polyT, n_reads = y$n_reads)
  gr <- convert_character2gr(y$processed, extra_info = extra_info)
  polyA_gr <- gr[gr$is_polyA]
  polyA_gr <- split(polyA_gr, polyA_gr$origin)
  # polyA_strands <- split(strand(polyA_gr), polyA_gr$origin)
  # polyA_strands[lengths(polyA_strands) != 1] <- RleList(factor("*", levels = c("+", "-", "*")))
  # polyA_strands <- as.character(polyA_strands)
  gr <- split(gr, gr$origin)
  elementMetadata(gr)$read_annotation <- split(x, x$origin)

  # Reverse complement the whole reads if there is no positive strand in the constructed read
  no_plus_strand <- !any(strand(non_polyA_or_insert(gr)) == "+")
  non_insert_gr <- non_polyA_or_insert(gr)
  has_polyA_reads <- lengths(polyA_gr) > 0
  needed_to_be_rc <- no_plus_strand
  gr[needed_to_be_rc] <- .rc_granges(gr[needed_to_be_rc])
  # polyA_strands[polyA_strands == "+" & needed_to_be_rc] <- "-"
  # polyA_strands[polyA_strands == "-" & needed_to_be_rc] <- "+"
  # assigned_direction <- rep("undetermined", length(gr))
  # assigned_direction[polyA_strands == "-"] <- "left"
  # assigned_direction[polyA_strands == "+"] <- "right"

  # Checking the presence of LINE1 insertion
  ins_gr <- gr[seqnames(gr) == "Hot_L1_polyA"]
  no_read_has_polyA <- sum(has_polyA_reads) == 0
  is_insert_body <- seqnames(gr) == "Hot_L1_polyA" & start(gr) < 6000
  has_insert_reads <- sum(sum(is_insert_body)) != 0
  used_origins <- 0L
  if(has_insert_reads){
    # If has, getting the one having the minimum LINE1 start
    ins_start <- start(ins_gr)
    has_min_insert <- any(ins_start == min(min(ins_start)))
    potential_start <- gr[has_min_insert]
    selected_ins_gr <- non_polyA_or_insert(potential_start, invert = TRUE)
    ins_width <- width(selected_ins_gr)
    selected_start <- potential_start[any(ins_width == max(unlist(ins_width)))][[1]]
    start_is_polyA <- seqnames(selected_start) == "polyA"
    start_is_polyA_only <- seqnames(selected_start) == "Hot_L1_polyA" & start(selected_start) > 6000
    selected_start <- selected_start[!(start_is_polyA|start_is_polyA_only)]
    used_origins <- c(used_origins, selected_start$origin)
  }else{
    if(!no_read_has_polyA){
      # If the reads group has no insert but polyA and having two genomic regions,
      # they are potentially transduced regions.
      potential_start <- gr[has_polyA_reads]
      selected_start <- .get_the_best(potential_start)[[1]] # This is the selected end
      n_region <- length(csaw::mergeWindows(selected_start[!selected_start$is_polyA], tol = search_range)$regions)
      if(n_region == 1){
        return(list(GRangesList(), storage, remaining = NULL))
      }
    }else{
      return(list(GRangesList(), storage, remaining = NULL))
    }
  }

  # Determining which way they locate (right = 5p, left = 3p)
  non_insert <- as.logical(non_polyA_or_insert(selected_start, return_TF = TRUE))
  insert_in_start <- selected_start[!non_insert]
  start_left <- non_insert[1]
  start_right <- tail(non_insert, n = 1)
  head_tail_is_overlap <- is_head_tail_overlap(selected_start)
  starter <- selected_start[non_insert]
  used_origins <- c(used_origins, selected_start$origin)
  if(length(csaw::mergeWindows(starter, tol = search_range)$regions) > 2 & !head_tail_is_overlap){
    return(list(GRangesList(), storage, remaining = elementMetadata(gr)[["read_annotation"]][-unique(used_origins)]))
  }

  if((start_left + start_right) == 1){
    # If the direction can be clearly determined
    start_direction <- ifelse(start_left, "left", "right")
    selected_start$type <- "undefined"
    if(start_direction == "left"){
      starter <- starter[1]
      selected_start_ins <- tail(insert_in_start, 1)
    }else{
      starter <- tail(starter, 1)
      selected_start_ins <- insert_in_start[1]
    }
    relative_direction <- .overlap_direction(gr, starter, search_range = 10000)

    # The end direction should be on the right hand side if start is on the left, vice verse.
    end_direction <- ifelse(start_direction == "left", "right", "left")
    # is_opposite_to_start <-  assigned_direction == end_direction & relative_direction == end_direction
    is_opposite_to_start <-  relative_direction == end_direction
    possible_end <- gr[is_opposite_to_start & has_polyA_reads & (!seq_along(gr) %in% selected_start$origin)]
    selected_start$type[selected_start == selected_start_ins] <- "insert"
    selected_start[selected_start == starter]$type <- "genomic"
    storage <- data_input(storage = storage, selected_start, direction = start_direction)

    if(length(possible_end) == 0){
      return(list(GRangesList(gr[selected_start$origin[1]]), storage, remaining = elementMetadata(gr)[["read_annotation"]][-unique(used_origins)]))
    }


    # Determining the paired end
    selected_end <- .get_the_best(possible_end)[[1]]
    selected_end$type <- "undefined"
    if(end_direction == "left"){
      selected_end[1]$type <- "genomic"
      which_end_insert <- min(which(non_polyA_or_insert(selected_end, return_TF = TRUE, invert = TRUE)))
    }else{
      selected_end[length(selected_end)]$type <- "genomic"
      which_end_insert <- max(which(non_polyA_or_insert(selected_end, return_TF = TRUE, invert = TRUE)))
    }
    selected_end[which_end_insert]$type <- "insert"
    storage <- data_input(storage = storage, selected_end, direction = end_direction)
    used_origins <- c(used_origins, selected_end$origin)
    end_not_insert <- non_polyA_or_insert(selected_end, return_TF = TRUE)
    current_transduced_region <- selected_end[ifelse(end_direction == "left",
                                                     max(which(end_not_insert)),
                                                     min(which(end_not_insert)))]
  }
  skip_transduction <- FALSE

  if(sum(start_left + start_right) == 2){
    # If the constructed reads do not have a cleared direction
    if(sum(has_polyA_reads) > 0){
      if(head_tail_is_overlap){
        # If head and tail is overlap, a constructed read contains a full insertion of LINE1
        orphan_result <- handle_orphan(selected_start, storage = storage)
        storage <- orphan_result$storage
        result <- gr[orphan_result$selected_gr]
        used_origins <- c(used_origins, orphan_result$selected_gr)
        skip_transduction <- TRUE # It will skip the transduction inference, sine the insertion was fully constructed.
      }else{
        # Generate all combination and choose the one makes sense
        # When the donor LINE1 contains transduction, then further jump to the other regions and result as orphan transudction
        # This is difficult to be determined. Simply generating all the combination and then chose the one make sense will be easier.
        perm <- permutations(n = length(gr), r = 2)
        perm <- perm[apply(perm ,1 , function(x)any(x %in% which(has_polyA_reads))), ]
        gr_perm <- apply(perm, 1, function(x)c(gr[[x[1]]], gr[[x[2]]]))
        is_ht_ol <- sapply(gr_perm, is_head_tail_overlap)
        gr_perm <- gr_perm[is_ht_ol]
        gr_perm <- gr_perm[sapply(gr_perm, function(x)x_is_on_the_left(x[1], tail(x, 1)))]
        gr_perm <- GRangesList(gr_perm)
        selected <- .get_the_best(gr_perm)
        orphan_result <- handle_orphan(selected[[1]], storage = storage)
        selected_gr <- gr[orphan_result$selected_gr]
        used_origins <- c(used_origins, orphan_result$selected_gr)
        storage <- orphan_result$storage
        end_direction <- orphan_result$end_direction
        selected_start <- if(end_direction=="left"){selected_gr[[2]]}else{selected_gr[[1]]}
        starter <- selected_start[1]
        selected_end <- if(end_direction=="left"){selected_gr[[1]]}else{selected_gr[[2]]}
        end_not_insert <- non_polyA_or_insert(selected_end, return_TF = TRUE)
        current_transduced_region <- selected_end[ifelse(end_direction == "left",
                                                         max(which(end_not_insert)),
                                                         min(which(end_not_insert)))]
      }
    }else{
      return(list(GRangesList(), storage, remaining = elementMetadata(gr)[["read_annotation"]][-unique(used_origins)]))
    }
  }

  if((start_left + start_right) == 0){
    return(list(GRangesList(), storage, remaining = elementMetadata(gr)[["read_annotation"]][-unique(used_origins)]))
  }



  if(!skip_transduction){
    # Getting reads (things) that can be potential transduced regions
    collected_impossible <- c(starter)
    thing2search <- gr[!(names(gr) %in% c(selected_start$origin, selected_end$origin))]
    current_thing2search <- subsetByOverlaps(thing2search, current_transduced_region + search_range, ignore.strand = TRUE)
    current_thing2search <- current_thing2search[sapply(non_polyA_or_insert(current_thing2search), has_more_than_one)]
    in_the_middle <- GRangesList()

    while(length(current_thing2search) > 0){
      # It recursively identifies transduced regions based on the selected end.
      combined <- unlist(current_thing2search, use.names = FALSE)
      ol_regions <- subsetByOverlaps(combined, current_transduced_region + search_range, ignore.strand = TRUE)
      ol_regions <- subsetByOverlaps(ol_regions, collected_impossible, invert = TRUE, ignore.strand = TRUE)
      if(length(ol_regions) == 0){
        break()
      }
      ol_regions <- split(ol_regions, as.character(ol_regions$origin))
      what_strands <- CharacterList(unique(strand(ol_regions)))
      is_one_strand <- lengths(what_strands) == 1
      ol_regions <- ol_regions[is_one_strand]
      current_thing2search <- current_thing2search[is_one_strand]
      what_strands <- as.character(what_strands[is_one_strand])
      strand_not_match <- what_strands != as.character(strand(current_transduced_region))
      current_thing2search[strand_not_match] <- .rc_granges(current_thing2search[strand_not_match])
      thing_direction <- .overlap_direction(current_thing2search, current_transduced_region, search_range = search_range)
      current_thing2search <- current_thing2search[thing_direction == end_direction]
      collected_impossible <- c(collected_impossible, current_transduced_region)

      if(length(current_thing2search) == 0){
        break
      }

      if(end_direction == "left"){
        in_the_middle <- c(in_the_middle, .get_the_best(current_thing2search))
      }else{
        in_the_middle <- c(.get_the_best(current_thing2search), in_the_middle)
      }
      current_transduced_region <- subsetByOverlaps(current_thing2search[[1]], current_transduced_region + search_range, invert = TRUE)
      current_transduced_region <- non_polyA_or_insert(current_transduced_region)
      current_thing2search <- subsetByOverlaps(thing2search, current_transduced_region + search_range, ignore.strand = TRUE)
      thing_origins <- unlist(current_thing2search)$origin
      used_origins <- c(used_origins, unlist(in_the_middle)$origin)
      if(all(thing_origins %in% used_origins) | length(current_thing2search) == 0){
        break
      }
    }

    transduced_regions <- subsetByOverlaps(unlist(c(in_the_middle, selected_end), use.names = FALSE), starter + search_range, invert = TRUE)
    transduced_regions <- non_polyA_or_insert(transduced_regions)
    transduced_regions <- as.character(range(transduced_regions))
    if(!length(transduced_regions) == 0){
      transduced_regions <- paste(transduced_regions, collapse = ",")
    }else{
      transduced_regions <- NA
    }

    if(end_direction == "left"){
      storage[["transduced_genomic_region_5p"]] <- transduced_regions
      result <- c(gr[unique(selected_end$origin)], in_the_middle, gr[unique(selected_start$origin)])
    }else{
      storage[["transduced_genomic_region_3p"]] <- transduced_regions
      result <- c(gr[unique(selected_start$origin)], in_the_middle, gr[unique(selected_end$origin)])
    }
  }
  storage <- .get_TSD(result, storage = storage) # getting the TSD
  return(list(result, storage, remaining = elementMetadata(gr)[["read_annotation"]][-unique(used_origins)]))
}


#' @export
#' @importFrom S4Vectors split
#' @importFrom IRanges CharacterList
#' @importFrom GenomicRanges GRangesList
#' @importFrom S4Vectors 'elementMetadata<-' elementMetadata split revElements
.rc_granges <- function(x){
  # subsetByOverlaps(usable_clusters, GRanges(8, IRanges(110581038, 110581082)))
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

#' @importFrom S4Vectors split
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
    p$annotation[!is_not_seq] <- string_reverseComplement(p$annotation[!is_not_seq])
    to_be_replaced <- p$annotation[is_not_seq]
    is_minus <- grepl(":-", to_be_replaced)
    is_plus <- grepl(":\\+", to_be_replaced)
    to_be_replaced[is_minus] <- sub(":-", ":+", to_be_replaced[is_minus])
    to_be_replaced[is_plus] <- sub(":\\+", ":-", to_be_replaced[is_plus])
    p$annotation[is_not_seq] <- to_be_replaced
    # RC the sequence
    p$seq <- string_reverseComplement(p$seq)
    return(p)
  })
  return(do.call(rbind, z))
}

data_input <- function(storage, data, direction){
    data_slot <- c("qname", "rname", "start", "end", "orientation")
    data_slot <- paste0(data_slot, ifelse(direction == "left", "_5p", "_3p"))
    n_reads_slot <- paste0("n_reads", ifelse(direction == "left", "_5p", "_3p"))
    for(i in c("genomic", "insert")){
      selected <- data[data$type == i]
      if(length(selected) == 0){
        next
      }
      selected_slot <- paste0(i, "_", data_slot)
      storage[[selected_slot[1]]] <- selected$QNAME
      storage[[selected_slot[2]]] <- as.character(seqnames(selected))
      storage[[selected_slot[3]]] <- start(selected)
      storage[[selected_slot[4]]] <- end(selected)
      storage[[selected_slot[5]]] <- as.character(strand(selected))
      if(i == "genomic"){
        storage[[n_reads_slot]] <- selected$n_reads[1]
      }
    }

  return(storage)
}

#' @importFrom S4Vectors split
#' @importFrom IRanges LogicalList
#' @importFrom GenomicRanges findOverlaps

is_overlap <- function(query, reference){
  ref_is_list <- class(reference) %in% c("list", "CompressedGRangesList")
  if(ref_is_list){
    reference <- unlist(reference)
  }
  is_list <- class(query) %in% c("list", "CompressedGRangesList")
  if(is_list){
    list_grp <- factor(rep(seq_along(query), lengths(query)), levels = seq_along(query))
    query <- unlist(query, use.names = FALSE)
  }
  is_overlap <- seq_along(query) %in% findOverlaps(query, reference, ignore.strand = TRUE)@from

  if(is_list){
    return(LogicalList(split(is_overlap, list_grp)))
  }else{
    return(is_overlap)
  }

}

.overlap_direction <- function(query, seed, search_range){
  if(!class(query) %in% c("GRanges", "CompressedGRangesList")){
    stop()
  }
  is_ol <- is_overlap(query, seed + search_range)
  ol_count <- sum(is_ol)
  is_unwanted <- (ol_count > 1 | ol_count == 0 | lengths(query) == 1)
  is_left <- sapply(is_ol, "[[", i = 1)
  is_right <- sapply(is_ol, tail, n = 1)
  is_unwanted <- is_left + is_right != 1 | is_unwanted
  direction <- ifelse(is_left, "left", "right")
  direction[is_unwanted] <- "undetermined"
  return(direction)
}

.get_the_best <- function(y, search_range = 10000){
  # Gathering information to judge the best guess
  # Priority: Most information > longest informative reads
  n_info <- lengths(lapply(y, function(x)csaw::mergeWindows(x, tol = search_range)$regions))
  most_informative <- n_info == max(n_info)
  longest <- seq_along(y) %in% which.max(sum(width(y)))
  score <- most_informative * 10 + longest
  y[which.max(score)]
}


convert_polyAT_to_character <- function(x, strand = "+"){
  if(strand == "+"){
    return(paste0("polyA", ":", 1, "-", nchar(x), ":",strand = "+"))
  }
  if(strand == "-"){
    return(paste0("polyA", ":", 1, "-", nchar(x), ":",strand = "-"))
  }
}

#' @importFrom stringr str_sub
#' @importFrom IRanges LogicalList
#' @importFrom GenomicAlignments GAlignments mapToAlignments
.get_TSD <- function(result, storage){
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
      storage[["duplicated_seq_5p"]] <- genomic_overlap_seq[1]
      storage[["duplicated_seq_3p"]] <- genomic_overlap_seq[2]
    }
  }
  return(storage)
}


#' @importFrom IRanges overlapsAny
x_is_on_the_left <- function(x, y, search_range = 10000){
  if(overlapsAny(x, y + search_range)){
    x_is_on_the_left <- start(x) < start(y)
    return(x_is_on_the_left)
  }else{
    stop("No overlap is found")
  }
}

has_more_than_one <- function(x, search_range = 10000){
  length(csaw::mergeWindows(x, tol = search_range)$regions) > 1
}

non_polyA_or_insert <- function(x, invert = FALSE, return_TF = FALSE){
  if(invert){
    test <- seqnames(x) %in% c("polyA", "Hot_L1_polyA")
  }else{
    test <- !seqnames(x) %in% c("polyA", "Hot_L1_polyA")
  }
  if(return_TF){
    return(test)
  }else{
    return(x[test])
  }
}

#' @importFrom IRanges overlapsAny
is_head_tail_overlap <- function(x, search_range = 10000){
  head_x <- x[1]
  tail_x <- x[length(x)]
  overlapsAny(head_x, tail_x + search_range)
}

#' @importFrom BiocGenerics unique strand '%in%'
handle_orphan <- function(selected_gr, storage, search_range = 10000){
  selected_5p <- selected_gr[1]
  selected_5p$type <- "genomic"
  storage <- data_input(storage = storage, selected_5p, direction = "left")
  selected_3p <- tail(selected_gr, 1)
  selected_3p$type <- "genomic"
  storage <- data_input(storage = storage, selected_3p, direction = "right")
  insert_strand <- unique(strand(selected_gr[selected_gr$is_polyA]))[1]
  end_direction <- ifelse(as.character(insert_strand[1]) == "+", "right", "left")
  transduced_regions <- selected_gr[-(c(1, length(selected_gr), which(seqnames(selected_gr) %in% c("polyA", "Hot_L1_polyA"))))]
  transduced_regions <- csaw::mergeWindows(transduced_regions, tol = search_range)$regions
  if(length(transduced_regions) > 0){
    transduced_regions <- paste(as.character(transduced_regions), collapse = ",")
  }else{
    transduced_regions <- NULL
  }

  if(end_direction == "right"){
    storage[["transduced_genomic_region_3p"]] <- transduced_regions
  }else{
    storage[["transduced_genomic_region_5p"]] <- transduced_regions
  }
  insert_in_gr <- non_polyA_or_insert(selected_gr, invert = TRUE)
  end_1 <- insert_in_gr[1]
  end_1$type <- "insert"
  storage <- data_input(storage = storage, end_1, direction = ifelse(end_direction=="right", "left", "right"))
  end_2 <- tail(insert_in_gr, 1)
  end_2$type <- "insert"
  storage <- data_input(storage = storage, end_2, direction = ifelse(end_direction=="right", "right", "left"))
  return(list(selected_gr = unique(selected_gr$origin), storage = storage, end_direction = end_direction))
}
