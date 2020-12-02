#' @export
#' @importFrom BiocGenerics start end
#' @importFrom GenomicRanges GRangesList
#' @importFrom IRanges subsetByOverlaps CharacterList
#' @importFrom S4Vectors elementMetadata 'elementMetadata<-' revElements
#' @importFrom Biostrings DNAStringSetList letterFrequency complement
#' @importFrom GenomeInfoDb seqnames

insert_inference <- function(x){
  x <- convert_character2gr(x)
  elementMetadata(x)$id <- seq_along(x)

  # Finding the Granges having insertion
  is_insert <- seqnames(x) == "Hot_L1_polyA"
  ins_gr <- x[is_insert]
  has_insert <- any(is_insert)

  if(sum(has_insert) == 0){
    message("No insert regions is found")
    return(GRangesList())
  }

  # Choosing the one with the earliest start site
  insert_gr <- x[is_insert]
  l1_start <- min(start(insert_gr))
  min_l1 <- min(l1_start)
  is_min_l1 <- l1_start == min_l1

  # If no anchor is found
  if(sum(is_min_l1) == 0){
    message("No anchor regions is found")
    return(GRangesList())
  }
  # If the anchor GRanges is only polyA
  if(min_l1 > 6000){
    message("All regions have polyA only")
    return(GRangesList())
  }

  # Determining which starter to use
  starter_gr <- x[which(is_min_l1)] # Getting the starting GRange
  non_insert <- seqnames(starter_gr) != "Hot_L1_polyA"
  non_insert <- lapply(non_insert, as.logical)

  if(length(which(is_min_l1)) > 1){
    is_min_l1 <- seq_along(x)%in%which(is_min_l1)[which.max(max(width(insert_gr[is_min_l1])))]
    # message("Too many anchors for searching")
    # return(GRangesList())
  }
  right_match <- sapply(non_insert, tail, n = 1)
  left_match <- sapply(non_insert, head, n = 1)
  has_match <- mapply(function(x, y)any(c(x, y)), x = left_match, y = right_match)

  if(sum(has_match) > 1){
    has_long <- elementMetadata(starter_gr)$is_long
    if(any(has_long)){
      starter_gr <- starter_gr[has_long]
      right_match <- right_match[has_long]
      left_match <- left_match[has_long]
    }
  }

  if(length(starter_gr) > 1){
    starter_gr <- starter_gr[1]
    left_match <- left_match[1]
    right_match <- right_match[1]
  }

  # Determining the direction of the anchor cluster
  if(any(c(left_match, right_match))){
    start_direction <- ifelse(left_match, "left", "right")
  }else{
    return(GRangesList())
  }

  # Selecting genomic ranges for searching next GRanges
  starter <- starter_gr[seqnames(starter_gr) != "Hot_L1_polyA"][[1]]
  next_gr <- subsetByOverlaps(x[!is_min_l1], starter + 1000)#, ignore.strand=TRUE)

  # Checking whether they are in the right orientation
  .determine_direction <- function(q, z){
    tot <- seq_along(q)
    ol <- findOverlaps(q, z + 1000, ignore.strand = TRUE)@from
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

  next_direction <- sapply(next_gr, .determine_direction, z=starter)
  is_sensible <- next_direction!=start_direction & next_direction != "undetermined"

  if(sum(is_sensible)==0){
    message("no sensible end")
    return(x[is_min_l1][1])
  }

  sensible <- next_gr[is_sensible]
  included_seq <- elementMetadata(sensible)$anno[!grepl(":", elementMetadata(sensible)$anno)]
  included_seq <- included_seq[nchar(included_seq) > 5]
  included_seq <- DNAStringSetList(included_seq)
  has_polyA_seq <- unlist(lapply(included_seq, function(x){
    j <- Biostrings::letterFrequency(x, letters = c("A", "T", "G", 'C'))
    tot <- rowSums(j)
    a_prop <- j[,1] / tot
    t_prop <- j[,2] / tot
    any(a_prop>0.8|t_prop > 0.8)
    }))
  has_polyA <- has_polyA_seq|any(end(sensible[seqnames(sensible)=="Hot_L1_polyA"]) > 6000)
  sensible <- sensible[has_polyA]

  if(!any(has_polyA)){
    message("no sensible end")
    return(x[is_min_l1][1])
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
    end_partner <- sensible[1]
  }

  # Finding 3' transduction
  # initiation
  is_thing4search <- !(elementMetadata(x)$id%in%elementMetadata(next_gr)$id|is_min_l1)

  if(end_direction == "left"){
    next_target <- end_partner[[1]][-1]
    next_target <- next_target[seqnames(next_target) != "Hot_L1_polyA"]
  }else{
    next_target <- head(end_partner[[1]], -1)
    next_target <- next_target[seqnames(next_target) != "Hot_L1_polyA"]
  }

  middle_regions <- GRangesList()
  last_time <- rep(FALSE, length(x))
  while(sum(is_thing4search) > 0){
    thing4search <- x[is_thing4search]
    possible_middle <- subsetByOverlaps(thing4search, next_target+1000, ignore.strand = TRUE)
    if(length(possible_middle) == 0){
      break
    }

    ol_region <- GRangesList(lapply(possible_middle, function(x)subsetByOverlaps(x, next_target + 1000, ignore.strand = TRUE)))
    strands <- lapply(strand(ol_region), function(x)unique(as.character(x)))
    ol_region <- ol_region[lengths(strands) == 1]
    possible_middle <- possible_middle[lengths(strands) == 1]

    # Flipping the strand if the strand does not match but overlap
    .rc <- function(x){
      if(length(x)==0){
        return(GRangesList(GRanges()))
      }
      # x needs to be a GRangeList
      x <- revElements(x)
      original <- CharacterList(strand(x))
      strand(x)[original=="+"] <- "-"
      strand(x)[original=="-"] <- "+"
      elementMetadata(x)$anno <- revElements(elementMetadata(x)$anno)
      elementMetadata(x)$anno[grepl(":", elementMetadata(x)$anno)] <- CharacterList(lapply(x, as.character))
      y <- elementMetadata(x)$anno[!grepl(":", elementMetadata(x)$anno)]
      complemented <- CharacterList(lapply(y, function(x)as.character(complement(DNAStringSet(x)))))
      elementMetadata(x)$anno[!grepl(":", elementMetadata(x)$anno)] <- complemented
      elementMetadata(x)$rc <- TRUE
      x
    }

    strand_not_match <- as.character(strand(ol_region)) != as.character(strand(next_target))
    possible_middle[strand_not_match] <- .rc(possible_middle[strand_not_match])

    middle_direction <- sapply(possible_middle, .determine_direction, z = next_target)
    sensible_middle <- possible_middle[middle_direction!=start_direction&middle_direction != "undetermined"]

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

    next_target <- subsetByOverlaps(middle[[1]], next_target+1000, invert = TRUE)
    next_target <- next_target[seqnames(next_target) != "Hot_L1_polyA"]
    if(length(next_target)==0){
      break
    }
    is_thing4search <- !(elementMetadata(x)$id %in% elementMetadata(sensible_middle)$id | !is_thing4search)
    if(all(last_time == is_thing4search)|all(is_thing4search == FALSE)){
      break
    }
    last_time <- is_thing4search
  }

  if(start_direction == "left"){
    result <- c(starter_gr, middle_regions, end_partner)
  }

  if(start_direction == "right"){
    result <- c(end_partner, middle_regions, starter_gr)
  }
  return(result)
}
