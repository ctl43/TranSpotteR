#' @export
infer_transposon <- function(a, tol = 10000, max_transduced = 50000, chromosome = chromosome){
  grl <- suppressWarnings(.preprocess_info(a, chromosome = chromosome, tol = tol))

  ## Filtering out annotation with polyA only
  insert_gr <- non_polyA_or_insert(grl)
  anno_is_polyA_only <- all(.check_polyA(insert_gr))
  grl <- grl[!anno_is_polyA_only]

  # Getting the best annotation if there is more than one annotation in a read cluster
  grl <- .get_the_best(grl, elementMetadata(grl)$cluster_region)

  # Finding annotation that fully cover the insertion and they will be handled separately
  first_gr <- unlist(first_element(grl))
  last_gr <- unlist(last_element(grl))
  first_last_gr <- c(first_gr, last_gr)
  first_last_gr <- split(first_last_gr, as.integer(first_last_gr$origin))
  head_tail_ol <- how_many_regions_in_range(first_last_gr, tol = tol) == 1
  fully_covered <- grl[head_tail_ol]
  used_origins <- as.character(unlist(fully_covered)$origin)

  # Checking whether they have polyA
  has_pa_1 <- .check_polyA(fully_covered)
  fully_covered <- fully_covered[any(has_pa_1)]
  has_pa_1 <- has_pa_1[any(has_pa_1)]

  # Checking whether there is only polyA
  has_multiple_1 <- how_many_regions_in_range(fully_covered[!has_pa_1], tol = tol) != 1
  fully_covered <- fully_covered[has_multiple_1]
  if(length(fully_covered) != 0){
    storage_1 <- .initialize_storage(length(fully_covered))
    storage_1 <- data_input(storage_1, fully_covered, direction = "left")
    storage_1 <- data_input(storage_1, fully_covered, direction = "right")
    storage_1 <- .get_TSD(unlist(first_element(fully_covered)), unlist(last_element(fully_covered)), storage_1)
  }
  # Need to handle transduction!
  head_tail_trimmed <- (last_element(first_element(fully_covered, invert = TRUE), invert = TRUE))
  head_is_pa <- unlist(first_element(.check_polyA(head_tail_trimmed)), use.names = FALSE)
  tail_is_pa <- unlist(last_element(.check_polyA(head_tail_trimmed)), use.names = FALSE)
  head_is_pa[head_is_pa & tail_is_pa] <- FALSE
  tail_is_pa[head_is_pa & tail_is_pa] <- FALSE
  fc_left_transduced <- merge_glist(non_polyA_or_insert(head_tail_trimmed[head_is_pa]), tol = max_transduced)
  fc_right_transduced <- merge_glist(non_polyA_or_insert(head_tail_trimmed[tail_is_pa]), tol = max_transduced)
  fc_left_transduced <- grl_to_character(fc_left_transduced)
  fc_right_transduced <- grl_to_character(fc_right_transduced)
  storage_1$transduced_genomic_region_5p[as.logical(head_is_pa)] <- fc_left_transduced
  storage_1$transduced_genomic_region_3p[as.logical(tail_is_pa)] <- fc_right_transduced

  # For reads that are not fully covered by a contig
  p_clusters <- first_gr[!head_tail_ol & first_gr$origin_direction == "+"]
  m_clusters <- last_gr[!head_tail_ol & last_gr$origin_direction == "-"]
  paired <- pairing_annotation(p_clusters, m_clusters)

  ## Filtering pairs that do not have LINE1 and polyA sequence
  paired_p <- grl[as.character(paired@first$origin)]
  paired_m <- grl[as.character(paired@second$origin)]
  grp_p <- rep(seq_along(paired_p), lengths(paired_p))
  grp_m <- rep(seq_along(paired_m), lengths(paired_m))
  combined_paired <- split(unlist(c(paired_p, paired_m)), c(grp_p, grp_m))
  is_selected <- any(match(seqnames(combined_paired), c("Hot_L1_polyA", "polyA"), nomatch = 0)>0)
  has_pa_p <- any(.check_polyA(paired_p))
  has_pa_m <- any(.check_polyA(paired_m))
  has_insert_p <- any(seqnames(paired_p)=="Hot_L1_polyA")
  has_insert_m <- any(seqnames(paired_m)=="Hot_L1_polyA")

  # Paired insertions
  has_pa_2 <- (has_pa_p | has_pa_m & !(has_pa_p & has_pa_m)) #does not allow to have polyA at both ends
  is_selected <- has_insert_p | has_insert_m & has_pa_2
  selected_pairs <- paired[is_selected]
  if(length(selected_pairs) != 0 ){
    selected_p <- paired_p[is_selected]
    selected_m <- paired_m[is_selected]
    end_on_the_right <- any(.check_polyA(selected_m))
    storage_2 <- .initialize_storage(length(selected_p))
    storage_2 <- data_input(storage_2, selected_p, direction = "left")
    storage_2 <- data_input(storage_2, selected_m, direction = "right")
    storage_2$n_reads_5p <- elementMetadata(selected_p)$nreads
    storage_2$n_reads_3p <- elementMetadata(selected_m)$nreads
    # Handling transduction on the right
    right_transduced <- non_polyA_or_insert(last_element(selected_m[end_on_the_right], invert = TRUE))
    right_transduced <- right_transduced[width(range(right_transduced)) < max_transduced]
    right_transduced <- grl_to_character(right_transduced)
    right_transduced[right_transduced == ""] <- NA
    storage_2$transduced_genomic_region_3p[end_on_the_right] <-  right_transduced

    # Handling transduction on the left
    left_transduced <- non_polyA_or_insert(first_element(selected_p[!end_on_the_right], invert = TRUE))
    left_transduced <- left_transduced[width(range(left_transduced)) < max_transduced]
    left_transduced <- grl_to_character(left_transduced)
    left_transduced[left_transduced == ""] <- NA
    storage_2$transduced_genomic_region_5p[!end_on_the_right] <-  left_transduced
    storage_2 <- .get_TSD(selected_pairs@first, selected_pairs@second, storage_2)
  }


  # Processing those without end
  used_origins <- c(used_origins, as.character(c(selected_pairs@first$origin, selected_pairs@second$origin)))
  unpaired <- grl[!(names(grl) %in% used_origins)]
  unpaired_has_pa <- .check_polyA(unpaired)
  unpaired_has_multiple <- how_many_regions_in_range(unpaired[!unpaired_has_pa], tol = tol) > 1
  unpaired_has_multiple <- any(unpaired_has_pa) & unpaired_has_multiple
  unpaired_has_insert <- any(seqnames(unpaired) == "Hot_L1_polyA" & start(unpaired) < 6000)
  is_selected_unpaired <- unpaired_has_insert|unpaired_has_multiple
  unpaired <- unpaired[is_selected_unpaired]

  ## For the unpaired positive clusters
  unpaired_is_p <- unlist(first_element(unpaired))$origin_direction == "+"
  unpaired_p <- unpaired[unpaired_is_p]
  if(length(unpaired_p) != 0){
    storage_3 <- .initialize_storage(sum(unpaired_is_p))
    storage_3 <- data_input(storage_3, unpaired_p, direction = "left")
    ### Transduction information on the 5' side
    unpaired_p_has_pa <- any(.check_polyA(unpaired_p))
    unpaired_p_has_multiple <- how_many_regions_in_range(non_polyA_or_insert(unpaired_p), tol = tol) > 1
    is_unpaired_p_transduced <- unpaired_p_has_pa & unpaired_p_has_multiple
    unpaired_p_transduced <- first_element(non_polyA_or_insert(unpaired_p[is_unpaired_p_transduced]), invert = TRUE)
    unpaired_p_transduced <- unpaired_p_transduced[width(range(unpaired_p_transduced)) < max_transduced]
    unpaired_p_transduced_grp <- rep(seq_along(unpaired_p_transduced), lengths(unpaired_p_transduced)) # in case there are more than 1 transduced regions
    unpaired_p_transduced <- CharacterList(split(as.character(unlist(unpaired_p_transduced, use.names = FALSE)), unpaired_p_transduced_grp))
    unpaired_p_transduced <- paste(unpaired_p_transduced, collapse = ",")
    storage_3$transduced_genomic_region_5p[is_unpaired_p_transduced] <- unpaired_p_transduced
  }



  ## For the unpaired negative clusters
  unpaired_m <- unpaired[!unpaired_is_p]
  if(length(unpaired_m) != 0){
    storage_4 <- .initialize_storage(sum(!unpaired_is_p))
    storage_4 <- data_input(storage_4, unpaired[!unpaired_is_p], direction = "right")

    ### Transduction information on the 3' side
    unpaired_m_has_pa <- any(.check_polyA(unpaired_m))
    unpaired_m_has_multiple <- how_many_regions_in_range(non_polyA_or_insert(unpaired_m), tol = tol) > 1
    is_unpaired_m_transduced <- unpaired_m_has_pa & unpaired_m_has_multiple
    unpaired_m_transduced <- range(last_element(non_polyA_or_insert(unpaired_m[is_unpaired_m_transduced]), invert = TRUE))
    unpaired_m_transduced_grp <- rep(seq_along(unpaired_m_transduced), lengths(unpaired_m_transduced))
    unpaired_m_transduced <- CharacterList(split(as.character(unlist(unpaired_m_transduced, use.names = FALSE)), unpaired_m_transduced_grp))
    unpaired_m_transduced <- paste(unpaired_m_transduced, collapse = ",")
    storage_4$transduced_genomic_region_3p[is_unpaired_m_transduced] <- unpaired_m_transduced
  }

  # Collecting all the origins used to predict the insertions
  collected_origins <- split(unlist(first_element(fully_covered))$origin, seq_along(fully_covered))
  collected_origins <- c(collected_origins, split(c(selected_pairs@first$origin, selected_pairs@second$origin), rep(seq_along(selected_pairs), 2)))
  collected_origins <- c(collected_origins, split(unlist(first_element(unpaired_p))$origin, seq_along(unpaired_p)))
  collected_origins <- c(collected_origins, split(unlist(first_element(unpaired_m))$origin, seq_along(unpaired_m)))
  collected_origins <- IntegerList(collected_origins)
  storage <- rbind(storage_1, storage_2, storage_3, storage_4)
  return(list(collected_origins, storage))
}

#' @export
#' @importFrom BiocGenerics match
.check_polyA <- function(x){
  tmp_1 <- match(seqnames(x), "Hot_L1_polyA", nomatch = 0) > 0
  tmp_2 <- end(x) > 6000
  tmp_3 <- match(seqnames(x), "polyA", nomatch = 0) > 0
  (tmp_1 & tmp_2)|tmp_3
}

#' @export
#' @importFrom BiocGenerics match
non_polyA_or_insert <- function(x, invert = FALSE, return_TF = FALSE){
  if(invert){
    test <- match(seqnames(x), c("polyA", "Hot_L1_polyA"), nomatch = 0) > 0
  }else{
    test <- !(match(seqnames(x), c("polyA", "Hot_L1_polyA"), nomatch = 0) > 0)
  }
  if(return_TF){
    return(test)
  }else{
    return(x[test])
  }
}

#' @export
.get_the_best <- function(grl, grl_grps, tol = 10000){
  # More informative > contig length
  ## Getting the number of information within a range
  n_regions <- unname(how_many_regions_in_range(grl, tol = tol))
  anno_len <- unname(sum(width(grl)))
  region_grp <- factor(grl_grps, levels = unique(grl_grps))
  max_n <- max(IntegerList(split(n_regions, region_grp)))
  max_n <- rep(max_n, table(region_grp))
  is_max_score <- (n_regions == max_n) * 10

  ## Combining all width of annotations
  max_len <- max(IntegerList(split(anno_len, region_grp)))
  max_len <- rep(max_len, table(region_grp))
  is_max_len_score <- (anno_len == max_len) * 1

  ## Computing the score
  scores <- unname(is_max_score + is_max_len_score)
  max_score <- max(IntegerList(split(scores, region_grp)))
  max_score <- unname(rep(max_score, table(region_grp)))
  get_which <- LogicalList(split(scores == max_score, region_grp))
  grl[unlist(!duplicated(get_which) & get_which)]
}

#' @export
pairing_annotation <- function(p_clusters, m_clusters, tol = 10000){
  combined_gr <- c(p_clusters, m_clusters)
  combined_gr <- combined_gr[order(unstrand(combined_gr))]
  combined_gr$order <- seq_along(combined_gr)
  splitted_gr <- split(combined_gr, combined_gr$origin_direction)
  p_start <- resize(unstrand(splitted_gr[[1]]), width = 1)
  m_end <- resize(unstrand(splitted_gr[[2]]), fix = "end", width = 1)
  preceded <- precede(p_start, m_end)
  anchor <- splitted_gr[[1]][!is.na(preceded)]
  preceded[is.na(preceded)] <- 0
  preceded <- splitted_gr[[2]][preceded]
  paired <- Pairs(anchor, preceded)
  genomic_dist <- distance(paired@first, paired@second)
  ordinal_dist <- preceded$order - anchor$order
  paired[genomic_dist < tol & ordinal_dist == 1]
}

#' @export
data_input <- function(storage, data, direction){
  if(nrow(storage) != length(data)){
    stop()
  }
  data_slot <- c("qname", "rname", "start", "end", "orientation")
  data_slot <- paste0(data_slot, ifelse(direction == "left", "_5p", "_3p"))

  for(type in c("insert", "genomic")){
    if(direction == "right"){
      if(type == "insert"){
        selected <- first_element(non_polyA_or_insert(data, invert = TRUE))
      }
      if(type == "genomic"){
        selected <- last_element(non_polyA_or_insert(data))
      }
    }else{
      if(type == "insert"){
        selected <- last_element(non_polyA_or_insert(data, invert = TRUE))
      }
      if(type == "genomic"){
        selected <- first_element(non_polyA_or_insert(data))
      }
    }
    selected_slot <- paste0(type, "_", data_slot)
    has_sth <- lengths(selected) != 0
    selected <- unlist(selected)
    storage[[selected_slot[1]]][has_sth] <- selected$QNAME
    storage[[selected_slot[2]]][has_sth] <- as.character(seqnames(selected))
    storage[[selected_slot[3]]][has_sth] <- start(selected)
    storage[[selected_slot[4]]][has_sth] <- end(selected)
    storage[[selected_slot[5]]][has_sth] <- as.character(unlist(strand(selected)))
  }

  storage[[paste0("n_reads", ifelse(direction == "left", "_5p", "_3p"))]] <- elementMetadata(data)$nreads
  return(storage)
}

#' @export
.initialize_storage <- function(n){
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
  setDT(storage)
  storage <- rep(list(storage), n)
  rbindlist(storage)
}

#' @export
#' @importFrom stringr str_sub
#' @importFrom IRanges LogicalList
#' @importFrom GenomicAlignments GAlignments mapToAlignments pmapToAlignments
.get_TSD <- function(left_gr, right_gr, storage){
  ol <- findOverlaps(left_gr, right_gr)
  ol <- ol[ol@from == ol@to]
  selected <- ol@from
  selected_left <- left_gr[selected]
  selected_right <- right_gr[selected]
  lr_intersect <- pintersect(selected_left, selected_right)
  if(length(lr_intersect) != 0){
    lr_gr <- c(selected_left, selected_right)
    # Translating the overlapping regions to the read position
    l_alignments <- GAlignments(seqnames(selected_left),
                              pos = start(selected_left),
                              cigar = c(selected_left$cigar),
                              strand = strand(selected_left),
                              names = c(paste0("left_", selected)))
    r_alignments <- GAlignments(seqnames(selected_right),
                                pos = start(selected_right),
                                cigar = c(selected_right$cigar),
                                strand = strand(selected_right),
                                names = c(paste0("right_", selected)))
    l_mapped_to_reads <- pmapToAlignments(lr_intersect, l_alignments)
    r_mapped_to_reads <- pmapToAlignments(lr_intersect, r_alignments)
    left_ok <- end(l_mapped_to_reads) == cigarWidthAlongQuerySpace(selected_left$cigar)
    right_ok <- start(r_mapped_to_reads) == 1
    both_ok <- left_ok & right_ok
    if(any(both_ok)){
      selected <- selected[both_ok]
      selected_left <- selected_left[both_ok]
      selected_right <- selected_right[both_ok]
      l_mapped_to_reads <- l_mapped_to_reads[both_ok]
      r_mapped_to_reads <- r_mapped_to_reads[both_ok]
      lr_intersect <- lr_intersect[both_ok]
      l_genomic_overlap_seq <- str_sub(c(selected_left$seq),
                                     start = start(l_mapped_to_reads), end = end(l_mapped_to_reads))
      r_genomic_overlap_seq <- str_sub(c(selected_right$seq),
                                       start = start(r_mapped_to_reads), end = end(r_mapped_to_reads))
      storage[["overlapping_genomic_region"]][selected] <- as.character(lr_intersect)
      storage[["duplicated_seq_5p"]][selected] <- l_genomic_overlap_seq
      storage[["duplicated_seq_3p"]][selected] <- r_genomic_overlap_seq
    }
  }
  return(storage)
}

#' @export
.preprocess_info <- function(a, chromosome = c(1:22, "X", "Y", "Hot_L1_polyA", "polyA"), tol){
  x <- a[[1]]
  n_anno <- sapply(x, nrow)
  origins <- factor(rep(seq_along(x), n_anno), levels = seq_along(x))
  cluster_gr <- convert_character2gr(a$cluster_region)
  origin_direction <- rep(strand(cluster_gr), n_anno)
  x <- rbindlist(x)
  x$origin <- origins
  x$origin_direction <- as.factor(origin_direction)
  selected <- grep(":", x$annotation)
  y <- x[selected]
  extra_info <- list(QNAME = y$QNAME, origin = y$origin,
                     origin_direction = y$origin_direction,
                     cigar = y$cigar, seq = y$seq)
  gr <- convert_character2gr(y$annotation, extra_info = extra_info)
  grl <- split(gr, gr$origin)
  elementMetadata(grl) <- a

  has_multiple <- how_many_regions_in_range(grl, tol = tol) > 1
  no_unwanted_rname <- all(match(split(seqnames(gr), gr$origin), chromosome, nomatch = 0) > 0)
  grl <- grl[has_multiple & no_unwanted_rname]

  # Getting annotations that fit the cluster orientation
  sensible_things <- first_things <- unlist(first_element(grl))
  last_things <- unlist(last_element(grl))
  is_m <- last_things$origin_direction == "-"
  sensible_things[is_m] <- last_things[is_m]
  gr4ol <- convert_character2gr(elementMetadata(grl)$cluster_region)
  strand(gr4ol) <- "+"
  ol <- findOverlaps(sensible_things, gr4ol)
  ol_idx <- ol@from[ol@from == ol@to]
  grl <- grl[ol_idx]
  grl
}

