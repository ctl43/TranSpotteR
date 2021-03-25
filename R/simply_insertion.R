#' @export
simply_insertion <- function(a, targets, tol = 10000){
  grl <- .preprocess_info(a)
  grl <- grl[any(seqnames(grl) %in% targets)]
  
  # Finding annotation that fully cover the insertion and they will be handled separately
  first_gr <- unlist(first_element(grl))
  last_gr <- unlist(last_element(grl))
  first_last_gr <- c(first_gr, last_gr)
  first_last_gr <- split(first_last_gr, as.integer(first_last_gr$origin))
  head_tail_ol <- how_many_regions_in_range(first_last_gr, tol = tol) == 1
  fully_covered <- grl[head_tail_ol]
  fully_covered_ranges <- split(fully_covered_ranges, strand(fully_covered_ranges))
  fully_covered_ranges$`-` <- subsetByOverlaps(fully_covered_ranges$`-`, fully_covered_ranges$`+` + tol, ignore.strand = TRUE, invert = TRUE)
  fully_covered_ranges <- unlist(fully_covered_ranges)
  fully_covered <- fully_covered[fully_covered_ranges$origin]
  storage_1 <- .simply_initialize_storage(length(fully_covered))
  storage_1 <- simply_data_input(storage_1, fully_covered, direction = "left")
  storage_1 <- simply_data_input(storage_1, fully_covered, direction = "right")
  
  # paired insertion
  p_clusters <- first_gr[!head_tail_ol & first_gr$origin_direction == "+"]
  m_clusters <- last_gr[!head_tail_ol & last_gr$origin_direction == "-"]
  paired <- pairing_annotation(p_clusters, m_clusters)
  storage_2 <- .initialize_storage(length(paired))
  storage_2 <- data_input(storage_2, paired@first, direction = "left")
  storage_2 <- data_input(storage_2, paired@second, direction = "right")
  storage_2$n_reads_5p <- elementMetadata(paired@first)$nreads
  storage_2$n_reads_3p <- elementMetadata(paired@second)$nreads
  
  # unpaired positive cluster
  used_origins <- c(used_origins, as.character(c(fully_covered$origin, paired@first$origin, paired@second$origin)))
  unpaired <- grl[!(names(grl) %in% used_origins)]
  unpaired_is_p <- unlist(first_element(unpaired))$origin_direction == "+"
  unpaired_p <- unpaired[unpaired_is_p]
  storage_3 <- .simply_initialize_storage(sum(unpaired_is_p))
  storage_3 <- simply_data_input(storage_3, unpaired_p, direction = "left")
  
  # unpaired negative cluster
  unpaired_m <- unpaired[!unpaired_is_p]
  storage_4 <- .initialize_storage(sum(!unpaired_is_p))
  storage_4 <- data_input(storage_4, unpaired[!unpaired_is_p], direction = "right")
  
  collected_origins <- split(unlist(first_element(fully_covered))$origin, seq_along(fully_covered))
  collected_origins <- c(collected_origins, split(c(paired@first$origin, paired@second$origin), rep(seq_along(paired), 2)))
  collected_origins <- c(collected_origins, split(unlist(first_element(unpaired_p))$origin, seq_along(unpaired_p)))
  collected_origins <- c(collected_origins, split(unlist(first_element(unpaired_m))$origin, seq_along(unpaired_m)))
  collected_origins <- IntegerList(collected_origins)
  
  storage <- rbind(storage_1, storage_2, storage_3, storage_4)
  return(list(collected_origins, storage))
}


simply_data_input <- function(storage, data, direction){
  if(nrow(storage) != length(data)){
    stop()
  }
  data_slot <- c("qname", "rname", "start", "end", "orientation")
  data_slot <- paste0(data_slot, ifelse(direction == "left", "_5p", "_3p"))
  
  for(type in c("insert", "genomic")){
    if(direction == "right"){
      if(type == "insert"){
        selected <- first_element(non_targets(data, targets = targets, invert = TRUE))
      }
      if(type == "genomic"){
        selected <- last_element(non_targets(data, targets = targets))
      }
    }else{
      if(type == "insert"){
        selected <- last_element(non_targets(data, targets = targets, invert = TRUE))
      }
      if(type == "genomic"){
        selected <- first_element(non_targets(data, targets = targets))
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
.simply_initialize_storage <- function(n){
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
                  "insert_break_5p" = NA,
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
                  "is_exact_3p" = NA,
                  "n_reads_3p" = NA)
  setDT(storage)
  storage <- rep(list(storage), n)
  rbindlist(storage)
}

#' @export
#' @importFrom BiocGenerics match
non_targets <- function(x, targets ,invert = FALSE, return_TF = FALSE){
  if(invert){
    test <- match(seqnames(x), targets, nomatch = 0) > 0
  }else{
    test <- !(match(seqnames(x), targets, nomatch = 0) > 0)
  }
  if(return_TF){
    return(test)
  }else{
    return(x[test])
  }
}