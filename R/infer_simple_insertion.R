#' @export
infer_simple_insertion <- function(a, chromosome, target, tol = 10000){
  grl <- .preprocess_info(a, chromosome = c(chromosome, target), tol = tol)
  grl <- grl[any(S4Vectors::match(seqnames(grl), target, nomatch = 0) > 0)]
  cluster_regions <- convert_character2gr(elementMetadata(grl)$cluster_region)
  grl <- grl[!S4Vectors::match(seqnames(cluster_regions), target, nomatch = 0) > 0]
  grl <- .get_the_best(grl, elementMetadata(grl)$cluster_region)
  # names(grl) <- elementMetadata(grl)$cluster_region

  # Finding annotation that fully cover the insertion and they will be handled separately
  first_gr <- unlist(first_element(grl))
  last_gr <- unlist(last_element(grl))
  first_last_gr <- c(first_gr, last_gr)
  first_last_gr <- split(first_last_gr, as.integer(first_last_gr$origin))
  head_tail_ol <- how_many_regions_in_range(first_last_gr, tol = tol) == 1
  fully_covered <- grl[head_tail_ol]
  storage_1 <- .simple_initialize_storage(length(fully_covered))
  if(length(fully_covered) != 0){
    storage_1 <- simple_data_input(storage_1, fully_covered, direction = "left", targets = target)
    storage_1 <- simple_data_input(storage_1, fully_covered, direction = "right", targets = target)
  }

  # paired insertion
  p_clusters <- first_gr[!head_tail_ol & first_gr$origin_direction == "+"]
  m_clusters <- last_gr[!head_tail_ol & last_gr$origin_direction == "-"]
  if(length(p_clusters) > 0 & length(m_clusters) > 0){
    paired <- pairing_annotation(p_clusters, m_clusters)
  }else{
    paired <- Pairs(first = NULL, second = NULL)
  }

  storage_2 <- .simple_initialize_storage(length(paired))
  grl_first <- grl[as.character(paired@first$origin)]
  grl_second <- grl[as.character(paired@second$origin)]

  if(length(paired) != 0){
    storage_2 <- simple_data_input(storage_2, grl_first, direction = "left", targets = target)
    storage_2 <- simple_data_input(storage_2, grl_second, direction = "right", targets = target)
    storage_2$n_reads_5p <- elementMetadata(grl_first)$nreads
    storage_2$n_reads_3p <- elementMetadata(grl_second)$nreads
  }

  # unpaired positive cluster
  used_origins <- as.character(c(unlist(fully_covered)$origin, paired@first$origin, paired@second$origin))
  unpaired <- grl[!(names(grl) %in% used_origins)]
  unpaired_is_p <- unlist(first_element(unpaired))$origin_direction == "+"
  unpaired_p <- unpaired[unpaired_is_p]
  storage_3 <- .simple_initialize_storage(sum(unpaired_is_p))
  if(length(unpaired_p) != 0){
    storage_3 <- simple_data_input(storage_3, unpaired_p, direction = "left", targets = target)
  }
  # unpaired negative cluster
  unpaired_m <- unpaired[!unpaired_is_p]
  storage_4 <- .simple_initialize_storage(sum(!unpaired_is_p))
  if(length(unpaired_m) != 0){
    storage_4 <- simple_data_input(storage_4, unpaired_m, direction = "right", targets = target)
  }

  collected_origins <- split(unlist(first_element(fully_covered))$origin, seq_along(fully_covered))
  if(length(paired) > 0){
    collected_origins <- c(collected_origins, split(c(paired@first$origin, paired@second$origin), rep(seq_along(paired), 2)))
  }
  collected_origins <- c(collected_origins, split(unlist(first_element(unpaired_p))$origin, seq_along(unpaired_p)))
  collected_origins <- c(collected_origins, split(unlist(first_element(unpaired_m))$origin, seq_along(unpaired_m)))
  collected_origins <- IntegerList(collected_origins)
  storage <- rbind(storage_1, storage_2, storage_3, storage_4)
  return(list(collected_origins, storage))
}


simple_data_input <- function(storage, data, direction, targets){
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
.simple_initialize_storage <- function(n){
  storage <- list("genomic_qname_5p" = NA,
                  "genomic_rname_5p" = NA,
                  "genomic_start_5p" = NA,
                  "genomic_end_5p" = NA,
                  "genomic_orientation_5p" = NA,
                  "insert_qname_5p" = NA,
                  "insert_rname_5p" = NA,
                  "insert_start_5p" = NA,
                  "insert_end_5p" = NA,
                  "insert_orientation_5p" = NA,
                  "is_exact_5p" = NA,
                  "n_reads_5p" = NA,
                  "genomic_qname_3p" = NA,
                  "genomic_rname_3p" = NA,
                  "genomic_start_3p" = NA,
                  "genomic_end_3p" = NA,
                  "genomic_orientation_3p" = NA,
                  "insert_qname_3p" = NA,
                  "insert_rname_3p" = NA,
                  "insert_start_3p" = NA,
                  "insert_end_3p" = NA,
                  "insert_orientation_3p" = NA,
                  "is_exact_3p" = NA,
                  "n_reads_3p" = NA)
  setDT(storage)
  if(n == 0){
    return(storage[-1,])
  }
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
