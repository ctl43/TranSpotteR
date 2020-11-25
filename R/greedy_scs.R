#' @export
#' @importFrom data.table data.table
#' @importFrom Biostrings consensusMatrix

# Finding the overlapping percentage and overlapping length between pair

overlap_reads <- function(vec, min_pid = 85, min_len = 10){
  if(is.null(names(vec))){
    names(vec) <- seq_along(vec)
  }
  if(length(vec)<=1){
    # return empty data.table
    return(data.table(
      "seq1"=character(),
      "seq2"=character(),
      "seq1_left_clipped_len" = integer(),
      "seq1_right_clipped_len" = integer(),
      "seq2_left_clipped_len" = integer(),
      "seq2_right_clipped_len" = integer(),
      "seq1_len" = integer(),
      "seq2_len" = integer(),
      "seq1_ol_start" = integer(),
      "seq1_ol_end" = integer(),
      "seq2_ol_start" = integer(),
      "seq2_ol_end" = integer(),
      "pid" = numeric(),
      "n_match" = integer(),
      "aln_ol_len" = integer()
      ))
  }
  n_seq <- length(vec)
  idx <- names(vec)
  info <- data.table(t(combn(idx, 2)))
  colnames(info) <- c("seq1", "seq2")
  info <- cbind(info, do.call(cbind, overlapper(vec[info[[1]]], vec[info[[2]]])))
  # Filtering abberent alignment
  info <- info[info$pid > min_pid & info$aln_ol_len > min_len & (info$seq1_ol_start==0 | info$seq2_ol_start==0), ]
  info <- info[order(info$aln_ol_len, decreasing=TRUE), ]

  # Converting to 1-based index
  info$seq1_ol_start <- info$seq1_ol_start + 1
  info$seq1_ol_end <- info$seq1_ol_end + 1
  info$seq2_ol_start <- info$seq2_ol_start + 1
  info$seq2_ol_end <- info$seq2_ol_end + 1
  info <- .reorder_seq_info(info)
  info
}

#' @export
greedy_scs <- function(vec, msa_result = FALSE, return_no_assembly = FALSE, add_id = TRUE){
  # Dealing with NULL vector
  if(is.null(vec)){
    consensus <- NULL
    msa_out <- list()
    if(msa_result){
      return(list(consensus = consensus, msa = msa_out))
    }else{
      return(consensus)
    }
  }

  names(vec) <- seq_along(vec)
  original_vec <- vec
  # This can be vastly improved by not calculating the overlaps that have been calculated before
  # Will be implemented into C++ in the future

  # Creating a drafted assembled sequence for MSA by greedy shortest supersequence (SCS)
  for(i in seq_along(vec)){ # Maximum number of loops
    ol <- overlap_reads(vec)
    if(nrow(ol)>0){
      selected <- ol[1, ]
      new_seq <- paste0(substr(vec[[selected$seq1]], 1, selected$seq1_ol_start - 1L),
                        substr(vec[[selected$seq2]], selected$seq2_ol_start, selected$seq2_len))
      names(new_seq) <- paste0(selected$seq1, "_", selected$seq2)
      vec <- c(new_seq, vec[!names(vec)%in%c(selected$seq1, selected$seq2)])
    }else{
      break
    }
  }
  # Picking out those read that did not assemble
  ass_id <- names(vec)
  solo_id <- ass_id[!grepl("_", ass_id)]
  non_solo_id <- grep("_", ass_id, value = TRUE)
  member_idx <- strsplit(non_solo_id, "_")
  msa_out <- mapply(function(x,y)t_coffee(c(original_vec[x],y)),
                    x = member_idx, y = vec[non_solo_id], SIMPLIFY = FALSE)
  consensus <- sapply(msa_out, .process_msa)
  names(consensus) <- non_solo_id
  names(msa_out) <- non_solo_id

  consensus <- unlist(consensus)

  if(add_id & length(consensus) > 0){
    names(consensus) <- seq_along(consensus)
  }

  if(return_no_assembly){
    consensus <- c(consensus, vec[solo_id])
    msa_out <- c(msa_out, split(vec[solo_id], solo_id))
  }

  if(msa_result){
    return(list(consensus = consensus, msa = msa_out))
  }else{
    return(consensus)
  }
}

# Function to determine a pair of sequences that which one is of the left hand side
.reorder_seq_info <- function(info){
  # Determining whether seq1 is on the left side
  seq1_on_right <- info$seq1_right_clipped_len!=0
  ordered_info <- info # Creating a copy
  info_name_1 <- c("seq1", "seq1_left_clipped_len", "seq1_right_clipped_len", "seq1_len", "seq1_ol_start", "seq1_ol_end")
  info_name_2 <- sub("seq1","seq2", info_name_1)

  for(i in seq_along(info_name_1)){
    ordered_info[[info_name_1[i]]][seq1_on_right] <- info[[info_name_2[i]]][seq1_on_right]
  }

  for(i in seq_along(info_name_2)){
    ordered_info[[info_name_2[i]]][seq1_on_right] <- info[[info_name_1[i]]][seq1_on_right]
  }
  ordered_info
}

.process_msa <- function(x){
  x <- replcae_space(x) # replacing the space at the beginning and at the end
  x <- consensusMatrix(x)
  x <- x[rownames(x)!=" ",,drop = FALSE]
  x <- paste(rownames(x)[apply(x, 2, which.max)], collapse = "")
  gsub("-","",x)
}
