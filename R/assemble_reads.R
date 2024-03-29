#' @export
#' @importFrom data.table data.table setDT setkeyv
#' @importFrom Biostrings consensusMatrix DNAStringSet
#' @importFrom uuid UUIDgenerate
assemble_reads <- function(vec, n_reads = NULL, msa_result = FALSE, add_id = TRUE, min_len = 8L,
                       min_pid = 85, consensus_min_len = 100, return_hidden_supported = FALSE)
  # It assembles reads (getting the shortest common superstring problem, scs) by overlap-layout-consensus method.
  # The overlapping part, it simply chooses the pair with longest overlapping length, so called a greedy way.
  # Written by Cheuk-Ting Law
{
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
  if(is.null(n_reads)){
    n_reads <- rep(1L, length(vec))
  }

  idx <- names(n_reads) <- names(vec) <- seq_along(vec)
  original_vec <- vec

  info <- data.table(t(combn(names(vec), 2)))
  colnames(info) <- c("seq1", "seq2")
  info <- data.table(info, setDT(overlapper(vec[info[[1]]], vec[info[[2]]])))
  info <- .reorder_seq_info(info)
  # Getting alignment with high similarity
  info <- info[info$pid > min_pid & info$aln_ol_len > min_len & (info$seq1_ol_start == 1L | info$seq2_ol_start == 1L), ]
  info <- info[order(info$aln_ol_len, decreasing = TRUE), ]

  # Filtering out reads that are fully covered by other reads
  fully_covered <- which_is_fully_covered(info)
  info <- info[!(seq1 %in% fully_covered[[1]] | seq2 %in% fully_covered[[1]])]
  vec <- vec[!names(vec) %in% fully_covered[[1]]]
  supported_by_fully_covered <- fully_covered[[2]]

  while(nrow(info) > 0){
    selected <- info[1L, ]
    left_idx <- selected$seq1
    right_idx <- selected$seq2
    new_seq <- paste0(substr(vec[left_idx], 1L, selected$seq1_ol_start - 1L),
                      substr(vec[right_idx], selected$seq2_ol_start, selected$seq2_len))
    names(new_seq) <- paste0(selected$seq1, "_", selected$seq2)
    vec <- vec[!names(vec) %in% c(selected$seq1, selected$seq2)]

    # Only overlapping the newly combined sequence and the total sequence
    ol <- overlapper(rep(new_seq, length(vec)), vec)
    suppressWarnings(ol <- data.table(data.table(seq1 = names(new_seq), seq2 = names(vec)), ol))
    ol <- .reorder_seq_info(ol)
    ol <- ol[ol$pid > min_pid & ol$aln_ol_len > min_len & (ol$seq1_ol_start == 1L | ol$seq2_ol_start == 1L), ]
    fully_covered <- which_is_fully_covered(ol)
    vec <- c(new_seq, vec)
    info <- info[!(seq1 %in% c(left_idx, right_idx, fully_covered[[1]]) | seq2 %in% c(left_idx, right_idx, fully_covered[[1]]))]
    supported_by_fully_covered <- c(supported_by_fully_covered, fully_covered[[2]])
    info <- rbindlist(list(ol, info))
    info <- info[order(info$aln_ol_len, decreasing = TRUE), ]
  }

  # Picking out those read that did not assemble
  ass_id <- names(vec)
  solo_id <- ass_id[!grepl("_", ass_id)]
  non_solo_id <- grep("_", ass_id, value = TRUE)
  member_idx <- strsplit(non_solo_id, "_")

  pairwise_aln <- mapply(function(x, y) overlapper(rep(y, length(x)), original_vec[x]), x = member_idx, y = vec[non_solo_id], SIMPLIFY = FALSE)
  msa_view_aln <- lapply(pairwise_aln, function(x)msa_view(x[["seq1_aln"]], x[["seq2_aln"]]))
  consensus <- unlist(sapply(msa_view_aln, .process_msa))
  names(msa_view_aln) <- non_solo_id

  if(length(consensus) > 0L){
    if(add_id ){
      names(msa_view_aln) <- names(consensus) <- UUIDgenerate(n = length(consensus))
    }else{
      names(consensus) <- non_solo_id
    }
  }


  # Computing the number of reads consisting the consensus sequence
  n_reads_out <- unlist(lapply(member_idx, function(x)sum(n_reads[x])))

  if(!is.null(consensus_min_len)){
    n_reads_out <- n_reads_out[nchar(consensus) >= consensus_min_len]
    msa_view_aln <- msa_view_aln[nchar(consensus) >= consensus_min_len]
    consensus <- consensus[nchar(consensus) >= consensus_min_len]
  }

  if(return_hidden_supported){
    rescued <- vec[names(vec) %in% supported_by_fully_covered]
    rescued_n_reads <- n_reads[names(rescued)]
    rescued_msa_view <- split(rescued, names(rescued))
    if(add_id){
      names(rescued_msa_view) <- names(rescued) <- UUIDgenerate(n = length(rescued))
    }
    msa_view_aln <- c(rescued_msa_view, msa_view_aln)
    consensus <- c(rescued, consensus)
    n_reads_out <- c(rescued_n_reads, n_reads_out)
  }

  if(msa_result){
    return(list(consensus = consensus, n_reads = n_reads_out, msa = msa_view_aln))
  }else{
    return(list(consensus = consensus, n_reads = n_reads_out))
  }
}

which_is_fully_covered <- function(info){
  fully_covered_1 <- info$seq1_left_clipped_len == 0L & info$seq1_right_clipped_len ==0L
  fully_covered_2 <- info$seq2_left_clipped_len == 0L & info$seq2_right_clipped_len ==0L
  fully_covered_1[fully_covered_1 & fully_covered_2] <- FALSE
  init <- rep(NA, sum(fully_covered_1, fully_covered_2))
  fully_covered <- data.frame(fully_covered = init, by = init)
  fully_covered_1 <- data.table(fully_covered = info$seq1[fully_covered_1], by = info$seq2[fully_covered_1])
  fully_covered_2 <- data.table(fully_covered = info$seq2[fully_covered_2], by = info$seq1[fully_covered_2])
  rbind(fully_covered_1, fully_covered_2)
}


.reorder_seq_info <- function(info){
  # Function to determine a pair of sequences that which one is of the left hand side
  # Determining whether seq1 is on the left side
  seq1_on_right <- info$seq1_right_clipped_len != 0L
  ordered_info <- info # Creating a copy
  info_name_1 <- c("seq1", "seq1_left_clipped_len", "seq1_right_clipped_len", "seq1_len", "seq1_ol_start", "seq1_ol_end")
  info_name_2 <- sub("seq1","seq2", info_name_1)

  for(i in seq_along(info_name_1)){
    ordered_info[[info_name_1[i]]][seq1_on_right] <- info[[info_name_2[i]]][seq1_on_right]
  }

  for(i in seq_along(info_name_2)){
    ordered_info[[info_name_2[i]]][seq1_on_right] <- info[[info_name_1[i]]][seq1_on_right]
  }

  # Converting to 1-based index
  ordered_info$seq1_ol_start <- ordered_info$seq1_ol_start + 1L
  ordered_info$seq1_ol_end <- ordered_info$seq1_ol_end + 1L
  ordered_info$seq2_ol_start <- ordered_info$seq2_ol_start + 1L
  ordered_info$seq2_ol_end <- ordered_info$seq2_ol_end + 1L
  setkeyv(ordered_info, c("seq1", "seq2"))
  ordered_info
}

.process_msa <- function(x){
  x <- replcae_space(x) # replacing the space at the beginning and at the end
  x <- consensusMatrix(x)
  x <- x[rownames(x)!=" ",,drop = FALSE]
  x <- paste(rownames(x)[apply(x, 2, which.max)], collapse = "")
  gsub("-","",x)
}
