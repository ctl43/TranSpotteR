#' @export
#' @importFrom data.table data.table setDT
#' @importFrom Biostrings consensusMatrix DNAStringSet
#' @importFrom uuid UUIDgenerate

#' @export
greedy_scs <- function(vec, n_reads = NULL, msa_result = FALSE,
                         return_no_assembly = FALSE, add_id = TRUE, min_len = 8L, min_pid = 85,
                         run_msa = TRUE)
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
  # vec <- sim_frag
  if(is.null(n_reads)){
    n_reads <- rep(1L, length(vec))
  }

  idx <- names(n_reads) <- names(vec) <- seq_along(vec)
  original_vec <- vec

  info <- data.table(t(combn(names(vec), 2)))
  colnames(info) <- c("seq1", "seq2")
  info <- data.table(info, setDT(overlapper(vec[info[[1]]], vec[info[[2]]])))
  info <- .reorder_seq_info(info)
  # Filtering abberent alignment
  info <- info[info$pid > min_pid & info$aln_ol_len > min_len & (info$seq1_ol_start == 0L | info$seq2_ol_start == 0L), ]
  info <- info[order(info$aln_ol_len, decreasing = TRUE), ]

  # Converting to 1-based index
  info$seq1_ol_start <- info$seq1_ol_start + 1L
  info$seq1_ol_end <- info$seq1_ol_end + 1L
  info$seq2_ol_start <- info$seq2_ol_start + 1L
  info$seq2_ol_end <- info$seq2_ol_end + 1L
  setkeyv(info, c("seq1", "seq2"))

  while(nrow(info) > 0){
    selected <- info[1L, ]
    left_idx <- selected$seq1
    right_idx <- selected$seq2
    new_seq <- paste0(substr(vec[left_idx], 1L, selected$seq1_ol_start - 1L),
                      substr(vec[right_idx], selected$seq2_ol_start, selected$seq2_len))
    names(new_seq) <- paste0(selected$seq1, "_", selected$seq2)
    vec <- vec[!names(vec) %in% c(selected$seq1, selected$seq2)]
    ol <- overlapper(rep(new_seq, length(vec)), vec)
    suppressWarnings(ol <- data.table(data.table(seq1 = names(new_seq), seq2 = names(vec)), ol))
    ol <- ol[ol$pid > min_pid & ol$aln_ol_len > min_len & (ol$seq1_ol_start == 0L | ol$seq2_ol_start == 0L), ]
    ol <- .reorder_seq_info(ol)
    ol$seq1_ol_start <- ol$seq1_ol_start + 1L
    ol$seq1_ol_end <- ol$seq1_ol_end + 1L
    ol$seq2_ol_start <- ol$seq2_ol_start + 1L
    ol$seq2_ol_end <- ol$seq2_ol_end + 1L
    vec <- c(new_seq, vec)
    info <- info[!(seq1 %in% c(left_idx, right_idx) | seq2 %in% c(left_idx, right_idx))]
    info <- rbindlist(list(ol, info))
    info <- info[order(info$aln_ol_len, decreasing = TRUE), ]
  }

  # Picking out those read that did not assemble
  ass_id <- names(vec)
  solo_id <- ass_id[!grepl("_", ass_id)]
  non_solo_id <- grep("_", ass_id, value = TRUE)
  member_idx <- strsplit(non_solo_id, "_")

  if(run_msa){
    msa_out <- mapply(function(x,y)t_coffee(c(original_vec[x],y)),
                      x = member_idx, y = vec[non_solo_id], SIMPLIFY = FALSE)
    consensus <- sapply(msa_out, .process_msa)
    names(msa_out) <- non_solo_id
  }else{
    consensus <- vec
  }
  names(consensus) <- non_solo_id
  consensus <- unlist(consensus)
  if(add_id & length(consensus) > 0L){
    names(consensus) <- UUIDgenerate(n = length(consensus))
  }

  # Computing the number of reads consisting the consensus sequence
  non_solo_n_reads <- unlist(lapply(member_idx, function(x)sum(n_reads[x])))

  if(return_no_assembly){
    solo <- CharacterList(vec[solo_id])
    n_reads <- c(non_solo_n_reads, n_reads[solo_id])
    consensus <- c(consensus, solo)
    msa_out <- c(msa_out, split(solo, solo_id))
  }else{
    n_reads <- non_solo_n_reads
  }

  if(msa_result){
    return(list(consensus = consensus, msa = msa_out))
  }else{
    return(list(consensus, n_reads))
  }
}

# Function to determine a pair of sequences that which one is of the left hand side
.reorder_seq_info <- function(info){
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
  ordered_info
}

.process_msa <- function(x){
  x <- replcae_space(x) # replacing the space at the beginning and at the end
  x <- consensusMatrix(x)
  x <- x[rownames(x)!=" ",,drop = FALSE]
  x <- paste(rownames(x)[apply(x, 2, which.max)], collapse = "")
  gsub("-","",x)
}
