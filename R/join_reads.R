single_join_reads <- function(seq_1, seq_2, seed_size = min_overlap, min_overlap = 10, avoid_polyAT=FALSE){
  # seed_size is to increase the speed (larger 5 is better)
  # Should I keep seed size???
  # if(min_overlap < 5){
  #   stop("min_overlap must be larger than 5")
  # }
  
  if(length(seq_1)==0|length(seq_2)==0){
    return(DNAStringSet())
  }
  
  seq_1 <- DNAString(as.character(seq_1))
  seq_2 <- DNAString(as.character(seq_2))
  len_1 <- nchar(seq_1)
  len_2 <- nchar(seq_2)
  
  seed_start <- len_1-(seed_size-1)
  seed_seq <- as.character(subseq(seq_1, len_1-(seed_size-1), len_1))
  matches <- matchPattern(as.character(seed_seq), seq_2)
  match_start <- start(matches)
  if(avoid_polyAT){
    polyA <- paste(rep("A", seed_size), collapse = "")
    polyT <- paste(rep("T", seed_size), collapse = "")
    if(length(matches)==0|seed_seq==polyA|seed_seq==polyT){
      return(DNAStringSet())
    }
  }else{
    if(length(matches)==0){
      return(DNAStringSet())
    }
  }
  
  .internal_fun <- function(seed_start, match_start, seq_1, seq_2){
    current_1 <- seed_start
    current_2 <- match_start
    extend_1 <- Views(seq_1, IRanges(current_1, width=1))
    extend_2 <- Views(seq_2, IRanges(current_2, width=1))
    
    while(extend_1==extend_2){
      current_1 <- current_1-1
      current_2 <- current_2-1
      if(current_1==0|current_2==0){
        break
      }
      extend_1 <- Views(seq_1, IRanges(current_1, width=1))
      extend_2 <- Views(seq_2, IRanges(current_2, width=1))
    }
    ol_range_1 <- c(start(extend_1), len_1)
    ol_range_2 <- c(start(extend_2), match_start+(seed_size-1))
    list(seq1_ol_range=ol_range_1, seq2_ol_range=ol_range_2)
  }
  
  ol_range <- lapply(match_start, function(q).internal_fun(seed_start, q, seq_1, seq_2))
  
  # Determining which range to use
  temp_chosen <- which(sapply(ol_range, function(q)q$seq2_ol_range[1])==1)
  ol_len <- sapply(ol_range[temp_chosen], function(q)diff(q$seq2_ol_range))
  if(any(ol_len>min_overlap)){
    chosen <- which.max(ol_len)
    chosen <- temp_chosen[chosen]
    selected <- ol_range[[chosen]]
    # Combining sequence
    trimmed_seq2 <- Views(seq_2, IRanges(selected$seq2_ol_range[2]+1, len_2))
    DNAStringSet(c(seq_1, DNAString(as.character(trimmed_seq2))))
  }else{
    DNAStringSet()
  }
}


join_reads <- function(p, q, seed_size = min_overlap, min_overlap = 10){
  if(any(c(length(p), length(q))==0)){ # to prevent that p has no sequence
    return(DNAStringSet())
  }
  
  # Grid search all
  out <- lapply(p, function(x)lapply(q, function(y)single_join_reads(x, y, seed_size = min_overlap, 
                                                              min_overlap = min_overlap)))
  out <- Reduce("c",unlist(out))
  names(out) <- seq_along(out)
  out
}
