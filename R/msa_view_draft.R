query_to_reference_position <- function(raln, qaln){
  r_pos <- 1L
  aln_len <- nchar(raln)
  q_in_rpos <- integer(aln_len)
  # stored_pos <- c()
  for(i in seq_len(aln_len)){
    if(raln[i] != "-"){
      q_in_rpos[i] <- r_pos
      r_pos <- r_pos + 1
    }else{
      q_in_rpos[i] <- r_pos
    }
  }
  q_in_rpos
}

convert_gap_to_ref_view <- function(y){
  if(length(grep("-", y)) == 0L){
    return(seq_along(y[[1]]))
  }
  original <- BString(gsub("-", "", y[[1]]))
  ref_len <- nchar(original)
  out <- gap_finder(y)
  idx <- order(out[[2]], decreasing = TRUE)
  out[[2]] <- out[[2]][idx]
  out[[1]] <- out[[1]][idx]
  idx <- order(out[[1]])
  out[[1]] <- out[[1]][idx]
  out[[2]] <- out[[2]][idx]
  is_dup <- duplicated(out[[1]])
  out[[1]] <- out[[1]][!is_dup]
  out[[2]] <- out[[2]][!is_dup]
  n_gap <- length(out[[1]])
  
  gapped_aln <- integer()
  gap_info_counter <- 1L
  
  n_loop <- nchar(original) + 1
  
  for(i in seq_len(n_loop)){
    
    cur_gap_before <- out$gap_before[gap_info_counter]
    
    if(i == cur_gap_before){
      cur_gap_len <- out$gap_count[gap_info_counter]
      
      if(i == ref_len + 1){ # if it is the right clipped len
        cur_n_gap <- seq_len(cur_gap_len)
      }else{
        cur_n_gap <- seq_len(cur_gap_len + 1)
      }
      
      for(j in cur_n_gap){
        gapped_aln <- c(gapped_aln, i)
      }
      
      if(gap_info_counter != n_gap){
        gap_info_counter <- gap_info_counter + 1L
      }
      
    }else{
      
      if(i != n_loop){
        gapped_aln <- c(gapped_aln, i)
      }
      
    }
  }
  return(gapped_aln)
}

gap_finder <- function(y){
  gap_pos_storage <- gap_count_storage <- integer()
  for(x in seq_along(y)){
    ref <- y[[x]]
    gap_pos <- integer()
    gap_counter <- 0L
    ungapped_counter <- 0L
    previous_is_gap <- FALSE
    for(i in 1:nchar(ref)){
      if(ref[i] != "-"){
        if(previous_is_gap){
          gap_pos_storage <- c(gap_pos_storage, gap_pos)
          gap_count_storage <- c(gap_count_storage, gap_counter)
          gap_pos <- integer()
          gap_counter <- 0L
        }
        ungapped_counter <- ungapped_counter + 1L
        previous_is_gap <- FALSE
      }else{
        if(previous_is_gap){
          gap_counter <- gap_counter + 1L
          previous_is_gap <- TRUE
        }else{
          gap_counter <- gap_counter + 1L
          gap_pos <- ungapped_counter + 1L
          previous_is_gap <- TRUE
        }
      }
    }
    
    #For the last bit in case it has clipped region at the end
    if(previous_is_gap){
      gap_pos_storage <- c(gap_pos_storage, gap_pos)
      gap_count_storage <- c(gap_count_storage, gap_counter)
    }
  }
  return(list(gap_before = gap_pos_storage, gap_count = gap_count_storage))
}

translate_to_extended <- function(cur_q, cur_q_aln, extended){
  # cur_q <- translated_pos[[2]]
  # cur_q_aln <- qaln[[2]]
  q_counter <- 1L
  extended_q <- BString(paste(rep("X", length(extended)), collapse = ""))
  cur_r_base <- extended[1]
  cur_q_base <- cur_q[q_counter]
  n_q <- length(cur_q)
  for(i in seq_len(length(extended) - 1)){
    if(q_counter != n_q){
      next_q_base <- cur_q[q_counter + 1]
      next_r_base <- extended[i + 1]
      if(cur_r_base == cur_q_base & next_r_base == next_q_base){
        extended_q[i] <- cur_q_aln[q_counter]
        extended_q[i + 1] <- cur_q_aln[q_counter + 1]
        q_counter <- q_counter + 1
        cur_q_base <- next_q_base
      }else{
        extended_q[i] <- "-"
      }
    }else{
      extended_q[i + 1L] <- "-"
    }
    cur_r_base <- next_r_base
  }
  extended_q
}

combined_aln <- function(raln, qaln){}
  extended <- convert_gap_to_ref_view(raln)
  translated_pos <- mapply(query_to_reference_position, raln = raln, qaln = qaln, SIMPLIFY = FALSE)
  translated <- mapply(function(x,y)translate_to_extended(x, y, extended = extended), x = translated_pos, y = qaln)
  translated
}