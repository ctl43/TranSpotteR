# This method has a problem, it may discard read that are useful
read_cap3_consensus <- function(cap3_out){
  # Reading in cap3 format output
  
  detail_loc <- grep("DETAILED DISPLAY OF CONTIGS", cap3_out)
  cap3_detail <- cap3_out[-(1:detail_loc)]
  
  if(detail_loc==length(cap3_out)){ # No contig
    return(DNAStringSet())
  }
  
  start <- grep("Contig", cap3_detail)
  if(length(start)==1){ # One contig only
    end <- length(cap3_detail)
  }else{ # More than one contig
    end <- c(start[-1]-1, length(cap3_detail))
  }
  
  splitted_detail <- mapply(function(x,y)cap3_detail[x:y], x=start, y=end, SIMPLIFY = FALSE)
  .internal_process <- function(x){
    paste0(gsub("consensus             ", "", x[grep("consensus", x)]), collapse = "")
  }
  
  out <- sapply(splitted_detail,.internal_process)
  out <- DNAStringSet(gsub("-","", out))
  names(out) <- seq_along(out)
  out
}

read_velvet_out <- function(velvet_out){ # velvet_out is a fasta file
  velvet_out <- readFasta(velvet_out)
  if(length(velvet_out)==0){
    return(velvet_out@sread)
  }else{
    # names(velvet_out@sread) <- paste0("velvet_", 1:length(velvet_out))
    names(velvet_out@sread) <- seq_along(velvet_out)
    return(velvet_out@sread)
  }
}

.assembly <- function(seq, working_dir=NULL,
                      #cap3_param = "-a 11 -b 16 -c 6 -d 21 -e 11 -f 2 -g 1 -h 3 -i 21 -j 31 -k 0 -m 1 -n -1 -o 16 -p 66 -r 0 -s 251 -t 31 -u 1 -v 1 -y 6 -z 1",
                      cap3_param = "-i 21 -j 31 -o 16 -s 251",
                      velvet_dir="~/tools/velvet-1.2.10", 
                      cap3_dir="~/tools/CAP3/",
                      method="both", joining_reads=TRUE){
  # Combining OLC and DBG graph as DBG method does not include full polyA into the final assembly
  if (is.null(working_dir)) {
    tempdir()
    working_dir <- tempfile()
    dir.create(working_dir)
    on.exit(unlink(working_dir, recursive = TRUE))
  }
  temp_fa1 <- file.path(working_dir, basename(working_dir))
  writeFasta(DNAStringSet(seq), temp_fa1)
  
  # Velvet
  if(method=="both"|method=="velvet"){
    system(paste(file.path(velvet_dir, "velveth"), working_dir, "21", temp_fa1, collapse = " "))
    system(paste(file.path(velvet_dir, "velvetg"),working_dir, "-min_contig_lgth 1", collapse = " "))
    velvet_contigs <- read_velvet_out(file.path(working_dir, "contigs.fa"))
  }else{
    velvet_contigs <- DNAStringSet()
  }
  
  # CAP3
  if(method=="both"|method=="cap3"){
    cap3_contigs <- system(paste(file.path(cap3_dir, "cap3"), temp_fa1, cap3_param, collapse = " "), intern=TRUE)
    cap3_contigs <- read_cap3_consensus(cap3_contigs)
    # cap3_contigs <- DNAStringSet(gsub("-","", cap3_contigs))
  }else{
    cap3_contigs <- DNAStringSet()
  }
  
  if(method=="cap3"){
    return(cap3_contigs)
  }
  cap3_is_null <- length(cap3_contigs)==0
  velvet_is_null <- length(velvet_contigs)==0
  
  # Only contig from CAP3
  if(!cap3_is_null&velvet_is_null){
    both_contigs <- c(cap3=cap3_contigs)
    final_contig <- cap3_contigs
  }
  
  # Only contig from velvet
  if(cap3_is_null&!velvet_is_null){
    both_contigs <- c(velvet=velvet_contigs)
    final_contig <- velvet_contigs
  }
  
  # No contig from both
  if(cap3_is_null&velvet_is_null){
    both_contigs <- DNAStringSet()
    final_contig <- DNAStringSet()
  }
  
  # Final CAP3 to combining both velvet and cap3 output
  if(!cap3_is_null&!velvet_is_null){
    both_contigs <- c(velvet_contigs, cap3_contigs)
    temp_fa2 <- basename(tempfile(tmpdir=""))
    temp_fa2 <- file.path(working_dir, temp_fa2)
    writeFasta(both_contigs, temp_fa2)
    cap3_final_contig <- system(paste(file.path(cap3_dir, "cap3"), temp_fa2, cap3_param, collapse = " "), intern=TRUE)
    cap3_final_contig <- read_cap3_consensus(cap3_final_contig)
    has_cap3_final <- length(cap3_final_contig)>0

    if(has_cap3_final){
      final_contig <- cap3_final_contig
    }else{
      system(paste(file.path(velvet_dir, "velveth"), working_dir, "21", temp_fa1, collapse = " "))
      system(paste(file.path(velvet_dir, "velvetg"),working_dir, "-min_contig_lgth 1", collapse = " "))
      velvet_final_contigs <- read_velvet_out(file.path(working_dir, "contigs.fa"))
      final_contig <- velvet_final_contigs
    }

    if(length(final_contig)==0){
      final_contig <- DNAStringSet()
    }
    if(joining_reads){
      if(length(final_contig)==2){
        test_1 <- join_reads(final_contig[1], final_contig[2])
        test_2 <- join_reads(final_contig[2], final_contig[1])
        test_1_len <- nchar(test_1)
        test_2_len <- nchar(test_2)
        
        if(length(test_1_len)==0){
          test_1_len <- 0
        }
        
        if(length(test_2_len)==0){
          test_2_len <- 0
        }
        
        if(test_1_len!=test_2_len){
          if(test_1_len>test_2_len){
            final_contig <- test_1
          }else{
            final_contig <- test_2
          }
        }
        
      }
    }
  }
  unlink(working_dir)
  list(both_contigs, final_contig)
}
