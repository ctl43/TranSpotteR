#' @export
#' @importFrom Biostrings DNAStringSet
#' @importFrom ShortRead writeFasta

bwa_alignment <- function(seq, ref="/home/ctlaw/reference/Homo_sapiens/hs37d5/hs37d5_KJ173426.fa",
                         working_dir=NULL, samtools_param="-F 4", call_bwa = "bwa mem ", threads = 1){
  if(length(seq)==0){
    return(data.frame("QNAME"= character(), "FLAG"= integer(),
                      "RNAME"= character(), "POS"= integer(),
                      "MAPQ"= integer(), "CIGAR"= character(),
                      "SEQUENCE"= character()))
  }
  read_names <- names(seq)
  if(is.null(read_names)){
    names(seq) <- paste0("read_", seq_along(seq))
  }

  if(any(duplicated(read_names))){
    stop("Duplicated sequence name is detected.")
  }

  if (is.null(working_dir)) {
    tempdir()
    working_dir <- tempfile()
    dir.create(working_dir)
    on.exit(unlink(working_dir, recursive = TRUE))
  }
  if(samtools_param==""){
    call_samtools <- "|samtools view /dev/stdin"
  }else{
    call_samtools <- paste0("|samtools view ", samtools_param, " /dev/stdin")
  }
  temp_fa <- file.path(working_dir, basename(working_dir))
  writeFasta(DNAStringSet(seq), temp_fa)
  cmd <- paste0(call_bwa, "-t ",threads, " ",ref, " ", temp_fa, call_samtools,"|cut -f 1,2,3,4,5,6,10")
  result <- system(cmd, intern = TRUE)
  if(length(result)==0){
    return(data.frame(
      "QNAME"= character(), "FLAG" = integer(),
      "RNAME"= character(), "POS" = integer(),
      "MAPQ"= integer(), "CIGAR" = character(),
      "SEQUENCE"= character()))
  }else{
    bwa_aln <- read.delim(text = result, stringsAsFactors = FALSE, header = FALSE)
    colnames(bwa_aln) <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "SEQUENCE")
    bwa_aln$QNAME <- as.character(bwa_aln$QNAME)
    bwa_aln
  }

}
