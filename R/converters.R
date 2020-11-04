sam2tb <- function(sam, samtools="samtools", include_unmapped=TRUE){
  # getting header
  header <- system(paste0(samtools, " view -H ", sam), intern=TRUE)
  header <- fread(paste(grep("^@SQ",header, value=TRUE), collapse = "\n"), header=FALSE)
  seq_names <- gsub("^SN:","", header[[2]])
  seq_lens <- as.integer(gsub("^LN:","", header[[3]]))
  seq_info <- Seqinfo(seq_names, seq_lens)

  #reading in the sam files
  txt <- paste0("samtools view ", sam, "|cut -f 1,2,3,4,5,6,10")
  txt <- paste(system(txt, intern=TRUE), collapse = "\n")
  what <- c("character","integer","character", "integer", "integer", "character","character")
  x <- fread(txt, colClasses = what, col.names = c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "SEQUENCE"))
  
  if(include_unmapped){
    seq_info@seqnames <- c(seq_info@seqnames, "UNMAPPED")
    seq_info@seqlengths <- c(seq_info@seqlengths, 99999L)
    seq_info@is_circular <- c(seq_info@is_circular, NA)
    seq_info@genome <- c(seq_info@genome, NA)
    x[x$CIGAR=="*"]$RNAME <- "UNMAPPED"
    x[x$CIGAR=="*"]$POS <- 1
    x[x$CIGAR=="*"]$CIGAR <- paste0(nchar(x[x$CIGAR=="*"]$SEQUENCE),"M")
  }
  x
  

}

tb2gr <- function(x, seqinfo=NULL,  mapq_filter=NULL){
  if(is.null(x)){
    return(GRanges())
  }
  gr <- GenomicRanges::GRanges(seqnames=x$RNAME, 
                               IRanges(x$POS, 
                                       width=cigarWidthAlongReferenceSpace(x$CIGAR)),
                               strand=ifelse(bitwAnd(x$FLAG, 0x10), "-", "+"), 
                               seqinfo=seqinfo)
  names(gr) <- x$QNAME
  elementMetadata(gr) <- x
  
  # Remove pairs with low quality
  if(!is.null(mapq_filter)){
    discard <- gr$QNAME[gr$MAPQ<mapq_filter]
    gr <- gr[!gr$QNAME%in%discard,]
    gr
  }else{
    gr
  }
}

sam2gr <- function(x){
  txt <- sam2tb(x)
  gr <- tb2gr(x)
}
