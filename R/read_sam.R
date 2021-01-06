#' @export
#' @importFrom data.table fread
#' @importFrom GenomeInfoDb Seqinfo

get_headers <- function(sam, return_seqinfo = FALSE, n_rows = 100){
  # Getting headers
  skip <- 0
  last <- 1
  while(length(last) > 0){
    keep <- last
    readin <- fread(sam, sep = "", skip = skip, nrows = n_rows)[[1]]
    last <- tail(grep("^@", readin), 1)
    skip <- skip + n_rows
  }
  nline <- skip - (n_rows*2) + keep
  headers <- fread(sam, sep = "", skip = 0, nrows = nline)[[1]]
  seqinfo <- fread(text = grep("^@SQ", headers, value=TRUE), header = FALSE)[,2:3]
  chr <- sub("^SN:", "", seqinfo[[1]])
  chr_size <- as.integer(sub("^LN:", "", seqinfo[[2]]))
  if(return_seqinfo){
    seqinfo <- Seqinfo(chr, chr_size)
    return(list(seqinfo, nline = nline))
  }else{
    return(nline)
  }
}

#' @export
read_sam <- function(sam_file, start, end, select = c(1:6, 10), nheader = NULL, nfields = 30){
  if(is.null(nheader)){
    nheaders <- get_headers(sam_file)
  }
  skip <- start + nheaders - 1
  nrow <- end - start + 1
  sam <- fread(sam_file, sep = "", skip = skip, nrow = nrow, select = list("character" = 1))[[1]]
  # nfields <- max(lengths(strsplit(sam, "\t")), na.rm = TRUE)
  col_names <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "SEQUENCE")
  if(is.null(nfields)){
    nfields <- max(count.fields(sam_file, skip = skip, sep = "\t", comment.char = "", quote = ""), na.rm = TRUE)
  }
  # colClasses <- c("character","integer","character", "integer", "integer", "character","character")
  out <- fread(text = c(paste(1 : nfields, collapse = "\t"), sam),
               header = TRUE, fill = TRUE, select = select, nThread = 1)#, colClasses = colClasses)#, skip = skip, nrow = nrow)
  names(out) <- col_names
  out
}
