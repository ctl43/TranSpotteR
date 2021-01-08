#' @export
#' @importFrom vroom vroom

get_header <- function(sam){
  N <- 1L
  curfile <- file(sam, open="r")
  collected <- c()
  repeat {
    curline <- readLines(curfile, n=1)
    if (!grepl("^@", curline)) {
      break
    }
    if (grepl("^@SQ", curline)) {
      collected[N] <- curline
    }
    N <- N + 1L
  }
  collected <- paste(collected, collapse = "\n")
  collected <- fread(collected,header = FALSE)
  collected$V2 <- sub("^SN:","", collected$V2)
  collected$V3 <- as.integer(sub("^LN:","", collected$V3))
  collected <- split(collected$V3, f=collected$V2)
  close(curfile)
  return(list(nheader = N - 1L, seq_info = unlist(collected)))
}

#' @export
read_sam <- function(sam, start = 1, end = NULL, select = c(1:6, 10),
                     nheader = NULL, header = NULL, return_header = FALSE){
  if(!(all(select%in%1:11))){
    stop()
  }

  if(is.null(nheader)|return_header){
    header <- get_header(sam)
    nheader <- header$nheader
  }

  sam_col_types  <- list(X1="c", X2="i", X3="c", X4="i", X5="i",
                         X6="c", X7="c", X8="i", X9="i",X10="c",
                         X11="c")
  col_types <- sam_col_types[select]

  skip <- nheader + start - 1
  if(is.null(end)){
    system.time(sam_readin <- vroom(sam, delim = "\t", skip = skip,
                                    col_names = FALSE,
                                    col_types = col_types,
                                    col_select = select,
                                    escape_double = FALSE, quote = ""))
  }else{
    nrow <- end - start + 1
    system.time(sam_readin <- vroom(sam, delim = "\t", skip = skip,
                                    col_names = FALSE, n_max = nrow,
                                    col_types = col_types,
                                    col_select = select,
                                    escape_double = FALSE, quote = ""))
  }
  sam_col_names <- c("QNAME","FLAG","RNAME","POS","MAPQ","CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL")
  colnames(sam_readin) <- sam_col_names[select]
  if(return_header){
    list(header = header$seq_info, mapping = sam_readin)
  }else{
    sam_readin
  }
}
