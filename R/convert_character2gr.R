#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors split
#' @importFrom IRanges IRanges


# Converting read annotation to genomic range
convert_character2gr <- function(y){
  y <- y[grepl(":", y)]
  grp <- factor(rep(seq_along(y), lengths(y)), levels = seq_along(y))
  y <- unlist(y)
  y <- do.call(rbind, strsplit(y, ":|-"))
  y[,4][y[,4]==""] <- "-" # adding back the negative strand
  y <- data.frame(y)
  y[,2] <- as.integer(y[,2])
  y[,3] <- as.integer(y[,3])
  gr <- GRanges(seqnames=y[,1], IRanges(start=y[,2], end=y[,3]), strand = y[,4])
  split(gr, grp)
}
