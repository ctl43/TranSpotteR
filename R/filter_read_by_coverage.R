filter_read_by_coverage <- function(x, min_read=2, min_peak=5, window_tol=500, rescue_tol=500, size_tol=0, rescue=FALSE){
  cvg <- coverage(x)
  sliced <- slice(cvg, min_read)
  max_count <- lapply(sliced, max)
  irange <- unlist(IRangesList(lapply(sliced, function(x)x@ranges)))
  collapsed <- GRanges(names(irange), irange, seqinfo=seqinfo(x))
  collapsed$max_count <- unlist(max_count)
  peaks <- collapsed[collapsed$max_count>=min_peak]
  
  if(rescue){
    temp <- collapsed[!collapsed$max_count>=min_peak] # Recusing reads that are around the main clusters
    rescued <- subsetByOverlaps(temp, peaks + rescue_tol)
    collapsed <- c(peaks, rescued)
  }else{
    collapsed <- peaks
  }
  
  collapsed <- csaw::mergeWindows(collapsed, tol=window_tol)$regions # Merge cluster within the tolerated range
  
  ol <- findOverlaps(x, collapsed)
  collapsed$count <- as.integer(table(ol@to))
  collapsed <- collapsed[width(collapsed)>size_tol] # small size usually indicate the overrepresentation of PCR duplicate
  
  # Identifying the peaks of reads clusters
  splitted <- split(collapsed, as.factor(seqnames(collapsed)))
  splitted <- mapply(function(x,y)Views(x,y@ranges),x=cvg, y=splitted)

  find_peak <- function(x){
    val <- runValue(x)
    len <- runLength(x)
    max_value <- max(val)
    max_loc <- which(val==max_value)
    max_loc <- max_loc[which.max(len[max_loc])]
    c(as.integer(len[max_loc]/2)+sum(len[1:(max_loc-1)]), max_value)
  }
  
  peaks <- lapply(splitted, function(x)t(sapply(as.vector(x), find_peak)))
  peaks <- peaks[lengths(peaks)>0]
  peaks <- do.call(rbind, peaks)
  collapsed$peaks <- peaks[,1]
  collapsed$peak_count <- peaks[,2]
  collapsed
}
