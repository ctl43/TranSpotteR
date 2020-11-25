#' @export
#' @importFrom IRanges IntegerList LogicalList CharacterList
#' @importFrom S4Vectors 'elementMetadata<-' elementMetadata split
#' @importFrom BiocGenerics start end strand
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocParallel bplapply MulticoreParam

# A wrapper for read annotation
annotate_constructed_reads <- function(clusters, partner_is_anchor = FALSE,  BPPARAM = MulticoreParam(workers = 3)){
  cluster_anno <- .annotate_reads(unlist(clusters$cluster_contigs), BPPARAM = BPPARAM)
  partner_anno <- .annotate_reads(unlist(clusters$partner_contigs), BPPARAM = BPPARAM)
  long_anno <- .annotate_reads(unlist(clusters$long_contigs), BPPARAM = BPPARAM)

  # Tidying up the data
  .tidyup <- function(anno, clusters, prefix = ""){
    grp <- lapply(anno, function(x)factor(gsub("\\..*","",names(x)), levels = names(clusters)))
    anno <- mapply(split,x=anno, f=grp)
    anno$left_clipped_length_in_aln_1 <- IntegerList(anno$left_clipped_length_in_aln_1)
    anno$unmapped_length_in_aln_1 <- IntegerList(anno$unmapped_length_in_aln_1)
    anno$right_clipped_length_in_aln_1 <- IntegerList(anno$right_clipped_length_in_aln_1)
    names(anno) <- paste0(prefix, "_", names(anno))
    for(i in names(anno)){
      elementMetadata(clusters)[[i]] <- anno[[i]]
    }
    clusters
  }
  clusters <- .tidyup(cluster_anno, prefix = "anchor", clusters=clusters)
  clusters <- .tidyup(partner_anno, prefix = "partner", clusters=clusters)
  clusters <- .tidyup(long_anno, prefix = "long", clusters=clusters)

  # Annotating read cluster with anchor regions
  info <- elementMetadata(clusters)
  if(partner_is_anchor){
    mapping2genome <- info[grep("partner.*mapping2genome|anchor.*mapping2genome|long.*mapping2genome", colnames(info))]
  }else{
    mapping2genome <- info[grep("anchor.*mapping2genome|long.*mapping2genome", colnames(info))]
  }

  left_clipped_length <- info[grep("left_clipped_length", colnames(info))]
  right_clipped_length <- info[grep("right_clipped_length", colnames(info))]
  clusters$has_long <- lengths(clusters$long_contigs) != 0

  # Combining the anchor regions
  combined_anchor <- unlist(mapping2genome, use.names = FALSE)
  combined_anchor <- unlist(combined_anchor, use.names = FALSE)
  combined_anchor <- split(combined_anchor, factor(gsub("\\..*","", names(combined_anchor)),
                                                   levels = rownames(mapping2genome)))
  clusters$combined_anchor <- combined_anchor

  # Checking whether the anchor regions are normal
  combined_anchor <- unlist(clusters$combined_anchor, use.names = FALSE)
  grp <- gsub("\\..*", "", names(combined_anchor))
  expanded <- clusters[grp]
  same_chr <- as.character(seqnames(expanded)) == as.character(seqnames(combined_anchor))
  start_dist <- start(expanded) - start(combined_anchor)
  abnormal <- start_dist > 50000 | !same_chr
  is_abnormal_anchor <- LogicalList(split(abnormal, factor(grp, levels = names(clusters))))
  clusters$is_abnormal_anchor <- is_abnormal_anchor

  # Combining bait regions
  mapping2bait <- info[grep("mapping2bait", colnames(info))]
  combined_bait <- unlist(mapping2bait, use.names = FALSE)
  combined_bait <- unlist(combined_bait, use.names = FALSE)
  combined_bait <- split(combined_bait, factor(gsub("\\..*","", names(combined_bait)), levels = rownames(info)))
  clusters$combined_bait <- combined_bait

  # Combining annotation
  anchor_anno <- unlist(clusters$anchor_read_annotation)
  partner_anno <- unlist(clusters$partner_read_annotation)
  total_name <- unique(c(names(anchor_anno), names(partner_anno)))
  combined_storage <- partner_storage <- anchor_storage <- rep(CharacterList(c()), length(total_name))
  names(combined_storage) <- names(partner_storage) <- names(anchor_storage) <- total_name
  anchor_storage[names(anchor_anno)] <- anchor_anno
  partner_storage[names(partner_anno)] <- partner_anno
  anchor_mapped <- lengths(anchor_storage) > 0
  partner_mapped <- lengths(partner_storage) > 0
  both_mapped <- anchor_mapped&partner_mapped
  combined_storage[!both_mapped&anchor_mapped] <- anchor_storage[(!both_mapped) & anchor_mapped]
  combined_storage[(!both_mapped)&partner_mapped] <- partner_storage[(!both_mapped) & partner_mapped]
  strand <- as.character(strand(clusters)[1])
  if(strand=="+"){
    combined_storage[both_mapped] <- (mapply(function(x,y)c(x,"NNNNN",y), x = anchor_storage[both_mapped],
                                             y = partner_storage[both_mapped], SIMPLIFY = FALSE))
  }else{
    combined_storage[both_mapped] <- (mapply(function(x,y)c(y,"NNNNN",x), x = anchor_storage[both_mapped],
                                             y = partner_storage[both_mapped], SIMPLIFY = FALSE))
  }
  names(combined_storage) <- sub("^[0-9]+\\.","", names(combined_storage))
  long_read_anno <- unlist(clusters$long_read_annotation)
  names(long_read_anno) <- sub("^[0-9]+\\.","", names(long_read_anno))
  origin <- factor(sub("\\..*","", c(names(long_read_anno), names(combined_storage))), levels = seq_along(clusters))
  clusters$combined_read_anno <- split(c(long_read_anno, combined_storage), origin)
  clusters
}

#' @importFrom BiocGenerics cbind do.call
.combined_rle <- function(x){
  if(length(x)==1){
    return(x[[1]])
  }
  x <- do.call(cbind, x)
  result <- x[,2]
  for(i in 3:ncol(x)){
    result[result=="S"] <- x[,i][result=="S"]
  }
  out <- Rle(result, x[,1])
  out
}

#' @importFrom GenomicAlignments cigarWidthAlongQuerySpace
#' @importFrom S4Vectors split runValue runLength nchar
#' @importFrom IRanges RleList CharacterList
#' @importFrom Biostrings reverseComplement BStringSet DNAStringSet
#' @importFrom BiocGenerics paste

.annotate_reads <- function(seq, ref = "~/dicky/reference/fasta/line1_reference/hot_L1_polyA.fa", BPPARAM = MulticoreParam(workers = 3)){
  if(length(seq) == 0){
    return(NULL)
  }
  # Mapping to LINE1 sequence
  aln_1 <- bwa_alignment(seq, ref = ref, samtools_param = "")
  aln_1 <- convertingHtoS(aln_1, unique_id = "QNAME")
  mapping_1 <- sam2gr(aln_1[aln_1$CIGAR!="*",])

  # Annotating the reads
  aln_1$unified_cigar <- unify_cigar_strand(aln_1$CIGAR, flag = aln_1$FLAG, to = "+", along_query=TRUE)

  # Getting clipped sequences
  clipped_seq_1 <- get_unmapped_clipped_read(aln_1, include_middle_unmapped = TRUE)

  # Mapping to genome
  middle <- clipped_seq_1$middle
  clipped_seq_1 <- clipped_seq_1$clipped

  aln_2 <- bplapply(clipped_seq_1, bwa_alignment, call_bwa = "bwa mem ", samtools_param = "-F 128 -F 4", BPPARAM = BPPARAM)
  aln_2 <- lapply(aln_2, convertingHtoS, unique_id = "QNAME")
  aln_2 <- lapply(aln_2, function(x)x[x$MAPQ>10,])
  mapping_2 <- lapply(aln_2, sam2gr)

  middle_aln_2 <- bwa_alignment(unlist(middle$seq), call_bwa = "bwa mem ", samtools_param = "-F 128 -F 4")
  middle_grp <- sub("\\.[0-9]+$","\\1",middle_aln_2$QNAME)
  middle_aln_2 <- convertingHtoS(middle_aln_2, unique_id = "QNAME")
  mid_mapping <- sam2gr(middle_aln_2)
  middle_coordinate <- as.character(mid_mapping)
  mid_start <- unlist(middle$start)[names(middle_coordinate)]
  mid_end <- unlist(middle$end)[names(middle_coordinate)]

  # Annotating the sequencing reads
  anno <- aln_1[!duplicated(aln_1$QNAME),]
  anno_seq <- DNAStringSet(anno$SEQUENCE)
  is_rc <- ifelse(bitwAnd(anno$FLAG, 0x10), "-", "+") == "-"
  anno_seq[is_rc] <- reverseComplement(DNAStringSet(anno$SEQUENCE[is_rc]))
  anno_seq <- RleList(strsplit(as.character(anno_seq), ""))
  names(anno_seq) <- anno$QNAME
  agreed_s <- anno_seq

  aln_1_q_pos <- .get_M_pos(aln_1$unified_cigar[aln_1$unified_cigar != "*"])
  match_start <- split(aln_1_q_pos$start, aln_1$QNAME[aln_1$unified_cigar != "*"])
  match_end <- split(aln_1_q_pos$end, aln_1$QNAME[aln_1$unified_cigar != "*"])
  if(length(mapping_1)!=0){
    aln_1_coordinate <- split(paste0(" ", as.character(mapping_1), " "), names(mapping_1))
    agreed_s[names(aln_1_coordinate)] <- RleList(mapply(function(q, x, y, z){
      for(i in seq_along(z)){
        q[x[i]:y[i]] <- z[i]
      }
      q
    }, q = agreed_s[names(aln_1_coordinate)], x = match_start, y = match_end, z = aln_1_coordinate))
  }

  for(i in c("left", "right", "unmapped")){
    current <- aln_2[[i]]
    if(nrow(current) == 0){
      next
    }
    unified_cigar <- unify_cigar_strand(current$CIGAR, flag=current$FLAG, to = "+", along_query = TRUE)
    genomic_loc <- paste0(" ", as.character(sam2gr(current)), " ")
    cigar_width <- cigarWidthAlongQuerySpace(unified_cigar)
    cigar_rle <- cigarToRleList(unified_cigar)
    cigar_rle <- mapply(function(x,y){x[x == "M"] <- y;x}, x = cigar_rle, y = genomic_loc)
    cigar_rle <- mapply(function(x,y){x[x == "I"] <- y;x}, x = cigar_rle, y = genomic_loc)
    cigar_rle <- split(cigar_rle, current$QNAME)
    cigar_rle <- lapply(cigar_rle, .combined_rle)

    if(i=="left"){
      agreed_s[current$QNAME] <- mapply(function(x,y,z){x[1:y] <- z;x}, x=agreed_s[current$QNAME], y = cigar_width, z = cigar_rle[current$QNAME])
    }

    if(i=="right"){
      l <- sum(runLength(agreed_s[current$QNAME]))
      agreed_s[current$QNAME] <- mapply(function(x,y,z,l){
        x[(l-y+1):l] <- z;x
      }, x=agreed_s[current$QNAME], y = cigar_width, z = cigar_rle[current$QNAME], l = l)
    }

    if(i=="unmapped"){
      agreed_s[current$QNAME] <- cigar_rle[current$QNAME]
    }
  }

  # for middle alignment
  mid_grp <- sub("\\.[0-9]+$", "", names(middle_coordinate))
  middle_coordinate <- split(middle_coordinate, mid_grp)
  mid_start <- split(mid_start, mid_grp)
  mid_end <- split(mid_end, mid_grp)
  agreed_s[names(middle_coordinate)] <- mapply(function(q, x, y, z){
    for(i in seq_along(z)){
      q[x[i]:y[i]] <- z[i]
    }
    q
  }, q=agreed_s[names(middle_coordinate)], x = mid_start, y = mid_end, z = middle_coordinate)
  agreed_s <- mapply(function(x,y){x[x=="S"] <- y[x=="S"]; x}, x = agreed_s, y=anno_seq)
  # agreed_s <- mapply(function(x,y){x[x=="I"] <- y[x=="I"]; x}, x = agreed_s, y=anno_seq)
  agreed_s <- RleList(agreed_s)

  # tidying up the annotated seq
  rv <- runValue(agreed_s)
  rl <- runLength(agreed_s)
  seq_grp <- rep(seq_along(rv), lengths(rv))
  rv <- unlist(rv, use.names = FALSE)
  rl <- unlist(rl, use.names = FALSE)
  is_seq <- !grepl(":", rv)
  sub_seq <- BStringSet(mapply(rep, x = BStringSet(rv[is_seq]), times = rl[is_seq]))
  rv[is_seq] <- sub_seq
  anno_out <- CharacterList(split(rv, seq_grp))
  anno_out <- paste(anno_out, collapse = "")
  anno_out <- sub("^ ", "", anno_out)
  anno_out <- CharacterList(strsplit(anno_out, " +"))
  names(anno_out) <- names(seq)

  clipped_l <- sapply(clipped_seq_1, nchar, simplify = FALSE)
  names(clipped_l$right) <-  names(clipped_l$left) <- names(clipped_l$unmapped) <- names(clipped_seq_1[[1]])
  list(mapping2bait = mapping_1, unmapped_mapping2genome = mapping_2$unmapped,
       left_mapping2genome = mapping_2$left,
       right_mapping2genome = mapping_2$right,
       mid_mapping2genome = mid_mapping,
       read_annotation = anno_out,
       left_clipped_length_in_aln_1 = clipped_l$left,
       right_clipped_length_in_aln_1 = clipped_l$right,
       unmapped_length_in_aln_1 = clipped_l$unmapped
  )
}



