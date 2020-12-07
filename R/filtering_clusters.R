#' @export
#' @importFrom S4Vectors elementMetadata 'elementMetadata<-'
#' @importFrom BiocGenerics start end grepl
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges LogicalList

filter_clusters <- function(x, insert_seqname = "Hot_L1_polyA"){
  # As the inference process is complicated and takes time, this function is to filter those obviously useless information
  relationship_grp <- factor(rep(x$group, lengths(x$read_annotation)), levels = sort(unique(x$group)))
  grp <- factor(rep(seq_along(x), lengths(x$read_annotation)), levels = seq_along(x))
  flat_anno <- unname(unlist(x$read_annotation))

  gr <- convert_character2gr(flat_anno)
  not_insert_only <- sum(grepl(insert_seqname, seqnames(gr))) != lengths(gr) # reads with insert sequence has not anchor information are useless


  # Filtering related read cluster group that has not insert infromation
  insert_gr <- gr[grepl(insert_seqname, seqnames(gr))]
  has_insert <- any(start(insert_gr) < 6000)
  group_has_insert <- any(LogicalList(split(has_insert, relationship_grp)))
  is_usable_grp <- sort(unique(x$group))[group_has_insert]
  is_from_usable_grp <- relationship_grp%in%is_usable_grp
  # good_insert <- any(width(insert_gr) > 100) # Insert length is longer than 100 nt

  # Extracting regions that are close together (does not indicate transduction event)
  anchor_gr <- gr[!grepl(insert_seqname, seqnames(gr))]
  gr_range <- range(start(anchor_gr))
  gr_range[,1][abs(gr_range[,1]) > 999999999] <- 0
  gr_range[,2][abs(gr_range[,2]) > 999999999] <- 0
  range_ok <- (gr_range[,2] - gr_range[,1]) > 50000
  chr_ok <- lengths(unique(seqnames(anchor_gr))) > 1

  # Getting regions has more than just genomic sequences
  tmp_anno <- CharacterList(strsplit(flat_anno, " "))
  tmp_anno <- tmp_anno[!grepl(":|NNNNN", tmp_anno)]
  has_seq <- sum(nchar(tmp_anno)) > 3
  orphan_ok <- lengths(anchor_gr) >= 1 & (has_seq | lengths(insert_gr) > 0) # Only contains genomic regions but with sequences (could possibly poly-A)
  region_ok <- lengths(anchor_gr) > 1 # Containing more than one 1 genomic regions

  # Combining all filters
  # is_informative <- not_insert_only & (orphan_ok | region_ok | range_ok | chr_ok | good_insert) & is_from_usable_grp
  is_informative <- not_insert_only & (orphan_ok | (region_ok & range_ok) | chr_ok) & is_from_usable_grp
  informative_anno <- split(flat_anno[is_informative], grp[is_informative])
  elementMetadata(x)$informative_anno <- informative_anno
  elementMetadata(x)$is_usable <- lengths(x$informative_anno) > 0
  x
}
