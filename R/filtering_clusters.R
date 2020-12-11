#' @export
#' @importFrom S4Vectors elementMetadata 'elementMetadata<-'
#' @importFrom BiocGenerics start end grepl
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges LogicalList

filter_clusters <- function(x, insert_seqname = "Hot_L1_polyA")
  # As the inference process is complicated and takes time,
  # this function is to filter those obviously useless information to reduce the computation time
  # Written by Cheuk-Ting Law
{
  # Getting the mapping information
  anno <- x$read_annotation
  anno_grp <- factor(rep(seq_along(x), lengths(anno)), levels = seq_along(x))
  flat_anno <- unlist(anno, recursive = FALSE, use.name = FALSE)
  anno <- CharacterList(lapply(flat_anno, "[[", i = "anno"))

  # Converting into a list of GRanges objects
  anno_ranges <- anno[S4Vectors::grepl(":", anno)]
  gr <- convert_character2gr(anno_ranges)
  gr <- unique(gr)

  # 1. Filtering reads with insert sequence but without anchor information
  not_insert_only <- sum(grepl(insert_seqname, seqnames(gr))) != lengths(gr)

  # 2. Filtering related read cluster group that has not insert information
  insert_gr <- gr[grepl(insert_seqname, seqnames(gr))]
  has_insert <- any(start(insert_gr) < 6000)
  tmp_grp <- rep(x$group, lengths(x$read_annotation))
  group <- factor(tmp_grp, levels = unique(tmp_grp))
  group_has_insert <- any(LogicalList(split(has_insert, group)))
  usable_grp <- as.integer(levels(group))[group_has_insert]
  is_from_usable_grp <- tmp_grp%in%usable_grp

  # 3. Extracting regions that are close together (does not indicate transduction event)
  anchor_gr <- gr[!grepl(insert_seqname, seqnames(gr))]
  gr_range <- range(start(anchor_gr))
  gr_range[,1][abs(gr_range[,1]) > 999999999] <- 0
  gr_range[,2][abs(gr_range[,2]) > 999999999] <- 0
  range_ok <- (gr_range[,2] - gr_range[,1]) > 50000
  chr_ok <- lengths(unique(seqnames(anchor_gr))) > 1

  # 4. Getting regions has more than just genomic sequences
  tmp_anno <- anno[!grepl(":", anno)]
  has_seq <- sum(nchar(tmp_anno)) > 3
  orphan_ok <- lengths(anchor_gr) >= 1 & (has_seq | lengths(insert_gr) > 0) # Only contains genomic regions but with sequences (could possibly poly-A)
  region_ok <- lengths(anchor_gr) > 1 # Containing more than one 1 genomic regions

  # Combining all filters
  # is_informative <- not_insert_only & (orphan_ok | region_ok | range_ok | chr_ok | good_insert) & is_from_usable_grp
  is_informative <- not_insert_only & (orphan_ok | (region_ok & range_ok) | chr_ok) & is_from_usable_grp
  informative_anno <- split(flat_anno[is_informative], anno_grp[is_informative])
  elementMetadata(x)$informative_anno <- informative_anno
  elementMetadata(x)$is_usable <- lengths(x$informative_anno) > 0
  x
}
