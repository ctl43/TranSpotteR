#' @export
#' @importFrom S4Vectors elementMetadata 'elementMetadata<-'
#' @importFrom BiocGenerics start end grepl
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges LogicalList

filter_clusters <- function(x, insert_seqname = "Hot_L1_polyA", additional_filter = c("has_polyA", "has_polyT"))
  # As the inference process is complicated and takes time,
  # this function is to filter those obviously useless information to reduce the computation time
  # Written by Cheuk-Ting Law
{

  # Getting the mapping information
  anno <- x$contig_detail
  anno <- CharacterList(lapply(anno, "[[", i = "annotation"))

  # Converting into a list of GRanges objects
  anno_ranges <- anno[S4Vectors::grepl(":", anno)]
  gr <- convert_character2gr(anno_ranges)
  gr <- unique(gr)

  # 1. Filtering reads with insert sequence but without anchor information
  is_insert <- grepl(insert_seqname, seqnames(gr))
  not_insert_only <- !all(is_insert)

  # 2. Filtering related read cluster group that has not insert information
  insert_gr <- gr[grepl(insert_seqname, seqnames(gr))]
  has_insert <- any(start(insert_gr) < 6010)
  group_has_insert <- any(LogicalList(split(has_insert, factor(x$group, levels = unique(x$group)))))
  usable_grp <- unique(x$group)[group_has_insert]
  is_from_usable_grp <- x$group %in% usable_grp

  # 3. Extracting regions that are close together (does not indicate transduction event)
  anchor_gr <- gr[!is_insert]
  gr_range <- range(start(anchor_gr))
  gr_range[,1][abs(gr_range[,1]) > 999999999] <- 0
  gr_range[,2][abs(gr_range[,2]) > 999999999] <- 0
  range_ok <- (gr_range[,2] - gr_range[,1]) > 50000
  chr_ok <- lengths(unique(seqnames(anchor_gr))) > 1

  # 4. Getting regions has more than just genomic sequences
  has_target <- lapply(additional_filter, function(x)any(LogicalList(lapply(x$contig_detail, "[[", i = additional_filter))))
  has_target <- Reduce("|", has_target)
  orphan_ok <- lengths(anchor_gr) >= 1 & (has_target | lengths(insert_gr) > 0) # Only contains genomic regions but with sequences (could possibly poly-A)
  region_ok <- lengths(anchor_gr) > 1 # Containing more than one 1 genomic regions

  # Combining all filters
  is_informative <- not_insert_only & (orphan_ok | (region_ok & range_ok) | chr_ok) & is_from_usable_grp
  x[is_informative]
}
