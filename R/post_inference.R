#' @export
#' @importFrom data.table rbindlist
get_somatic <- function(t_df, nt_total, excluded_range = 5000){
  conf_t <- get_confident(t_df)
  nt_total <- lapply(nt_total, unlist)
  nt_anno <- unlist(lapply(, as.character))
  nt_anno <- nt_anno[!grepl("Hot_L1_polyA:", nt_anno)]
  nt_anno <- convert_character2gr(nt_anno)
  somatic <- subsetByOverlaps(conf_t, nt_anno + excluded_range, invert = TRUE, ignore.strand = TRUE)
  somatic
}

#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
convert_dt_to_gr <- function(df){
  non_na_3p <- !is.na(df[["3p_insert_start"]])
  non_na_5p <- !is.na(df[["5p_insert_end"]])
  non_na_both <- non_na_3p & non_na_5p
  have_both_end <- GRanges(df[non_na_both]$`5p_chr`, IRanges(df[non_na_both]$`5p_genomic_break`,
                                                             df[non_na_both]$`5p_genomic_break`),
                           strand = df[non_na_both]$`3p_insert_Orientation`)
  elementMetadata(have_both_end) <- df[non_na_both]
  only_5p <- df[non_na_5p&!non_na_both]
  only_5p_end <- GRanges(only_5p$`5p_chr`,
                         IRanges(only_5p$`5p_genomic_break`, only_5p$`5p_genomic_break`),
                         strand = only_5p$`5p_insert_Orientation`)
  elementMetadata(only_5p_end) <- only_5p
  only_3p <- df[non_na_3p&!non_na_both]
  only_3p_end <- GRanges(only_3p$`3p_chr`,
                         IRanges(only_3p$`3p_genomic_break`, only_3p$`3p_genomic_break`),
                         strand = only_5p$`3p_insert_Orientation`)
  elementMetadata(only_3p_end) <- only_3p

  combined <- list(both = have_both_end, only_5p = only_5p_end, only_3p = only_3p_end)
  combined
}

get_exact <- function(df){
  exact_3p <- df$"3p_is_exact"
  exact_5p <- df$"5p_is_exact"
  exact_3p[is.na(exact_3p)] <- FALSE
  exact_5p[is.na(exact_5p)] <- FALSE
  selected <- df[exact_3p & exact_5p, ]
  selected
}

#' @export
#' @importFrom Biostrings getSeq
get_motif <- function(df){
  Hsapiens <- BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
  selected <- get_exact(df)
  start_at_5p <- selected$"5p_insert_start" < selected$"3p_insert_start"
  end_at_5p <- !start_at_5p
  cur_strand <- selected[end_at_5p & selected$"3p_insert_orientation"=="-", ]
  if(length(cur_strand)>0){
    cur_seq_1 <- getSeq(Hsapiens, cur_strand$"3p_chr",
                        start = cur_strand$"3p_genomic_break" - 10,
                        end = cur_strand$"3p_genomic_break" + 10)
  }else{
    cur_seq_1 <- NULL
  }

  cur_strand <- selected[(!end_at_5p) & selected$"5p_insert_orientation"=="+", ]
  if(length(cur_strand)>0){
    cur_seq_2 <- getSeq(Hsapiens, cur_strand$"5p_chr",
                        start = cur_strand$"5p_genomic_break" - 10,
                        end = cur_strand$"5p_genomic_break" + 10)
  }else{
    cur_seq_2 <- NULL
  }
  list(cur_seq_1, cur_seq_2)
}

get_full <- function(df){
  non_na <- df[!is.na(df$has_polyA),]
  non_na <- non_na[!is.na(non_na[["3p_insert_start"]])&!is.na(non_na[["5p_insert_end"]]),]
  full <- non_na[abs(non_na$`3p_insert_start`-non_na$`5p_insert_end`) > 5500, ]
  out <- GRanges(full$`5p_chr`, IRanges(full$`5p_genomic_regions_start`, full$`5p_genomic_regions_end`), strand=full$`3p_insert_Orientation`)
  elementMetadata(out) <- full
  out
}
