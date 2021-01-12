#' @export
#' @importFrom vroom vroom
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom data.table rbindlist data.table

extract_info_reads <- function(bam, readin = 2.5E6, tmp_dir, threads = 8L,
                               out_dir, chromosome = c(1:22, "X", "Y", "KJ173426")){
  # bam <- "/home/ctlaw/synology/HKU_88/bwa_output/HKU_liver_43/HKU_liver_tumor_43/ERR093456.bam"
  # tmp_dir <- "/home/ctlaw/dicky/analysis/Enhancer_hijack/temp"
  tag <- gsub(".bam$", "", basename(bam))

  if (is.null(tmp_dir)) {
    tempdir()
    tmp_dir <- tempfile()
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE))
  }

  sorted_sam <- file.path(tmp_dir, paste0(tag, ".sam"))

  print(paste(Sys.time(), "Sorting reads"))
  cmd <- paste("(samtools sort -@", threads,"-n -m 1500M -o -", bam, "|samtools view -h -", ")>", sorted_sam, collapse = " ")
  system(cmd)
  # txt <- "/home/ctlaw/dicky/analysis/Enhancer_hijack/temp/ERR093456.sam"
  n_row <- system(paste0("wc -l ", sorted_sam), intern = TRUE)
  n_row <- as.integer(gsub(" .*", "", n_row))
  n_skip <- seq(0, n_row, by = readin)

  print(paste(Sys.time(), "Importing reads and extracting informative reads"))
  info <- bplapply(n_skip, function(x){
    reads_txt <- read_sam(sorted_sam, start = x + 1L, nrow = readin, select = c(1:6, 10))
    info <- get_info(reads_txt, chromosome = chromosome)
  }, BPPARAM =  MulticoreParam(workers = threads))

  # Extracting information from head and tail reads
  extra <- do.call(rbind, lapply(info, "[[", "head_tail_reads"))
  extra_info <- .extract_informative_reads(extra, chromosome = chromosome)

  # Exporting discordant reads
  discordant <- rbindlist(lapply(info, "[[", "discordant"))
  discordant <- rbind(discordant, extra_info$discordant)
  # discordant$strand <- ifelse(bitwAnd(discordant$FLAG, 0x10), "-", "+")
  # pos_diff <- discordant[, list(range = range(POS)), by = QNAME]
  # pos_diff <- pos_diff[, list(diff = diff(range)), by = QNAME]
  # pos_diff_q <- pos_diff$QNAME[pos_diff$diff > 5000]
  # rname_strand <- discordant[, list(n_rname = length(unique(RNAME)),
  #                                   n_strand = length(unique(strand))), by = QNAME]
  # rs_diff_q <- rname_strand$QNAME[rname_strand$n_rname != 1]#|rname_strand$n_strand != 1]
  # discordant[discordant$QNAME %in% c(rs_diff_q, pos_diff_q)]
  disc_out <- file.path(out_dir, paste0(tag, "_disc.txt"))
  write.table(discordant, disc_out, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

  # Exporting reads mapped to multiple regions
  multi_mapped <- rbindlist(lapply(info, "[[", "multi_mapped"))
  multi_mapped <- rbind(multi_mapped, extra_info$multi_mapped)
  mm_out <- file.path(out_dir, paste0(tag, "_mm.txt"))
  write.table(multi_mapped, mm_out, quote = FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

  # Cleaning up
  system(paste("rm", sorted_sam, collapse = " "))
}


#' @export
get_info <- function(input, chromosome = NULL){
  readin <- N <- nrow(input)
  # Getting tail reads
  last_qname <- input[readin, "QNAME"]
  while(last_qname == input[N - 1, "QNAME"]){
    N <- N - 1
  }
  tail_reads <- input[N:readin,]

  # Getting head reads
  first_qname <- input[1, "QNAME"]
  Q <- 1
  while(first_qname == input[Q + 1, "QNAME"]){
    Q <- Q + 1
  }
  head_reads <- input[1:Q,]
  input <- input[-c(1:Q, N:readin),]
  system.time(informative_reads <- .extract_informative_reads(input, chromosome = chromosome))
  c(head_tail_reads=list(rbind(head_reads, tail_reads)), informative_reads)
}

#' @importFrom data.table rbindlist data.table
#' @importFrom Rsamtools bamFlagAsBitMatrix
.extract_informative_reads <- function(seq_info, chromosome = NULL){
  # Converting sam flag to readable matrix
  mat <- bamFlagAsBitMatrix(seq_info[["FLAG"]])

  # Gathering information for statisitc analysis
  if(!is.null(chromosome)){
    non_target_RNAME <- !seq_info[["RNAME"]] %in% chromosome
  }

  mat <- data.table(QNAME = seq_info[["QNAME"]], count = 1,
                    highMAPQ = seq_info[["MAPQ"]] > 5,
                    lowMAPQ = seq_info[["MAPQ"]] <= 5,
                    non_target_RNAME = non_target_RNAME,
                    mat)
  mat$targetIsHigh <- (mat$highMAPQ & !mat$non_target_RNAME)

  # Gathering flag statistics in a paired reads way
  paired_flag <- mat[, list(count = sum(count), highMAPQ = sum(highMAPQ),
                            lowMAPQ = sum(lowMAPQ), isProperPair = sum(isProperPair),
                            isSupplementaryAlignment = sum(isSupplementaryAlignment),
                            haveNonTargetRname = sum(non_target_RNAME),
                            isUnmappedQuery = sum(isUnmappedQuery),
                            targetIsHigh = sum(targetIsHigh)), by = "QNAME"]

  # Getting not interested reads
  all_low <- paired_flag$"count" == paired_flag$"lowMAPQ" # reads with a pair of low mapq
  all_proper <- paired_flag$"count" == paired_flag$"highMAPQ" & # sensible reads (not split reads, no unmapped reads and all high mapq)
    paired_flag$"isSupplementaryAlignment" == 0 &
    paired_flag$"isProperPair" == paired_flag$count &
    paired_flag$"isUnmappedQuery" == 0
  both_non_target <- paired_flag[["haveNonTargetRname"]] == paired_flag$count # both reads are on the non target chromosom
  unwanted <- all_proper | all_low | both_non_target

  # Getting reads with at least one with high mapq and one with low mapq
  is_high_low_mapq <- paired_flag[["highMAPQ"]] != 0 &
    (paired_flag[["lowMAPQ"]] != 0 | paired_flag[["isUnmappedQuery"]] != 0) &
    paired_flag[["targetIsHigh"]] > 0 & (!unwanted)
  high_low_mapq <- paired_flag[["QNAME"]][is_high_low_mapq]

  # Getting reads with long Clipped sequence
  sh_loc <- grep("S|H",seq_info$CIGAR, perl = TRUE)
  sh_data <- seq_info[sh_loc,]
  is_left_long_clip <- .get_clip_length(sh_data$CIGAR) > 10
  is_right_long_clip <- .get_clip_length(sh_data$CIGAR, start = FALSE) > 10
  sh_q <- sh_data$QNAME[is_left_long_clip|is_right_long_clip]
  interested <- c(high_low_mapq, sh_q)
  interested <- interested[!interested %in% paired_flag[["QNAME"]][unwanted]]
  high_low_mapq <- seq_info[seq_info[["QNAME"]] %in% interested,]

  # Getting discor reads
  is_high_high_discor <- paired_flag[["count"]] == paired_flag[["highMAPQ"]] &
    (paired_flag[["isProperPair"]] != paired_flag[["count"]] |
       paired_flag[["isSupplementaryAlignment"]] > 0) &
    !unwanted
  high_high_discor <- paired_flag[["QNAME"]][is_high_high_discor]
  high_high_discor <- seq_info[seq_info[["QNAME"]] %in% high_high_discor,]
  list(multi_mapped = high_low_mapq, discordant = high_high_discor)
}

