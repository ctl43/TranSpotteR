#' @export
#' @importFrom Rsamtools bamFlagAsBitMatrix
extract_reads <- function(bam, tmp_dir = NULL,
                          out_dir, chromosome = c(1:22, "X", "Y", "KJ173426"),
                          samtools = "samtools",
                          interested_region = NULL,
                          interested_seq = NULL,
                          keep_intermediate = FALSE,
                          threads = 5){
# A wrapper to extract discordant reads by calling samtools and filtering those obvious non-informatic read pairs
  if (is.null(tmp_dir)) {
    tmp_dir <- tempdir()
    if(!dir.exists(tmp_dir)){
      dir.create(tmp_dir)
    }
  }
  tmp_file <- basename(tempfile("tmp_"))
  tag <- gsub("\\..*", "", basename(bam))

  # tmp_dir <- "/home/ctlaw/temp"
  region_sam <- paste0(tag, "_", tmp_file, "_region.sam")
  region_sam <- file.path(tmp_dir, region_sam)
  region_cmd <- paste(c("(",samtools, "view", "-L", interested_region, bam,"|cut -f 1,2,3,4,5,6,10)>", region_sam), collapse = " ")

  # Getting discordant read pairs
  disc_sam <- paste0(tag, "_", tmp_file, "_disc.sam")
  disc_sam <- file.path(tmp_dir, disc_sam)
  disc_cmd <- paste(c("(", samtools, "view", "-F 14", bam,"|cut -f 1,2,3,4,5,6,10)>", disc_sam), collapse = " ")

  print(paste(Sys.time(), "Extracting reads overlapping with the interested regions and discordant reads"))
  if(threads > 1){
    parallel_cmd <- paste(c(region_cmd, disc_cmd, "wait"), collapse = "&")
    system(parallel_cmd)
  }else{
    system(region_cmd)
    system(disc_cmd)
  }
  Sys.sleep(0.5)

  # Checking whether their partner is here and only getting those with partner
  interested <- read_aln_info(region_sam)
  interested$QNAME_id <- paste0(interested$QNAME, "_", bamFlagAsBitMatrix(interested$FLAG, bitnames = "isSecondMateRead")[, 1] + 1)
  interested_pair_check <- data.table(QNAME = interested$QNAME,
                                      QNAME_id = interested$QNAME_id)[!grepl("H", interested$CIGAR),]
  interested_pair_check <- unique(interested_pair_check)
  interested_pair_check <- interested_pair_check[, list(QNAME_id_COUNT = length(QNAME_id)), by = QNAME]
  interested_pair_read <- interested_pair_check$QNAME[interested_pair_check$QNAME_id_COUNT == 2]
  setkey(interested, QNAME)
  interested <- interested[J(interested_pair_read), ]

  # Getting meaningful discordant reads
  discordant <- read_aln_info(disc_sam)
  discordant <- discordant[!discordant$QNAME %in% interested$QNAME]
  discordant$QNAME_id <- paste0(discordant$QNAME, "_", bamFlagAsBitMatrix(discordant$FLAG, bitnames = "isSecondMateRead")[, 1] + 1)
  disc_check <- discordant[, list(RANGE = diff(range(POS)), N_RNAME = length(unique(RNAME))), by = QNAME]
  pass_disc_check <- disc_check[(disc_check$RANGE > 5000 | disc_check$N_RNAME > 1), "QNAME"]
  setkey(discordant, QNAME)
  discordant <- discordant[J(pass_disc_check)]

  # Combining the reads and only getting reads that has at least one read that is not located at potential insert/having long clipped end
  combined <- rbind(interested, discordant)
  no_hard_clipped <- combined[!grepl("H", combined$CIGAR), ]
  seq <- no_hard_clipped$SEQUENCE
  names(seq) <- no_hard_clipped$QNAME_id
  aln <- setDT(bwa_alignment(seq, ref = interested_seq, samtools_param = "-F 2308", threads = threads))
  aln$ORIGIN_READ <- gsub("_.*", "", aln$QNAME)
  left_clipped_len <- .get_clip_length(aln$CIGAR)
  right_clipped_len <- .get_clip_length(aln$CIGAR, start = FALSE)
  aln$HAS_LONG_CLIPPED <- left_clipped_len > 20 | right_clipped_len > 20
  insert_check <- aln[, list(HAS_LONG_CLIPPED = sum(HAS_LONG_CLIPPED), QNAME_id_COUNT = length(QNAME)), by = ORIGIN_READ]
  failed_insert_check <- insert_check[insert_check$QNAME_id_COUNT == 2 & insert_check$HAS_LONG_CLIPPED == 0][[1]]
  setkey(combined, QNAME)
  combined <- combined[-combined[failed_insert_check, which = TRUE]]

  # At least has one read with high mapq
  mapq_check <- data.table(QNAME = combined$QNAME, HIGH_MAPQ = (combined$RNAME %in% chromosome & combined$MAPQ > 5))
  mapq_check <- mapq_check[, list(HIGH_MAPQ = sum(HIGH_MAPQ)), by = QNAME]
  has_high_mapq <- mapq_check$QNAME[!mapq_check$HIGH_MAPQ == 0]
  combined <- combined[J(has_high_mapq), ]
  write.table(combined, file.path(out_dir, paste0(tag, "_extracted.txt")),
              sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Cleaning up
  if(!keep_intermediate){
    on.exit(unlink(c(region_sam, disc_sam)))
  }
}

#' @export
read_aln_info <- function(x){
  what <- c("character","integer","character","integer","integer","character","character")
  dt <- fread(x, colClasses = what)
  colnames(dt) <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "SEQUENCE")
  dt
}

