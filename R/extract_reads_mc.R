#' @export
#' @importFrom vroom vroom
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom data.table rbindlist data.table

extract_info_reads <- function(bam, sorted_sam = NULL,
                               readin = 2.5E6, tmp_dir, threads = 8L,
                               out_dir, chromosome = c(1:22, "X", "Y", "KJ173426"),
                               samtools = "samtools",
                               interested_region = "/home/ctlaw/dicky/reference/RepeatMasker/L1PA1_2_3.txt",
                               cleanup = FALSE){
  if(!is.null(interested_region)){
    interested_region <- fread(interested_region, header = FALSE)
    interested_region <- interested_region[, 1:3]
    colnames(interested_region) <- c("RNAME", "START", "END")
  }
  if (is.null(tmp_dir)) {
    tempdir()
    tmp_dir <- tempfile()
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE))
  }

  if(is.null(sorted_sam)){
    tag <- gsub(".bam$", "", basename(bam))
    sorted_sam <- file.path(tmp_dir, paste0(tag, ".sam"))
    print(paste(Sys.time(), "Sorting reads"))
    cmd <- paste("(", samtools, "sort -@", threads,"-n -m 1500M -o -", bam, "|", samtools, "view -h -", ")>", sorted_sam, collapse = " ")
    system(cmd)
  }else{
    tag <- gsub(".sam$", "", basename(sorted_sam))
  }
  print(paste(Sys.time(), "Counting the number of total alignment records"))
  n_row <- system(paste0("wc -l ", sorted_sam), intern = TRUE)
  n_row <- as.integer(gsub(" .*", "", n_row))
  nheader <- get_header(sorted_sam)$nheader
  n_row <- n_row - nheader
  n_skip <- seq(0, n_row, by = readin)

  print(paste(Sys.time(), "Importing reads and extracting informative reads"))
  info <- bplapply(n_skip, function(x){
    reads_txt <- read_sam(sorted_sam, start = x + 1L, nrow = readin, select = c(1:6, 10))
    info <- get_info(reads_txt, chromosome = chromosome, interested_region = interested_region)
  }, BPPARAM =  MulticoreParam(workers = threads))

  # Extracting information from head and tail reads
  extra <- do.call(rbind, lapply(info, "[[", "head_tail_reads"))
  extra_info <- .extract_informative_reads(extra, chromosome = chromosome, interested_region = interested_region)

  print(paste(Sys.time(), "Exporting result"))
  # Exporting discordant reads
  discordant <- rbindlist(lapply(info, "[[", "discordant"))
  discordant <- rbind(discordant, extra_info$discordant)
  disc_out <- file.path(out_dir, paste0(tag, "_disc.txt"))
  write.table(discordant, disc_out, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

  # Exporting reads mapped to multiple regions
  interested <- rbindlist(lapply(info, "[[", "interested"))
  interested <- rbind(interested, extra_info$interested)
  interested_out <- file.path(out_dir, paste0(tag, "_interested.txt"))
  write.table(interested, interested_out, quote = FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

  # Cleaning up
  if(cleanup){
    system(paste("rm", sorted_sam, collapse = " "))
  }
}


#' @export
get_info <- function(input, chromosome = NULL, interested_region){
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
  system.time(informative_reads <- .extract_informative_reads(input, chromosome = chromosome, interested_region = interested_region))
  c(head_tail_reads=list(rbind(head_reads, tail_reads)), informative_reads)
}

#' @export
#' @importFrom data.table rbindlist data.table like
#' @importFrom Rsamtools bamFlagAsBitMatrix
.extract_informative_reads <- function(seq_info, chromosome = c(1:22, "X", "Y"), interested_region){
  # Converting sam flag to readable matrix
  selected_flag <- c("hasUnmappedMate", "isSupplementaryAlignment", "isProperPair")
  label <- bamFlagAsBitMatrix(seq_info$"FLAG", bitnames = selected_flag)
  colnames(label) <- c("HAS_UNMAPPED_MATE", "IS_SUPP_ALN", "IS_PROPER_PAIR")
  dt <- data.table(QNAME = seq_info$"QNAME", RNAME = seq_info$"RNAME",
                   START = seq_info$"POS", END = seq_info$"POS", MAPQ = seq_info$"MAPQ")
  dt <- data.table(dt, label)

  # Gathering information for statisitc analysis
  if(!is.null(chromosome)){
    NON_TARGET_RNAME <- !dt$"RNAME" %in% chromosome
  }else{
    NON_TARGET_RNAME <- FALSE
  }

  label <- data.table(QNAME = dt$QNAME, COUNT = 1,
                      HIGH_MAPQ = dt$"MAPQ" > 5,
                      LOW_MAPQ = dt$"MAPQ" <= 5,
                      NON_TARGET_RNAME = NON_TARGET_RNAME,
                      label)
  label$TARGET_IS_HIGH <- (label$HIGH_MAPQ & !label$"NON_TARGET_RNAME")

  # Gathering flag statistics in a paired reads way
  label_count <- label[, list(COUNT = sum(COUNT),
                              HIGH_MAPQ = sum(HIGH_MAPQ),
                              LOW_MAPQ = sum(LOW_MAPQ),
                              IS_PROPER_PAIR = sum(IS_PROPER_PAIR),
                              IS_SUPP_ALN = sum(IS_SUPP_ALN),
                              NON_TARGET_RNAME = sum(NON_TARGET_RNAME),
                              HAS_UNMAPPED_MATE = sum(HAS_UNMAPPED_MATE),
                              TARGET_IS_HIGH = sum(TARGET_IS_HIGH)),
                       by = "QNAME"]

  # Getting not interested reads
  all_low <- label_count$"COUNT" == label_count$"LOW_MAPQ" # reads with a pair of low mapq
  all_proper <- label_count$"COUNT" == label_count$"HIGH_MAPQ" & # sensible reads (not split reads, no unmapped reads and all high mapq)
    label_count$"IS_PROPER_PAIR" == label_count$"COUNT" &
    label_count$"IS_SUPP_ALN" == 0 &
    label_count$"HAS_UNMAPPED_MATE" == 0
  both_non_target <- label_count$"NON_TARGET_RNAME" == label_count$"COUNT" # both reads are on the non target chromosom
  unwanted <- all_proper | all_low | both_non_target

  # Getting reads with at least one with high mapq and one with low mapq
  is_high_low_mapq <- label_count$"TARGET_IS_HIGH" > 0 & (!unwanted) &
    (label_count$"LOW_MAPQ" != 0 | label_count$"HAS_UNMAPPED_MATE" != 0)
  high_low_mapq <- label_count$"QNAME"[is_high_low_mapq]

  # Getting reads with long Clipped sequence
  has_sh <- like(seq_info$"CIGAR", "S|H")
  sh_data <- seq_info[has_sh, ]
  is_long_left_clip <- .get_clip_length(sh_data$"CIGAR") > 10
  is_long_long_clip <- .get_clip_length(sh_data$"CIGAR", start = FALSE) > 10
  sh_q <- sh_data$"QNAME"[is_long_left_clip | is_long_long_clip]

  # Getting reads in interested regions
  if(!is.null(interested_region)){
    setkey(interested_region, "RNAME", "START", "END")
    region <- dt$"QNAME"[unique(foverlaps(dt, interested_region, type = "any", nomatch=NULL, which = TRUE)[[1]])]
  }else{
    region <- NULL
  }

  # Combining them into the interested RNAME
  interested <- c(high_low_mapq, sh_q)
  interested <- interested[!interested %in% label_count$"QNAME"[unwanted]]
  interested <- seq_info[seq_info$"QNAME" %in% c(interested, region), ]

  # Getting discor reads
  is_high_high_discor <- label_count$"COUNT" == label_count$"TARGET_IS_HIGH" &
    (label_count$"IS_PROPER_PAIR" != label_count$"COUNT" |
       label_count$"IS_SUPP_ALN" > 0) & !unwanted
  high_high_discor <- label_count$"QNAME"[is_high_high_discor]
  high_high_discor <- seq_info[seq_info$"QNAME" %in% high_high_discor,]

  list(interested = interested, discordant = high_high_discor)
}
