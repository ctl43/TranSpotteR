#' @export
#' @importFrom vroom vroom
#' @importFrom BiocParallel bplapply MulticoreParam SnowParam SerialParam
#' @importFrom data.table rbindlist data.table

extract_info_reads <- function(bam, sorted_sam = NULL,
                               readin = 2.5E6, tmp_dir, BPPARAM = SerialParam(),
                               out_dir, chromosome = c(1:22, "X", "Y", "KJ173426"),
                               samtools = "samtools",
                               interested_region = "/home/ctlaw/dicky/reference/RepeatMasker/L1PA1_2_3.txt",
                               interested_seq = "~/dicky/reference/fasta/line1_reference/hot_L1_polyA.fa",
                               cleanup = FALSE){

  if(class(BPPARAM)=="SerialParam"){
    threads <- 1
  }else{
    threads <- BPPARAM$workers
  }

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
    # cmd <- paste("(", samtools, "sort -@", threads,"-n -m 1500M -o -", bam, "|", samtools, "view -h -", ")>", sorted_sam, collapse = " ")
    cmd <- paste("(", samtools, "sort -@", threads,"-n -m 1500M -o -", bam, "|", samtools, "view -|cut -f 1,2,3,4,5,6,10", ")>", sorted_sam, collapse = " ")
    system(cmd)
  }else{
    tag <- gsub("\\..*", "", basename(sorted_sam))
  }
  print(paste(Sys.time(), "Counting the number of total alignment records"))
  n_row <- system(paste0("wc -l ", sorted_sam), intern = TRUE)
  n_row <- as.integer(gsub(" .*", "", n_row))
  # nheader <- get_header(sorted_sam)$nheader
  # n_row <- n_row - nheader
  n_skip <- seq(0, n_row, by = readin)

  print(paste(Sys.time(), "Importing reads and extracting informative reads"))
  # tmp_files <- replicate(length(n_skip), tempfile())
  # info <- bpmapply(function(x, sam, tmp_file, readin, chromosome, interested_region, get_info){
  #   reads_txt <- TranSpotteR::read_sam(sam, start = x + 1L, nrow = readin, select = c(1:6, 10))
  #   TranSpotteR::get_info(reads_txt, chromosome = chromosome, interested_region = interested_region, tmp_file = tmp_file)
  # }, BPPARAM =  BPPARAM, x = n_skip, sam = sorted_sam, tmp_file = tmp_files, readin = readin,
  # chromosome = list(chromosome), interested_region = list(interested_region), SIMPLIFY = FALSE)


  info <- bpmapply(function(x, sam, tmp_file, readin, chromosome, interested_region, get_info){
    what <- list(X1 = "c", X2 = "i", X3 = "c", X4 = "i", X5 = "i", X6 = "c", X10 = "c")
    reads_txt <- vroom(sam, delim = "\t", skip = x, col_names = FALSE, n_max = readin, col_types = what)
    # what <- c("character","integer","character","integer","integer","character","character")
    # reads_txt <- data.table::fread(sam, skip = x, nrows = readin, colClasses = what)
    colnames(reads_txt) <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "SEQUENCE")
    TranSpotteR::get_info(reads_txt, chromosome = chromosome, interested_region = interested_region, tmp_file = tmp_file)
  }, BPPARAM =  BPPARAM, x = n_skip, sam = sorted_sam, readin = readin,
  chromosome = list(chromosome), interested_region = list(interested_region), SIMPLIFY = FALSE)

  # Extracting information from head and tail reads
  extra <- do.call(rbind, lapply(info, "[[", "head_tail_reads"))
  extra_info <- .extract_informative_reads(extra, chromosome = chromosome, interested_region = interested_region)

  print(paste(Sys.time(), "Exporting result"))
  # Exporting discordant reads
  discordant <- rbindlist(lapply(info, "[[", "discordant"))
  discordant <- rbind(discordant, extra_info$discordant)
  has_multiple_seqnames <- lengths(unique(CharacterList(split(discordant$RNAME, discordant$QNAME)))) > 1
  start_range <- range(IntegerList(split(discordant$POS, discordant$QNAME)))
  keep <- rownames(start_range)[(start_range[,2] - start_range[,1] > 5000) | has_multiple_seqnames] # Extracting read that are far apart enough
  discordant <- discordant[discordant$QNAME %in% keep,]
  # to_keep <- discordant$QNAME[discordant$MAPQ > disc_min_mapq]
  # discordant <- discordant[discordant$QNAME %in% to_keep, ]
  disc_out <- file.path(out_dir, paste0(tag, "_disc.txt"))
  write.table(discordant, disc_out, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

  # Exporting reads mapped to multiple regions
  interested <- rbindlist(lapply(info, "[[", "interested"))
  interested <- rbind(interested, extra_info$interested)
  is_first <- !!bitwAnd(interested$FLAG, 0x40)
  interested$QNAME_id <- paste0(interested$QNAME,"_",as.integer(!is_first) + 1)
  non_hard_clipped <- interested[!like(interested$"CIGAR", "H")]

  # Aligning the read to interested sequence
  reads <- non_hard_clipped$SEQUENCE
  names(reads) <- non_hard_clipped$QNAME_id
  aln <- bwa_alignment(reads, ref = interested_seq, samtools_param = "-F 2308")
  has_long_clip <- .get_clip_length(aln$CIGAR) > 10|.get_clip_length(aln$CIGAR, start = FALSE) > 10
  q <- sub("_.*", "", aln$QNAME)
  both_here <- lengths(CharacterList(split(aln$QNAME, q)))>1
  both_here_q <- unique(q)[both_here]
  good_q <- q[!(q%in%both_here_q)|has_long_clip]
  interested <- interested[interested$QNAME %in% good_q]

  interested_out <- file.path(out_dir, paste0(tag, "_interested.txt"))
  write.table(interested, interested_out, quote = FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

  # Cleaning up
  if(cleanup){
    system(paste("rm", sorted_sam, collapse = " "))
  }
}


#' @export
get_info <- function(input, chromosome = NULL, interested_region, tmp_file){
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
  informative_reads <- .extract_informative_reads(input, chromosome = chromosome, interested_region = interested_region)#, tmp_file = tmp_file)
  # .extract_informative_reads(input, chromosome = chromosome, interested_region = interested_region, tmp_file = tmp_file)
  c(head_tail_reads = list(rbind(head_reads, tail_reads)), informative_reads)
  # rbind(head_reads, tail_reads)
}

#' @export
#' @importFrom data.table rbindlist data.table like setkey foverlaps
#' @importFrom Rsamtools bamFlagAsBitMatrix
.extract_informative_reads <- function(seq_info, chromosome = c(1:22, "X", "Y", "KJ173426"), interested_region){#, tmp_file){
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
    ol <- foverlaps(dt, interested_region, type = "any", nomatch=NULL, which = TRUE)
    region <- dt$"QNAME"[unique(ol[[1]])]
    region <- region[!region%in%(label_count$QNAME[all_low])]
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
  # tmp_disc <- paste0(tmp_file, "_disc.txt")
  # tmp_interested <- paste0(tmp_file, "_interested.txt")
  # write.table(interested, tmp_interested, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
  # write.table(high_high_discor, tmp_disc, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
  list(interested = interested, discordant = high_high_discor)
}
