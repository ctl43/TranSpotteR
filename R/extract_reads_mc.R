#' @export
#' @importFrom vroom vroom
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom data.table rbindlist data.table

extract_info_reads <- function(bam, readin = 2.5E6, tmp_dir, mc.cores = 12L, out_dir){
  # bam <- "/home/ctlaw/synology/HKU_88/bwa_output/HKU_liver_43/HKU_liver_tumor_43/ERR093456.bam"
  # tmp_dir <- "/home/ctlaw/dicky/analysis/Enhancer_hijack/temp"
  tag <- gsub(".bam", "", basename(bam))
  txt <- file.path(tmp_dir, paste0(tag, ".txt"))
  print(paste(Sys.time(), "Sorting reads"))
  cmd <- paste("(samtools sort -@", mc.cores,"-n -m 1500M -o /dev/stdout", bam, "|samtools view -h /dev/stdin", ")>", txt, collapse = " ")
  system(cmd)
  n_row <- system(paste0("wc -l ", txt), intern = TRUE)
  n_row <- as.integer(gsub(" .*", "", n_row))
  n_skip <- seq(0, n_row, by = readin)

  print(paste(Sys.time(), "Importing reads and extracting informative reads"))
  system.time(info <- bplapply(n_skip[1], function(x){
    reads_txt <- read_sam(txt, start = x + 1L, nrow = 2.5e6, select = c(1:6, 10))
    info <- get_info(reads_txt)
  }, BPPARAM =  MulticoreParam(workers = 12L)))

  # Extracting information from head and tail reads
  extra <- do.call(rbind, lapply(info, "[[", "head_tail_reads"))
  extra_info <- .extract_informative_reads(extra)

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
  write.table(discordant, disc_out, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

  # Exporting reads mapped to multiple regions
  multi_mapped <- rbindlist(lapply(info, "[[", "multi_mapped"))
  multi_mapped <- rbind(multi_mapped, extra_info$multi_mapped)
  mm_out <- file.path(out_dir, paste0(tag, "_mm.txt"))
  write.table(multi_mapped, mm_out, quote = FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

  print(paste(Sys.time(), "Mapping multi mapped reads to line1 and polyA"))
  # # Identifying reads mapped to line1
  # non_supp_second <- multi_mapped[ !bitwAnd(multi_mapped$FLAG, 0x100) & !bitwAnd(multi_mapped$FLAG, 0x800),]
  # multi_mapped_seq <- DNAStringSet(non_supp_second$SEQUENCE)
  # is_first <- !bitwAnd(non_supp_second$FLAG, 0x40)
  # names(multi_mapped_seq) <- paste0(non_supp_second$QNAME)
  # mm_fa_out <- file.path(out_dir, paste0(tag, "_multi_mapped.fasta"))
  # writeFasta(multi_mapped_seq, mm_fa_out)
  # cmd <- paste("bwa mem -t 5", "/home/ctlaw/dicky/reference/fasta/line1_reference/L1complete.fa", mm_fa_out,
  #        "|samtools view -F 4 -|cut -f 1,2,3,4,5,6,10")
  # mm_aln <- read.delim(text=system(cmd, intern = TRUE),stringsAsFactors = FALSE, header=FALSE)
  # line1_mm <- multi_mapped[multi_mapped$QNAME%in%mm_aln$V1,]
  # line1_mm_out <- file.path(out_dir, paste0(tag, "_line1_mm.txt"))
  # write.table(line1_mm, line1_mm_out, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
  #
  # # Exporting line1 anchor (coverting it to bed format)
  # line1_anchor <- line1_mm[!bitwAnd(line1_mm$FLAG, 0x100)&!bitwAnd(line1_mm$FLAG, 0x800)&line1_mm$MAPQ>=5,]
  # line1_anchor_out <- file.path(out_dir, paste0(tag, "_line1_anchor.txt"))
  # write.table(line1_anchor, line1_anchor_out, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")


  # Cleaning up
  system(paste("rm",txt, collapse = " "))
}


#' @export
get_info <- function(input){
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
  informative_reads <- .extract_informative_reads(input)
  c(head_tail_reads=list(rbind(head_reads, tail_reads)), informative_reads)
}

#' @importFrom data.table rbindlist data.table
#' @importFrom Rsamtools bamFlagAsBitMatrix
.extract_informative_reads <- function(seq_info){
  # Converting sam flag to readable matrix
  mat <- bamFlagAsBitMatrix(seq_info[["FLAG"]])

  # Gathering information for statisitc analysis
  mat <- data.table(QNAME = seq_info[["QNAME"]], count = 1,
                    highMAPQ = seq_info[["MAPQ"]] > 5,
                    lowMAPQ = seq_info[["MAPQ"]] <= 5,
                    # non_target_RNAME =! seq_info[["RNAME"]] %in% c(1:22, "X", "Y", "KJ173426"),
                    mat)

  # Gathering flag statistics in a paired reads way
  paired_flag <- mat[, list(count = sum(count), highMAPQ = sum(highMAPQ),
                            lowMAPQ = sum(lowMAPQ), isProperPair = sum(isProperPair),
                            isSupplementaryAlignment = sum(isSupplementaryAlignment),
                            # haveNonTargetRname = sum(non_target_RNAME),
                            isUnmappedQuery = sum(isUnmappedQuery)), by = "QNAME"]

  # Getting not interested reads
  all_low <- paired_flag$"count" == paired_flag$"lowMAPQ" # reads with a pair of low mapq
  all_proper <- paired_flag$"count" == paired_flag$"highMAPQ" & # sensible reads (not split reads, no unmapped reads and all high mapq)
    paired_flag$"isSupplementaryAlignment" == 0 &
    paired_flag$"isProperPair" == paired_flag$count &
    paired_flag$"isUnmappedQuery" == 0
  # have_non_target <- paired_flag[["haveNonTargetRname"]] > 0
  unwanted <- all_proper | all_low | have_non_target

  # Getting reads with at least one with high mapq and one with low mapq
  is_high_low_mapq <- paired_flag[["highMAPQ"]] != 0 &
    (paired_flag[["lowMAPQ"]] != 0 | paired_flag[["isUnmappedQuery"]] != 0) &
    !unwanted
  high_low_mapq <- paired_flag[["QNAME"]][is_high_low_mapq]
  high_low_mapq <- seq_info[seq_info[["QNAME"]] %in% high_low_mapq,]

  # Getting discor reads
  is_high_high_discor <- paired_flag[["count"]] == paired_flag[["highMAPQ"]] &
    (paired_flag[["isProperPair"]] != paired_flag[["count"]] | paired_flag[["isSupplementaryAlignment"]] > 0) &
    !unwanted
  high_high_discor <- paired_flag[["QNAME"]][is_high_high_discor]
  high_high_discor <- seq_info[seq_info[["QNAME"]] %in% high_high_discor,]
  list(multi_mapped = high_low_mapq, discordant = high_high_discor)
}

