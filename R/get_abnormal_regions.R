# library(ShortRead)
# #Getting false-positive prone regions
# cluster_files <- dir("/home/ctlaw/dicky/analysis/Enhancer_hijack/temp/TranSpotteR_output/previous",
#                      pattern = "_TranSpotteR_clusters.rds", full.names = TRUE)
#
# question_regions <- mclapply(cluster_files, function(x){
#   filtered <- gsub("previous", "", gsub("clusters", "annotated", x))
#   y <- readRDS(x)
#   z <- readRDS(filtered)
#   q <- y[!names(y) %in% z$cluster_origin]
#   GRanges(seqnames(q), IRanges(start(q), end(q)), strand(q))
# }, mc.cores = 10L)
#
#
# question_regions <- GRangesList(question_regions)
# flat_question_regions <- unlist(question_regions)
# regions_coverage <- coverage(flat_question_regions)
# saveRDS(regions_coverage, "/home/ctlaw/dicky/tools/TranSpotteR/data/problematic_count.rds")
#
# abnormal <- IRangesList(lapply(regions_coverage, function(x)IRanges(which(x > 30))))
# abnormal <- GRanges(rep(names(regions_coverage), lengths(abnormal)), unlist(abnormal))
# saveRDS(abnormal, "/home/ctlaw/dicky/tools/TranSpotteR/data/false_positive_region.rds")
