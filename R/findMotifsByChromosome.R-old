library(RPostgreSQL)
library(RUnit)
library(GenomicRanges)
require(BSgenome.Hsapiens.UCSC.hg38)

if(!exists("parseChromLocString"))
   source("~/github/trena/R/utils.R")
#------------------------------------------------------------------------------------------------------------------------
recallDataBaseTableStructure <- function()
{
   db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="lymphoblast", host="whovian")
   db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="brain_hint_20", host="khaleesi")
   dbListTables(db)
   tbl.hits <- dbGetQuery(db, "select * from hits limit 2")
   hits.coi <<- colnames(tbl.hits)
   tbl.regions <- dbGetQuery(db, "select * from regions limit 2")
   regions.coi <<- colnames(tbl.regions)

     # tbl.regions
     #                  loc chrom  start endpos
     # 1 chr9:111122-111132  chr9 111122 111132
     # 2 chr9:425141-425151  chr9 425141 425151

     # tbl.hits - transposed
     #
     # as.data.frame(t(dbGetQuery(db, "select * from hits limit 2")))
     #                                             1                                  2
     #  1) loc                        chr9:111122-111132                 chr9:111122-111132
     #  2) fp_start                               111111                             111120
     #  3) fp_end                                 111124                             111139
     #  4) type                       motif.in.footprint                 motif.in.footprint
     #  5) name       Hsapiens-jaspar2016-RUNX1-MA0002.1 Hsapiens-jaspar2016-RUNX1-MA0002.1
     #  6) length                                     11                                 11
     #  7) strand                                      -                                  -
     #  8) sample_id                         ENCSR000DBW                        ENCSR000DBW
     #  9) method                                   HINT                               HINT
     # 10) provenance                brain_hint_20.minid                brain_hint_20.minid
     # 11) score1                                     15                                 15
     # 12) score2                                12.4045                            12.4045
     # 13) score3                               3.48e-05                           3.48e-05
     # 14) score4                                   <NA>                               <NA>
     # 15) score5                                   <NA>                               <NA>
     # 16) score6                                   <NA>                               <NA>

} # recallDataBaseTableStructure
#------------------------------------------------------------------------------------------------------------------------
createFootprintTablesForRegion <- function(chromosome, start.loc, end.loc)
{
   if(!exists("hits.coi")) recallDataBaseTableStructure()
   print(hits.coi)
       #------------------------------------------------------------
       # first read the footprint calls
       #------------------------------------------------------------

    searchString <- sprintf('find ../whole_genome -name "*%s_fp.bed"', chromosome)
    files <- system(searchString, intern=TRUE)
    printf("file count: %d", length(files))
    fps <- list()
    for(f in files){
       f.path <- f # file.path(dir, f)
       printf("--- %s", f.path)
       tbl <- read.table(f.path, sep="\t", as.is=TRUE, nrow=-1)
       print(nrow(tbl))
       colnames(tbl) <- c("chrom", "start", "end", "name", "score", "strand", "ignore")
       print(head(tbl))
       tbl <- subset(tbl, chrom==chromosome & start >= start.loc & end <= end.loc)
       print(nrow(tbl))
       fps[[f]] <- tbl
       } # for f

       #------------------------------------------------------------
       # combine the footprint lists into a single data.frame
       #------------------------------------------------------------

    # browser()
    tbl.fps <- do.call(rbind, fps)[, c("chrom", "start", "end", "score")]
    rownames(tbl.fps) <- NULL
    dim(tbl.fps)
    fp.order <- order(tbl.fps$chrom, tbl.fps$start)
    tbl.fps <- tbl.fps[fp.order,]

       #------------------------------------------------------------
       # expand the footprints by <shoulder> in both directions
       #------------------------------------------------------------

    shoulder <- 10
    tbl.fpsExpanded <- tbl.fps
    tbl.fpsExpanded$start <- tbl.fpsExpanded$start - shoulder
    tbl.fpsExpanded$end <- tbl.fpsExpanded$end + shoulder

       #------------------------------------------------------------
       # collapse overlapping footprints
       #------------------------------------------------------------

    tbl.fpsExpanded <- as.data.frame(union(GRanges(tbl.fpsExpanded), GRanges(tbl.fpsExpanded)))
    colnames(tbl.fpsExpanded) <- c("chrom", "start", "end", "width", "strand")
    dim(tbl.fpsExpanded)
    filename <- sprintf("%s:%d-%d.footprints-combined.RData", chromosome, start.loc, end.loc)
    printf("saving %d combined expanded footprints to %s", nrow(tbl.fpsExpanded), filename)
    save(tbl.fpsExpanded, file=filename)
    write.table(tbl.fpsExpanded, sep="\t", quote=FALSE

} # createFootprintTablesForRegion
#------------------------------------------------------------------------------------------------------------------------
fimoPrep <- function(tbl.fpsExpanded, )
{
   sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg38, tbl.fpsExpanded$chrom, tbl.fpsExpanded$start, tbl.fpsExpanded$end))
   names(sequences) <- sprintf("%s:%d-%d", tbl.fpsExpanded$chrom, tbl.fpsExpanded$start, tbl.fpsExpanded$end)
   writeXStringSet(sequences, "tmp2.fa")


} # fimoPrep
#------------------------------------------------------------------------------------------------------------------------
# get fimo motif hits to these regions at the specified threshold
callMotifs <- function()
{
   library(FimoClient)

   FIMO_HOST <- "khaleesi"
   FIMO_HOST <- "localhost"
   FIMO_PORT <- 60000
   fc <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)
   tbl.motifs <- requestMatchForRegions(fc, tbl.fpsExpanded, "hg38", 1e-4)

   dim(tbl.motifs)
   colnames(tbl.motifs)[grep("start", colnames(tbl.motifs))] <- "motifStart"
   colnames(tbl.motifs)[grep("stop", colnames(tbl.motifs))] <- "motifEnd"

        #-----------------------------------------------------------------------------------------------
        # expand "sequence_name" to reproduce footprint start/stop locs, and then motif start/stop locs
        #-----------------------------------------------------------------------------------------------

   locs <- lapply(tbl.motifs$sequence_name, parseChromLocString)
   tbl.locs <- do.call(rbind, lapply(lapply(tbl.motifs$sequence_name, parseChromLocString), as.data.frame))
   tbl.motifs$chrom <- tbl.locs$chrom
   tbl.motifs$motifStart <- tbl.motifs$motifStart + tbl.locs$start
   tbl.motifs$motifEnd <- tbl.motifs$motifEnd + tbl.locs$start
   tbl.motifs$loc <- with(tbl.motifs, sprintf("%s:%d-%d", chrom, motifStart, motifEnd))

   colnames(tbl.locs)[grep("start", colnames(tbl.locs))] <- "fp_start"
   colnames(tbl.locs)[grep("end", colnames(tbl.locs))] <- "fp_end"

   tbl.motifs <- cbind(tbl.locs, tbl.motifs)

   fp.offset <- tbl.motifs$fp_start
        # our database scheme improperly uses fp_start and fp_end for motif start and end.
        # we accomdate that old choice here:
   tbl.motifs$fp_start <- tbl.motifs$motifStart + fp.offset
   tbl.motifs$fp_end  <- tbl.motifs$motifEnd + fp.offset
   # tbl.motifs$motifStart <- tbl.motifs$motifStart + tbl.motifs$fp_start
   # tbl.motifs$motifEnd <- tbl.motifs$motifEnd + tbl.motifs$fp_start

   colnames(tbl.motifs)     # "chrom", "fp_start", "fp_end", "motif", "sequence_name", "motifStart",
                            # "motifEnd", "strand", "score", "pValue", "qValue", "matched_sequence"

   tbl.newMotifs <- tbl.motifs[, c("loc", "chrom", "motifStart", "motifEnd", "motif", "strand", "score",
                                  "pValue", "qValue", "matched_sequence")]
   tbl.newMotifs$type <- "motif.in.footprint"
   tbl.newMotifs$method <- "HINT"
   tbl.newMotifs$provenance <- "placental_hint"
   tbl.newMotifs$length <- 1 + tbl.newMotifs$motifEnd - tbl.newMotifs$motifStart
   # tbl.newMotifs$loc <- sprintf("%s:%d-%d", tbl.newMotifs$chrom, tbl.newMotifs$motifStart, tbl.newMotifs$motifEnd)
   colnames(tbl.newMotifs)[grep("^motif$", colnames(tbl.newMotifs))] <- "name"
   colnames(tbl.newMotifs)[grep("motifStart", colnames(tbl.newMotifs))] <- "fp_start"
   colnames(tbl.newMotifs)[grep("motifEnd", colnames(tbl.newMotifs))] <- "fp_end"

   tbl.newMotifs$score1 <- tbl.newMotifs$score
   tbl.newMotifs$score2 <- tbl.newMotifs$qValue
   tbl.newMotifs$score3 <- tbl.newMotifs$pValue
   tbl.newMotifs$score4 <- NA
   tbl.newMotifs$score5 <- NA
   tbl.newMotifs$score6 <- NA
   tbl.newMotifs$sample_id <- NA

   tbl.newMotifs <- tbl.newMotifs[,hits.coi]
   new.order <- order(tbl.newMotifs$loc)
   tbl.newMotifs <- tbl.newMotifs[new.order,]

   locs <- sort(unique(tbl.newMotifs$loc))
   x <- lapply(locs, parseChromLocString)
   tbl.newRegions <- do.call(rbind, lapply(lapply(locs, parseChromLocString), as.data.frame))
   tbl.newRegions$chrom <- as.character(tbl.newRegions$chrom)
   tbl.newRegions$loc <- locs
   tbl.newRegions <- tbl.newRegions[, c("loc", "chrom", "start", "end")]
   colnames(tbl.newRegions)[grep("^end$", colnames(tbl.newRegions))] <- "endpos"

   invisible(list(regions=tbl.newRegions, hits=tbl.newMotifs))


} # callMotifs
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_createFootprintTablesForRegion_small()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_createFootprintTablesForRegion_small <- function()
{
   printf("--- test_createFootprintTablesForRegion_small")

   # hbb.upstream <- 5,219,057
   # tbl.lcr <- data.frame(chrom="chr11", start=5269925, end=5304186, stringsAsFactors=FALSE)
   chrom <- "chr11"
   start.loc <- 5200000
   end.loc   <- 5400000
   x <- createFootprintTablesForRegion(chrom, start.loc, end.loc)
   checkTrue(is.list(x))
   checkEquals(sort(names(x)), c("hits", "regions"))
   checkTrue(all(x$hits$loc %in% x$regions$loc))
      # the regions should be unique
   tbl.regions <- x$regions
   tbl.hits <- x$hits

   checkEquals(length(tbl.regions$loc), length(unique(tbl.regions$loc)))
   checkTrue(nrow(tbl.hits) > nrow(tbl.regions))

   checkEquals(colnames(tbl.regions), c("loc", "chrom", "start", "endpos"))
   smallest.region <- min(1 + tbl.regions$endpos - tbl.regions$start)
   biggest.region <- max(1 + tbl.regions$endpos - tbl.regions$start)
   checkEquals(smallest.region, 8)
   checkEquals(biggest.region, 30)
   
   expected <- c("loc", "fp_start", "fp_end", "type", "name", "length", "strand", "sample_id",
                 "method", "provenance", "score1", "score2", "score3", "score4", "score5", "score6")
   checkEquals(colnames(tbl.hits), expected)

   smallest.region <- min(1 + tbl.hits$fp_end - tbl.hits$fp_start)
   biggest.region <- max(1 + tbl.hits$fp_end - tbl.hits$fp_start)
   checkEquals(smallest.region, 8)
   checkEquals(biggest.region, 30)

   write.table(tbl.regions, file="regions-chr11-small.tsv", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
   write.table(tbl.hits, file="motifs-chr11-small.tsv", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


} # test_createFootprintTablesForRegion_small
#------------------------------------------------------------------------------------------------------------------------
test_createFootprintTablesForRegion_chr11_all <- function()
{
   printf("--- test_createFootprintTablesForRegion_chr11_all")

   x <- createFootprintTablesForRegion("chr11", 1, 500*10e6)
   checkTrue(is.list(x))
   checkEquals(sort(names(x)), c("hits", "regions"))
   checkTrue(all(x$hits$loc %in% x$regions$loc))
      # the regions should be unique
   tbl.regions <- x$regions
   tbl.hits <- x$hits

   checkEquals(length(tbl.regions$loc), length(unique(tbl.regions$loc)))
   checkTrue(nrow(tbl.hits) > nrow(tbl.regions))

   checkEquals(colnames(tbl.regions), c("loc", "chrom", "start", "endpos"))
   expected <- c("loc", "fp_start", "fp_end", "type", "name", "length", "strand", "sample_id",
                 "method", "provenance", "score1", "score2", "score3", "score4", "score5", "score6")
   checkEquals(colnames(tbl.hits), expected)

   save(x, file="chr11.all.regionsAndHits.RData")
   load("chr11.all.regionsAndHits.RData")
   tbl.regions <- x$regions
   tbl.hits <- x$hits
   write.table(tbl.regions, file="regions-chr11.tsv", sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
   write.table(tbl.hits, file="motifs-chr11.tsv", sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)

} # test_createFootprintTablesForRegion_chr11_all
#------------------------------------------------------------------------------------------------------------------------

