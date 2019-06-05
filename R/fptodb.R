library(RPostgreSQL)
library(RUnit)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("parseChromLocString"))
   source("~/github/trena/R/utils.R")
#------------------------------------------------------------------------------------------------------------------------
motif.hits.coi <- function()
{
   return( c("loc", "fp_start", "fp_end", "type", "name", "length", "strand", "sample_id",
             "method", "provenance", "score1", "score2", "score3", "score4", "score5", "score6"))
}
#------------------------------------------------------------------------------------------------------------------------
createFootprintTablesForRegion <- function(fpDirectory, resultsDirectory, chromosome, start.loc, end.loc,
                                           shoulder=10)
{
       #------------------------------------------------------------
       # first read the footprint calls
       #------------------------------------------------------------

    message(sprintf("--- creating footpring tables for specified region"))

    searchString <- sprintf('find %s -name "*%s_fp.bed"', fpDirectory, chromosome)
    files <- system(searchString, intern=TRUE)
    message(sprintf("file count: %d", length(files)))
    stopifnot(length(files) > 0)
    fps <- list()
    for(f in files){
       f.path <- f # file.path(dir, f)
       message(sprintf("--- %s", f.path))
       tbl <- read.table(f.path, sep="\t", as.is=TRUE, nrow=-1)
       colnames(tbl) <- c("chrom", "start", "end", "name", "score", "strand", "ignore")
       tbl <- subset(tbl, chrom==chromosome & start >= start.loc & end <= end.loc)
       fps[[f]] <- tbl
       } # for f

       #------------------------------------------------------------
       # combine the footprint lists into a single data.frame
       #------------------------------------------------------------

    tbl.fps <- do.call(rbind, fps)[, c("chrom", "start", "end", "score")]
    rownames(tbl.fps) <- NULL
    dim(tbl.fps)
    fp.order <- order(tbl.fps$chrom, tbl.fps$start)
    tbl.fps <- tbl.fps[fp.order,]

       #------------------------------------------------------------
       # expand the footprints by <shoulder> in both directions
       #------------------------------------------------------------

    tbl.fpsExpanded <- tbl.fps
    tbl.fpsExpanded$start <- tbl.fpsExpanded$start - shoulder
    tbl.fpsExpanded$end <- tbl.fpsExpanded$end + shoulder

       #------------------------------------------------------------
       # remove any regions which are an exact duplicate of another
       # sort order includes 
       #------------------------------------------------------------

    tbl.fpsExpanded$signature <- with(tbl.fpsExpanded, sprintf("%s:%d-%d", chrom, start, end))
    tbl.fpsExpanded <- tbl.fpsExpanded[with(tbl.fpsExpanded, order(signature, -score)),]
    dups <- which(duplicated(tbl.fpsExpanded$signature))
    if(length(dups) > 0)
       tbl.fpsExpanded <- tbl.fpsExpanded[-dups,]
    filename <- sprintf("%s:%d-%d.footprints-combined.RData", chromosome, start.loc, end.loc)
    message(sprintf("saving %d combined, expanded footprints to %s", nrow(tbl.fpsExpanded), filename))
    full.filename <- file.path(resultsDirectory, filename)
    tbl.fpsExpanded <- tbl.fpsExpanded[, c("chrom", "start", "end", "score")]
    save(tbl.fpsExpanded, file=full.filename)
    return(full.filename)

} # createFootprintTablesForRegion
#------------------------------------------------------------------------------------------------------------------------
createFastaFileForFimo <- function(serializedFileName, fastaFileName)
{
   printf("--- creating fasta file for FIMO")
   tbl.fps <- get(load(serializedFileName))
   sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg38, tbl.fps$chrom, tbl.fps$start, tbl.fps$end)
   names(sequences) <- sprintf("%s:%d-%d", tbl.fps$chrom, tbl.fps$start, tbl.fps$end)
   writeXStringSet(sequences, fastaFileName)

} # createFastaFileForFimo
#------------------------------------------------------------------------------------------------------------------------
runFimo <- function(fastaFileName, resultsDirectory, threshold)
{
   printf("--- running FIMO")
   FIMO <- "/users/pshannon/meme/bin/fimo"
   MOTIFS <- "~/github/fimoService/pfms/human-jaspar2018-hocomoco-swissregulon.meme"
   cmd <- sprintf("%s --oc %s --thresh %f --verbosity 1 %s %s",
                  FIMO, resultsDirectory, threshold, MOTIFS, fastaFileName)
   system(cmd)
   return(file.path(resultsDirectory, "fimo.tsv"))    
   # /users/pshannon/meme/bin/fimo --oc . --thresh 1e-4 ~/github/fimoService/pfms/human-jaspar2018-hocomoco-swissregulon.meme chr11-small.fa 

} # runFimo
#------------------------------------------------------------------------------------------------------------------------
parseLocStrings <- function(locStrings)
{
   match <- regexpr("(?<chromosome>chr.*):(?<startPos>\\d+)-(?<endPos>\\d+)", locStrings, perl=TRUE)
   columns <-     attr(match, "capture.names")
   starts <-  as.data.frame(attr(match, "capture.start"))
   lengths <- as.data.frame(attr(match, "capture.length"))
   chroms <- substring(locStrings, starts$chromosome, starts$chromosome + lengths$chromosome -1)
   start.locs <- as.integer(substring(locStrings, starts$startPos, starts$startPos + lengths$startPos -1))
   end.locs <- as.integer(substring(locStrings, starts$endPos, starts$endPos + lengths$endPos -1))
   data.frame(chrom=chroms, start=start.locs, end=end.locs, stringsAsFactors=FALSE)

} # parseLocStrings
#------------------------------------------------------------------------------------------------------------------------
transformFimoResultsForDatabaseFill <- function(fimoResultsFile)
{
   printf("--- transforming FIMO results for database fill")

   tbl.motifs <- read.table(fimoResultsFile, sep="\t", header=TRUE, stringsAsFactors=FALSE)
   colnames(tbl.motifs)[grep("start", colnames(tbl.motifs))] <- "motifStart"
   colnames(tbl.motifs)[grep("stop", colnames(tbl.motifs))] <- "motifEnd"

        #-----------------------------------------------------------------------------------------------
        # expand "sequence_name" to reproduce footprint start/stop locs, and then motif start/stop locs
        #-----------------------------------------------------------------------------------------------

   printf("    --- parseLocStrings")
   tbl.locs <- parseLocStrings(tbl.motifs$sequence_name)
   tbl.motifs$chrom <- tbl.locs$chrom
   tbl.motifs$motifStart <- tbl.motifs$motifStart + tbl.locs$start
   tbl.motifs$motifEnd <- tbl.motifs$motifEnd + tbl.locs$start
   tbl.motifs$loc <- with(tbl.motifs, sprintf("%s:%d-%d", chrom, motifStart, motifEnd))

   colnames(tbl.locs)[grep("start", colnames(tbl.locs))] <- "fp_start"
   colnames(tbl.locs)[grep("end", colnames(tbl.locs))] <- "fp_end"

   
   printf("    --- cbind")
   tbl.motifs <- cbind(tbl.locs, tbl.motifs)

   #fp.offset <- tbl.motifs$fp_start
        # our database scheme improperly uses fp_start and fp_end for motif start and end.
        # we accomdate that old choice here:
   #tbl.motifs$fp_start <- tbl.motifs$motifStart + fp.offset
   #tbl.motifs$fp_end  <- tbl.motifs$motifEnd + fp.offset
   # tbl.motifs$motifStart <- tbl.motifs$motifStart + tbl.motifs$fp_start
   # tbl.motifs$motifEnd <- tbl.motifs$motifEnd + tbl.motifs$fp_start
   colnames(tbl.motifs)[grep("p.value", colnames(tbl.motifs))] <- "pValue"
   colnames(tbl.motifs)[grep("q.value", colnames(tbl.motifs))] <- "qValue"
   colnames(tbl.motifs)[grep("motif_id", colnames(tbl.motifs))] <- "motif"

   colnames(tbl.motifs)     # "chrom", "fp_start", "fp_end", "motif", "sequence_name", "motifStart",
                            # "motifEnd", "strand", "score", "pValue", "qValue", "matched_sequence"

   new.colnames <- c("loc", "chrom", "motifStart", "motifEnd", "motif", "strand", "score",
                     "pValue", "qValue", "matched_sequence")
   tbl.newMotifs <- tbl.motifs[, new.colnames]
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

   tbl.newMotifs <- tbl.newMotifs[, motif.hits.coi()]
   new.order <- order(tbl.newMotifs$loc)
   tbl.newMotifs <- tbl.newMotifs[new.order,]

   printf("    --- parse newRegions locs")
   locs <- sort(unique(tbl.newMotifs$loc))
   tbl.newRegions <- parseLocStrings(locs)
   tbl.newRegions$chrom <- as.character(tbl.newRegions$chrom)
   tbl.newRegions$loc <- locs
   tbl.newRegions <- tbl.newRegions[, c("loc", "chrom", "start", "end")]
   colnames(tbl.newRegions)[grep("^end$", colnames(tbl.newRegions))] <- "endpos"

   invisible(list(regions=tbl.newRegions, hits=tbl.newMotifs))

} # transformFimoResultsForDatabaseFill
#------------------------------------------------------------------------------------------------------------------------
