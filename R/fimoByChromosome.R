library(RUnit)
library(GenomicRanges)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(MotifDb)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
  test_createFastaFileForFimo()
  test_runFimo()
  
} # runTests
#------------------------------------------------------------------------------------------------------------------------
createFastaFileForFimo <- function(tbl.regions, fastaFileName)
{
   sequences <- with(tbl.regions, getSeq(BSgenome.Hsapiens.UCSC.hg38, chrom, start, end))
   if(is(sequences, "DNAString"))
       sequences <- DNAStringSet(sequences)
   names(sequences) <- with(tbl.regions, sprintf("%s:%d-%d", chrom, start, end))
   message(sprintf("creating fasta file '%s' for %d sequences, for FIMO", fastaFileName, length(sequences)))

   writeXStringSet(sequences, fastaFileName)

} # createFastaFileForFimo
#------------------------------------------------------------------------------------------------------------------------
runFimo <- function(fastaFileName, resultsDirectory, threshold=5e-4)
{
   printf("--- running FIMO")
   FIMO <- "/users/pshannon/meme/bin/fimo"
   MOTIFS <- "~/github/fimoService/pfms/human-jaspar2018-hocomoco-swissregulon.meme"
   cmd <- sprintf("%s --oc %s --thresh -%f --text --verbosity 1 %s %s",
                  FIMO, resultsDirectory, threshold, MOTIFS, fastaFileName)
   print(cmd)
   system(cmd)
   return(file.path(resultsDirectory, "fimo.tsv"))    
   # /users/pshannon/meme/bin/fimo --oc . --thresh 1e-4 ~/github/fimoService/pfms/human-jaspar2018-hocomoco-swissregulon.meme chr11-small.fa 

} # runFimo
#------------------------------------------------------------------------------------------------------------------------
test_createFastaFileForFimo <- function()
{
   message(noquote(sprintf("--- test_createFastaFileForFimo")))
    
     # a < 1kb region in the promoter of GATA2, where TBX15 hits may be found
   tbl.regions <- data.frame(chrom="chr3", start=128497569, end=128498329, stringsAsFactors=FALSE)
   createFastaFileForFimo(tbl.regions, "smallTest.fa")
   checkTrue(file.exists("smallTest.fa"))
    
} # test_createFastaFileForFimo
#------------------------------------------------------------------------------------------------------------------------
test_runFimo <- function()
{
   message(noquote(sprintf("--- test_runFimo")))

   fastaFileName <- "smallTest.fa"  # created in test_createFastaFileforFimo
   checkTrue(file.exists(fastaFileName))
   resultsDirectory <- "./fimoResults"
   runFimo(fastaFileName, resultsDirectory, threshold=1e-2)   
   checkTrue(file.exists(file.path(resultsDirectory, "fimo.tsv")))

} # test_runFimo
#------------------------------------------------------------------------------------------------------------------------
runFimoGATA2.big <- function()
{
      # GATA2 and the full span of all enhancers
   tbl.regions <- data.frame(chrom="chr3", start=128013674, end=128712841, stringsAsFactors=FALSE)
      # only 263kb around GATA2's enhancers
   tbl.regions <- data.frame(chrom="chr3", start=128383794, end=128647775, stringsAsFactors=FALSE)
   with(tbl.regions, printf("span: %d", end-start))
   fastaFilename <- tempfile(fileext=".fa")
   createFastaFileForFimo(tbl.regions, fastaFilename)
   checkTrue(file.exists(fastaFilename))
   checkTrue(file.size(fastaFilename) > with(tbl.regions, end-start))  # bigger than the base count

   resultsDirectory <- tempdir()
   resultsDirectory <- "./fimoResults-10e-3-chr3-128383794-128647775-v2"
      #  4 minutes for 1e-4, 263k    201090 lines
      #  5 minutes for 1e-3, 263k   1510766 lnes
      # 18 minutes for 1e-2, 263k  12033436 lines

   system.time(runFimo(fastaFilename, resultsDirectory, threshold=1e-3))
   fimo.results.file <- file.path(resultsDirectory, "fimo.tsv")
   checkTrue(file.exists(fimo.results.file))
    
   tbl <- read.table(fimo.results.file, sep="\t", as.is=TRUE, nrow=-1, header=TRUE)  # two chopped names
   tbl.fixed <- fixMotifNamesTruncatedAt100characters(tbl)
   fimo.results.file.fixed <- file.path(resultsDirectory, "fimo-fixed.tsv")
   write.table(tbl.fixed, fimo.results.file.fixed, sep="\t", quote=FALSE, row.names=FALSE)

} # runFimoGATA2.big
#------------------------------------------------------------------------------------------------------------------------
fixMotifNamesTruncatedAt100characters <- function(tbl)
{
    char.100.lines <- which(nchar(tbl$motif_id) == 100)
    motifDb.indices <- lapply(tbl$motif_id[char.100.lines], function(id) grep(id, names(MotifDb)))
    motifDb.names <- names(MotifDb)[as.integer(motifDb.indices)]
    tbl$motif_id[char.100.lines] <- motifDb.names

    invisible(tbl)

} # fixMotifNamesTruncatedAt100characters
#------------------------------------------------------------------------------------------------------------------------
test_fixMotifNamesTruncatedAt100characters <- function()
{
   printf("--- test_fixMotifNamesTruncatedAt100characters")

   resultsDirectory <- "fimoResults-10e-3-chr3-128383794-128647775"
   filename <- file.path(resultsDirectory, "fimo.tsv")
   stopifnot(file.exists(filename))
   tbl <- read.table(filename, sep="\t", as.is=TRUE, nrow=1500, header=TRUE)  # two chopped names

   char.100.lines <- which(nchar(tbl$motif_id) == 100)   # 1092 1477

   checkTrue(length(char.100.lines)  > 0)
   tbl.fixed <- fixMotifNamesTruncatedAt100characters(tbl)
   char.100.lines.after <- which(nchar(tbl.fixed$motif_id) == 100)
   checkTrue(length(char.100.lines.after) == 0)
   checkTrue(all(nchar(tbl.fixed$motif_id[char.100.lines]) > 100))
    
} # test_fixMotifNamesTruncatedAt100characters
#------------------------------------------------------------------------------------------------------------------------
# if(!interactive()) test_runFimoBig()
