source("~/github/fptodb/R/fptodb.R")
#------------------------------------------------------------------------------------------------------------------------
test_small <- function()
{
   betaGlobinLocus <- list(chrom="chr11", start=5223904, end=5304186)
   with(betaGlobinLocus, printf("span: %d", 1 + end - start))
   fpDirectory <- "/local/trenaPlacentaProject/whole_genome"
   f <- with(betaGlobinLocus, createFootprintTablesForRegion(fpDirectory, chrom, start, end))
   resultsDirectory <- tempdir()
   fasta.filename <- "chr11-small.fa"
   createFastaFileForFimo(f, fasta.filename)
   threshold <- 1e-6
   fimoResultsFile <- runFimo(fasta.filename, resultsDirectory, threshold)
   tbl.fimo <- read.table(fimoResultsFile, sep="\t", header=TRUE)
   x <- transformFimoResultsForDatabaseFill(fimoResultsFile)
   lapply(x, dim)
   tbl.regions <- x$regions
   tbl.motifs <- x$hits

   checkTrue(all(tbl.motifs$loc %in% tbl.regions$loc))
   checkTrue(all(tbl.motifs$score3 <= threshold))
   
} # test_small
#------------------------------------------------------------------------------------------------------------------------
test_chr19 <- function()
{
   printf("--- test_chr19")

   chr19locus <- list(chrom="chr19", start=0, end=300000000)
   with(chr19locus, printf("span: %d", 1 + end - start))
   fpDirectory <- "/local/trenaPlacentaProject/whole_genome"
   f <- with(chr19locus, createFootprintTablesForRegion(fpDirectory, chrom, start, end))
   resultsDirectory <- tempdir()
   fasta.filename <- "chr19.fa"
   createFastaFileForFimo(f, fasta.filename)
   threshold <- 1e-6
   fimoResultsFile <- runFimo(fasta.filename, resultsDirectory, threshold)
   tbl.fimo <- read.table(fimoResultsFile, sep="\t", header=TRUE)
   x <- transformFimoResultsForDatabaseFill(fimoResultsFile)
   lapply(x, dim)
   tbl.regions <- x$regions
   tbl.motifs <- x$hits

   checkTrue(all(tbl.motifs$loc %in% tbl.regions$loc))
   checkTrue(all(tbl.motifs$score3 <= threshold))

   save(x, file="chr19-regions-motifs.RData")

} # test_chr19
#------------------------------------------------------------------------------------------------------------------------
test_parseLocStrings <- function()
{
    printf("--- test_parseLocStrings")
    locStrings.3 <- c("chr1:181559-181629", "chr1:181133-181224", "chr1:181267-181349")
    locStrings.20 <- c("chr1:185730-185768", "chr1:10478-10560", "chr1:181005-181123",
                       "chr1:267936-268092", "chr1:186471-186508", "chr1:183290-183357",
                       "chr1:181352-181386", "chr1:183231-183288", "chr1:138493-138527",
                       "chr1:17464-17516", "chr1:184396-184454", "chr1:184473-184510",
                       "chr1:201156-201188", "chr1:10102-10145", "chr1:263670-263703",
                       "chr1:16197-16304", "chr1:183774-183836", "chr1:182816-182852",
                       "chr1:264654-264693", "chr1:10414-10449")
    parseLocStrings(locStrings.3)
    parseLocStrings(locStrings.20)

} # test_parseLocStrings
#------------------------------------------------------------------------------------------------------------------------
test_chr1 <- function()
{
   printf("--- test_chr1")

   resultsDirectory <- "chr1"
   chr1locus <- list(chrom="chr1", start=0, end=300000000)
   with(chr1locus, printf("span: %d", 1 + end - start))

   fpDirectory <- "/local/trenaPlacentaProject/whole_genome"
   f <- with(chr1locus, createFootprintTablesForRegion(fpDirectory, resultsDirectory, chrom, start, end))
   fasta.filename <- file.path(resultsDirectory, "chr1.fa")
   createFastaFileForFimo(f, fasta.filename)
   threshold <- 1e-4
   fimoResultsFile <- runFimo(fasta.filename, resultsDirectory, threshold)
   tbl.fimo <- read.table(fimoResultsFile, sep="\t", header=TRUE)
   dim(tbl.fimo)
   x <- transformFimoResultsForDatabaseFill(fimoResultsFile)
   lapply(x, dim)
   tbl.regions <- x$regions
   tbl.motifs <- x$hits

   checkTrue(all(tbl.motifs$loc %in% tbl.regions$loc))
   checkTrue(all(tbl.motifs$score3 <= threshold))

   save(x, file=file.path("chr1", "chr1-regions-motifs.RData"))
   write.table(x$regions, file=file.path("chr1", "regions.tsv"), col.names=FALSE, row.names=FALSE,
               quote=FALSE, sep="\t")
   write.table(x$hits, file=file.path("chr1", "hits.tsv"), col.names=FALSE, row.names=FALSE,
               quote=FALSE, sep="\t")
   


} # test_chr1
#------------------------------------------------------------------------------------------------------------------------
test_chr2 <- function()
{
   printf("--- test_chr2")

   resultsDirectory <- "chr2"
   chr2locus <- list(chrom="chr2", start=0, end=300000000)
   with(chr2locus, printf("span: %d", 1 + end - start))

   fpDirectory <- "/local/trenaPlacentaProject/whole_genome"
   f <- with(chr2locus, createFootprintTablesForRegion(fpDirectory, resultsDirectory, chrom, start, end))
   fasta.filename <- file.path(resultsDirectory, "chr2.fa")
   createFastaFileForFimo(f, fasta.filename)
   threshold <- 1e-4
   fimoResultsFile <- runFimo(fasta.filename, resultsDirectory, threshold)
   tbl.fimo <- read.table(fimoResultsFile, sep="\t", header=TRUE)
   dim(tbl.fimo)
   x <- transformFimoResultsForDatabaseFill(fimoResultsFile)
   lapply(x, dim)
   tbl.regions <- x$regions
   tbl.motifs <- x$hits

   checkTrue(all(tbl.motifs$loc %in% tbl.regions$loc))
   checkTrue(all(tbl.motifs$score3 <= threshold))

   save(x, file=file.path("chr2", "chr2-regions-motifs.RData"))
   write.table(x$regions, file=file.path("chr2", "regions.tsv"), col.names=FALSE, row.names=FALSE,
               quote=FALSE, sep="\t")
   write.table(x$hits, file=file.path("chr2", "hits.tsv"), col.names=FALSE, row.names=FALSE,
               quote=FALSE, sep="\t")
   

} # test_chr2
#------------------------------------------------------------------------------------------------------------------------
test_chr3 <- function()
{
   printf("--- test_chr3")

   resultsDirectory <- "chr3"
   chr3locus <- list(chrom="chr3", start=0, end=300000000)
   with(chr3locus, printf("span: %d", 1 + end - start))

   fpDirectory <- "/local/trenaPlacentaProject/whole_genome"
   f <- with(chr3locus, createFootprintTablesForRegion(fpDirectory, resultsDirectory, chrom, start, end))
   fasta.filename <- file.path(resultsDirectory, "chr3.fa")
   createFastaFileForFimo(f, fasta.filename)
   threshold <- 1e-4
   fimoResultsFile <- runFimo(fasta.filename, resultsDirectory, threshold)
   tbl.fimo <- read.table(fimoResultsFile, sep="\t", header=TRUE)
   dim(tbl.fimo)
   x <- transformFimoResultsForDatabaseFill(fimoResultsFile)
   lapply(x, dim)
   tbl.regions <- x$regions
   tbl.motifs <- x$hits
   checkTrue(all(colnames(tbl.motifs) == motif.hits.coi()))

   checkTrue(all(tbl.motifs$loc %in% tbl.regions$loc))
   checkTrue(all(tbl.motifs$score3 <= threshold))

   save(x, file=file.path("chr3", "chr3-regions-motifs.RData"))
   write.table(x$regions, file=file.path("chr3", "regions.tsv"), col.names=FALSE, row.names=FALSE,
               quote=FALSE, sep="\t")
   write.table(x$hits, file=file.path("chr3", "hits.tsv"), col.names=FALSE, row.names=FALSE,
               quote=FALSE, sep="\t")

} # test_chr3
#------------------------------------------------------------------------------------------------------------------------
runChromosome <- function(chromName, length=300000000)
{
   printf("--- runChromosome %s", chromName)

   resultsDirectory <- chromName
   if(!file.exists(resultsDirectory))
       dir.create(resultsDirectory)
   
   chromlocus <- list(chrom=chromName, start=0, end=length)
   with(chromlocus, printf("span: %d", 1 + end - start))

   fpDirectory <- "/local/trenaPlacentaProject/whole_genome"
   f <- with(chromlocus, createFootprintTablesForRegion(fpDirectory, resultsDirectory, chrom, start, end))
   fasta.filename <- file.path(resultsDirectory, sprintf("%s.fa", chromName))
   
   createFastaFileForFimo(f, fasta.filename)
   threshold <- 1e-4
   fimoResultsFile <- runFimo(fasta.filename, resultsDirectory, threshold)
   tbl.fimo <- read.table(fimoResultsFile, sep="\t", header=TRUE)
   dim(tbl.fimo)
   x <- transformFimoResultsForDatabaseFill(fimoResultsFile)
   lapply(x, dim)
   tbl.regions <- x$regions
   tbl.motifs <- x$hits
   checkTrue(all(colnames(tbl.motifs) == motif.hits.coi()))

   checkTrue(all(tbl.motifs$loc %in% tbl.regions$loc))
   checkTrue(all(tbl.motifs$score3 <= threshold))

   save(x, file=file.path(resultsDirectory, sprintf("%s-regions-motifs.RData", chromName)))
   write.table(x$regions, file=file.path(resultsDirectory, "regions.tsv"), col.names=FALSE, row.names=FALSE,
               quote=FALSE, sep="\t")
   write.table(x$hits, file=file.path(resultsDirectory, "hits.tsv"), col.names=FALSE, row.names=FALSE,
               quote=FALSE, sep="\t")

} # runChromosome
#------------------------------------------------------------------------------------------------------------------------
