source("fimoByChromosome.R")
chromosome.lengths <- org.Hs.egCHRLENGTHS[c(as.character(1:22), "X", "Y")]
names(chromosome.lengths) <- sprintf("chr%s", names(chromosome.lengths))
tbl.regions <- data.frame(chrom="chr3", start=1, end=chromosome.lengths[["chr3"]], stringsAsFactors=FALSE)
tbl.regions <- data.frame(chrom="chr3", start=128484511, end=128492184, stringsAsFactors=FALSE)
fastaFileName <- "gata2.fa"
createFastaFileForFimo(tbl.regions, fastaFileName)
resultsDirectory <- "./"
runFimo(fastaFileName, resultsDirectory, threshold=5e-4)
