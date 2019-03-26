source("findMotifsByChromosome.R")
#runTests()
chromosome <- "chr1"
x <- createFootprintTablesForRegion(chromosome, 0, 300000000)
# checkTrue(all(x$hits$loc %in% x$regions$loc))
# tbl.regions <- x$regions
# tbl.motifs <- x$hits
# locs.file.name <- sprintf("%s-locs.tsv", chromosome)
# motifs.file.name <- sprintf("%s-motifs.tsv", chromosome)
# write.table(tbl.regions, file=locs.file.name,   sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
# write.table(tbl.motifs,  file=motifs.file.name, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
# 
