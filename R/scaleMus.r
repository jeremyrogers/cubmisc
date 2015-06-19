mu <- "../jeremyYeast/codon/Mut_rate/PNAS2011_S.cer.mut.tsv"
output <- "../jeremyYeast/codon/Mut_rate/PNAS2011_scaled.tsv"

orig <- read.table(logmu, sep="\t")
orig
m <- mean(orig[ ,3])

scaled <- orig[ ,3]/m


cat(mean(scaled));cat("\n")
#cat(as.character(orig[ ,1]));cat("\n");

mat <- matrix(c(as.character(orig[ ,1]), as.character(orig[ ,2]), scaled), ncol=3)
write.table(mat, file=output, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
