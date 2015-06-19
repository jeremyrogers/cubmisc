omega <- "../jeremyYeast/codon/Elong_rate/modded_elong.tsv"
output <- "../jeremyYeast/codon/Elong_rate/modded_elong_scaled.tsv"

orig <- read.table(omega, sep="\t")

m <- mean(orig[ ,3])

scaled <- orig[ ,3]/m


cat(mean(scaled));cat("\n")
#cat(as.character(orig[ ,1]));cat("\n");

mat <- matrix(c(as.character(orig[ ,1]), as.character(orig[ ,2]), scaled), ncol=3)
write.table(mat, file=output, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
