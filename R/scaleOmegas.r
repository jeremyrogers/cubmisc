omega <- "../prestonYeast/scaled_omega.csv"
output <- "../jeremyYeast/codon/Elong_rate/scaled_omega.tsv"

orig <- read.table(omega, sep="\t")

m <- mean(orig[ ,3])

scaled <- orig[ ,3]/m


cat(mean(scaled));cat("\n")
#cat(as.character(orig[ ,1]));cat("\n");

mat <- matrix(c(as.character(orig[ ,1]), as.character(orig[ ,2]), scaled), ncol=3)
write.table(mat, file=output, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
