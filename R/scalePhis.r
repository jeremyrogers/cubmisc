phi <- "../jeremyYeast/jyeast.phi.tsv"
output <- "../jeremyYeast/beyer.scaled.phi.tsv"

orig <- read.table(phi, sep=',', header=FALSE)
m <- mean(orig[ ,2])

scaled <- orig[ ,2]/m

#cat(as.character(orig[ ,1]));cat("\n");

mat <- matrix(c(as.character(orig[ ,1]), scaled), ncol=2)
cat(mat)
write.table(mat, file=output, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
