nsepr <- read.table("../jeremyYeast/codon/NSE_pr/modded_NsePr.tsv", sep="\t", header=FALSE)
output <- "../jeremyYeast/codon/Elong_rate/modded_elong.tsv"

B <- 0.0025

omega <- (B / nsepr[ ,3]) - B

omega_mat <- matrix(c(as.character(nsepr[ ,1]), as.character(nsepr[ ,2]), omega), ncol=3)
write.table(omega_mat, file=output, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) 
