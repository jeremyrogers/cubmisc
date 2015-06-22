source("config.r")

input <- "../jeremyYeast/codon/Elong_rate/modded_elong.tsv"
output <- "../jeremyYeast/codon/Elong_rate/modded_delta_omega.tsv"

t <- read.table(input, sep="\t", header=FALSE)

aminoacids <- t[ ,1]

new <- c()

for(aa in config$aa) {
	indices <- grep(aa, aminoacids)
	if(length(indices) != 1) {
		value <- t[indices[length(indices)],3]
		for(index in 1:(length(indices)-1)) {
			new[[length(new)+1]] <- aa
			new[[length(new)+1]] <- as.character(t[indices[index],2])
			new[[length(new)+1]] <- t[indices[index],3]-value
		}
	}
}

mat <- matrix(new, ncol=3, byrow=TRUE)

write.table(mat, file=output, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
