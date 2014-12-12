
#copied directly from the cubfits package data.io.r
read.seq <- function(file.name, forceDNAtolower = FALSE,
    convertDNAtoupper = TRUE){
  ret <- seqinr::read.fasta(file.name, forceDNAtolower = forceDNAtolower)
        
  ### Make sure everything is in upper case.
  if(convertDNAtoupper){ 
    ret <- lapply(ret, function(x){ dna.low2up(x) })
  } 

  ret
} # End of read.seq().



#copied directly from the cubfits package data.io.r
write.seq <- function(seq.data, file.name){
  seqinr::write.fasta(seq.data, names(seq.data), file.name)
  invisible()
} # End of write.seq().




#copied directly from the cubfits package data.codon.convert.r
### This converts DNA lower case to upper case.
dna.low2up <- function(x){
  n.set <- c("a", "c", "g", "t")
  N.set <- c("A", "C", "G", "T")
  for(i in 1:4){
    x[x == n.set[i]] <- N.set[i]
  }
  x
} # End of dna.low2up().

