##### Script to divide the genome into (ratio) distinct parts ######
random5ths <- function(fasta, phifile, ratio = 5, howMany=NULL, filename=NULL, sort=TRUE){
source("writing.r");
genome <- read.seq(fasta);
phi <- read.csv(phifile)
if(is.null(howMany)){ howMany <- ratio }
if(is.null(filename)){ filename <- strsplit(fasta, ".fasta")[[1]] }

sectionlength <- round(length(genome)/ratio);

for(j in 1:howMany)
{

section <- list();
sectionphi <- matrix(ncol = 2, nrow = sectionlength);
colnames(sectionphi) <- c("Gene", "Xobs")

for(i in 1:(sectionlength) ) #genome loop
{
  rand <- round(runif(1,1,length(genome)));
  
  section[i] <- genome[rand];
  names(section)[i] <- names(genome[rand])
  
  sectionphi[i,1] <- names(genome[rand]);
  sectionphi[i,2] <- phi$Xobs[rand];
  
  genome <- genome[-rand];
  phi <- phi[-rand,]
}

if(sort){
  sectionSort <- order(names(section))
  section <- section[sectionSort]
  sectionphi <- sectionphi[sectionSort,]
}

write.seq(section, paste(filename,"-section",j,"of",ratio,".fasta",sep=""))
write.csv(sectionphi, paste(filename,"-section",j,"of",ratio,".csv",sep=""), quote=FALSE, row.names=FALSE)
}

}