##### Script to divide the genome into (ratio) distinct parts ######
random5ths <- function(fasta, phifile, ratio = 5){
source("writing.r");
genome <- read.seq(fasta);
phi <- read.csv(phifile)

sectionlength <- length(genome)/ratio;

for(j in 1:ratio)
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
}

write.seq(section, paste("section",j,"of",ratio,".fasta",sep=""))
write.csv(sectionphi, paste("section",j,"of",ratio,".csv",sep=""), quote=FALSE, row.names=FALSE)
}

}