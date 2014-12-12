##### Script to divide the genome into (ratio) distinct parts ######
source("writing.r");
filename <- "S.cerevisiae.S288c.REU.sim.b-0.001.ces.fasta";
ratio <- 5;
genome <- read.seq(filename);

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
  names(section)[i] <- unlist(strsplit(x = names(genome[rand]),split = '\t'))[1];
  
  sectionphi[i,1] <- names(section)[i];
  sectionphi[i,2] <- unlist(strsplit(x = names(genome[rand]),split = '\t'))[3];
  
  genome <- genome[-rand];
}

write.seq(section, paste("section",j,"of",ratio,".fasta",sep=""))
write.csv(sectionphi, paste("section",j,"of",ratio,".csv",sep=""), quote=FALSE, row.names=FALSE)
}
