phi <- read.csv("pyeast.phi.tsv")

phi.noise10 <- phi;
phi.noise25 <- phi;
phi.noise50 <- phi;
phi.unoise100 <- phi;
phi.unoise200 <- phi;
phi.unoise300 <- phi;
phimean <- mean(phi[,2]);

for(i in 1:length(phi[,2])){
  phi.noise10[i,2] <- rnorm(1, phi[i,2], phi[i,2]/10)
  phi.noise25[i,2] <- rnorm(1, phi[i,2], phi[i,2]/4)
  phi.noise50[i,2] <- rnorm(1, phi[i,2], phi[i,2]/2)
  phi.unoise100[i,2] <- rnorm(1, phi[i,2], phimean*1)
  phi.unoise200[i,2] <- rnorm(1, phi[i,2], phimean*2)
  phi.unoise300[i,2] <- rnorm(1, phi[i,2], phimean*3)
}

write.table(phi.noise10, file="lyeast.phi.noise10.tsv", quote=FALSE, row.names=FALSE, sep=",")
write.table(phi.noise25, file="lyeast.phi.noise25.tsv", quote=FALSE, row.names=FALSE, sep=",")
write.table(phi.noise50, file="lyeast.phi.noise50.tsv", quote=FALSE, row.names=FALSE, sep=",")
write.table(phi.unoise100, file="lyeast.phi.unoise100.tsv", quote=FALSE, row.names=FALSE, sep=",")
write.table(phi.unoise200, file="lyeast.phi.unoise200.tsv", quote=FALSE, row.names=FALSE, sep=",")
write.table(phi.unoise300, file="lyeast.phi.unoise300.tsv", quote=FALSE, row.names=FALSE, sep=",")