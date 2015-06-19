#library(cubfits, quietly = TRUE, lib.loc = "~/cubfitsNSEdebug/")
library(cubfits, quietly = TRUE, lib.loc = "~/cubfitsBuild/")
#setwd("~/cubmisc/R/")
source("visualize_utility.r")
source("run_utility.r")
source("config.r")

result.folder <- "../results/jeremyyeast/scaledphi3/"
results <- "../results/jeremyyeast/scaledphi3/jeremyyeastwithPhi.dat"
results.phi <- "../results/jeremyyeast/scaledphi3/jeremyyeastwithPhi.phi"
read.genome <- "../jeremyYeast/Genome/jeremyYeast.fasta"
empirical.data <- "../jeremyYeast/jyeast.scaled.phi.tsv"
simulated_data <- TRUE;

if(simulated_data) {
	omega_true <- read.table("../jeremyYeast/codon/Elong_rate/modded_elong_scaled.tsv", sep="\t", header=FALSE)
	mu_true <- read.table("../jeremyYeast/codon/Mut_rate/PNAS2011_scaled.tsv", sep="\t", header=FALSE)
}

model <- 'nse';              #nse \ roc

genome <- 'jeremyYeast';       #ecoli \ REUyeast \ pYeast \ brewYeast \ loganYeast
prefix <- "06-19"            #generally, date the run started in "MM-DD" format
suffix <- "Scaled_phi"            #What was special about this run?

delta.a_12 <- 0
a_2 <- 4

load(results)
estm.phi <- read.csv(results.phi)
seq.string <- readGenome(read.genome, config$rm.short, config$rm.first.aa)
emp <- read.empirical.data(empirical.data, seq.string, config$selected.env, th=0)
cat("read data\n")
#names(chain$b.Mat[[1]]) <- read.table("names.txt")[[1]]

Prefix <- suffix;
suffix <- prefix;
prefix <- Prefix;
rm(Prefix)

emp <- emp$empirical[names(emp$empirical) %in% estm.phi[, 1]]

if(simulated_data){
  require(package = plotrix, lib.loc = "~/cubfits/Dependencies/")
  chunksize <- floor(length(chain$b.Mat)/100)*10
  ci <- .95 #Confidence interval, on a scale of (0,1). e.g. 0.95 means 95%CI
  
  
  #average converged omega
{
	write.table(omega_true, file="", quote=FALSE)
  omega <- matrix(0, nrow=length(omega_true[[3]]), ncol=4)
  colnames(omega) <- c("'true' deltaomega", "estimated deltaomega"
                       ,paste("Lower ", ci, "%ci", sep="") , paste("Upper ", ci, "%ci", sep="")
  )
  omega[,1] <- omega_true[[3]];
  #vec <- grep("eta", names(chain$b.Mat[[1]]))
  vec <- grep("omega", names(chain$b.Mat[[1]]))
  tmpomega <- matrix(ncol=chunksize, nrow=length(vec));
  
  for(ii in 1:length(vec)){
    for(jj in 1:chunksize){
      tmpomega[ii,jj] <- chain$b.Mat[[jj+length(chain$b.Mat)-chunksize]][vec[ii]]
    }##finish current codon
    omega[ii,2] <- mean(tmpomega[ii,])
    omega[ii,3:4] <- quantile(tmpomega[ii,], probs = c( (1-ci)/2 , (1+ci)/2))
  }##finish all codons
  
}
#omega[,1] is the 'true' omega and omega[,2] is the estimated omega
#omega[,3] and [,4] are the lower and upper 95% confidence interval

#average converged mu
{
	write.table(mu_true, file="");
  mu <- matrix(0, nrow=length(mu_true[[3]]), ncol=4)
  colnames(mu) <- c("'true' deltaM", "estimated deltaM"
                    ,paste("Lower ", ci, "%ci", sep="") , paste("Upper ", ci, "%ci", sep="")
  )
#  mu[,1] <- -1 * mu_true[[3]];
  mu[,1] <- mu_true[[3]];
  #The mu values used to generate preston's yeast are flipped in the code due to the deltaEta switch, which he wasn't there for  
  vec <- grep("mu", names(chain$b.Mat[[1]]))
  tmpmu <- matrix(ncol=chunksize, nrow=length(vec))
  for(ii in 1:length(vec)){
    for(jj in 1:chunksize){
      #mu[ii,2] <- mu[ii,2]+ chain$b.Mat[[jj]][vec[ii]]
      tmpmu[ii,jj] <- chain$b.Mat[[jj+length(chain$b.Mat)-chunksize]][vec[ii]]
    }##finish current codon
    mu[ii,2] <- mean(tmpmu[ii,])
    mu[ii,3:4] <- quantile(tmpmu[ii,], probs = c( (1-ci)/2 , (1+ci)/2))
  }##finish all codons
  
}


#### PLOT RESULTS VS 'TRUE' VALUES
# ?
#Plot Omega
pdf(paste(result.folder, prefix, "_codonParameters_", suffix, ".pdf", sep = ""), width = 8, height = 10)
par(mfrow=c(2,1))
reg <- lm(omega[,2]~omega[,1]);
plot(omega, main=expression(paste("'true' ",Delta,omega," vs Estimated ",Delta,omega))
     #,ylim=range(omega)
     ,ylim=range(omega[,3:4])
     ,xlab=expression(paste("'true' ", Delta, omega))
     ,ylab=expression(paste("estimated ", Delta, omega))
     ,sub=paste("R^2=", format(summary(reg)$r.squared, digits = 3), sep="") 
)
abline(reg, lwd=2, col="red")
abline(0,1, lwd=2, col="lightslategrey")
plotCI(x = omega[,1], y=omega[,2], ui=omega[,4], li=omega[,3], add=TRUE, scol = "darkslategrey");
legend("bottomright", legend = c("1-1 line", "Linear Fit Line", paste(100*ci, "% confidence", sep=""))
       ,fill=c("lightslategrey", "red", "darkslategrey")
)

###Plot Mu
reg <- lm(mu[,2]~mu[,1]);
plot(mu, main=expression(paste("'true' ",Delta,"M vs Estimated ",Delta,"M"))
     #,ylim=range(mu)
     ##,ylim=range(mu[,3:4])
     ,xlab=expression(paste("'true' ", Delta, "M"))
     ,ylab=expression(paste("estimated ", Delta,"M"))
     ,sub=paste("R^2=", format(summary(reg)$r.squared, digits = 3), sep="")
)
abline(0,1, lwd=2, col="lightslategrey")
abline(reg, lwd=2, col="red")
plotCI(x = mu[,1], y=mu[,2], ui=mu[,4], li=mu[,3], add=TRUE, scol = "darkslategrey");
legend("bottomright", legend = c("1-1 line", "Linear Fit Line", paste(100*ci, "% confidence", sep=""))
       ,fill=c("lightslategrey", "red", "darkslategrey")
)

dev.off()
}

#### PLOT TRACES
{
pdf(paste(result.folder, prefix, "_logmu_", suffix, ".pdf", sep = ""), width = 12, height = 11)
#plotTraces(data.set$res[[1]]$b.Mat, config$aa, param="logmu", main=paste("AA parameter trace ", prefix, sep=""))
plotTraces(chain$b.Mat, config$aa, param="logmu", main=paste("AA parameter trace ", prefix, sep=""))
dev.off()
 
if(model=="nse"){ symbol <- "omega" } else{ symbol <- "eta" }
pdf(paste(result.folder, prefix, "_delta", symbol, "_", suffix, ".pdf", sep = ""), width = 12, height = 11)
#plotTraces(data.set$res[[1]]$b.Mat, config$aa, param="deltaeta", main=paste("AA parameter trace ", prefix, sep=""))
plotTraces(chain$b.Mat, config$aa, param=paste("delta", symbol, sep=""), main=paste("AA parameter trace ", prefix, sep=""))
dev.off()

interval <- (length(chain$b.Mat)-config$use.n.samples):length(chain$b.Mat)
pdf(paste(result.folder, prefix, "_delta", symbol, "_histogram_", suffix, ".pdf", sep = ""), width = 12, height = 11)
plotBMatrixPosterior(chain$b.Mat, config$aa, interval, param=paste("delta", symbol, sep=""), main=paste("AA parameter histogram", prefix, sep=""), nclass=100, center=T)
dev.off()

pdf(paste(result.folder, prefix, "_logmu_histogram_", suffix, ".pdf", sep = ""), width = 12, height = 11)
plotBMatrixPosterior(chain$b.Mat, config$aa, interval, param="logmu", main=paste("AA parameter histogram ", prefix, sep=""), nclass=100, center=T)
dev.off()


pdf(paste(result.folder, prefix, "_pMat_trace_", suffix, ".pdf", sep = ""), width = 6, height = 4)
#plotPTraces(data.set$res[[1]]$p.Mat)
plotPTraces(chain$p.Mat)
dev.off()

pdf(paste(result.folder, prefix, "_expPhi_trace_", suffix, ".pdf", sep = ""), width = 6, height = 4)
if(prefix == "w_o" || prefix=="NoPhi")
{
  #plotExpectedPhiTrace(data.set$res[[1]]$phi.pred.Mat)
  plotExpectedPhiTrace(chain$phi.pred.Mat)
}else{
  #plotExpectedPhiTrace(data.set$res[[1]]$phi.Mat)
  plotExpectedPhiTrace(chain$phi.Mat)
}
dev.off()

### Log Likelihood 
pdf(paste(result.folder, prefix, "_logL_trace_", suffix, ".pdf", sep = ""), width = 12, height = 11)
logL <- chain$logL.Mat[[length(chain$logL.Mat)]]
par(mfrow=c(2,2))
for(ii in 0:3){
  sec <- round(length(unlist(chain$logL.Mat))*(ii/4)+1):length(unlist(chain$logL.Mat))
  plot(x = sec-1,
      y = unlist(chain$logL.Mat)[sec]
       #, type="l"
       , xlab="Iteration", ylab="logL"
       #, main=paste("logL:", logL)
       , main=paste(ii+1, "of 4")
       )
  
}
#ll <- plot.likelihood.trace(chain, data, config$use.n.samples)
  dev.off()
}

#### PLOT CUB
{
pdf(paste(result.folder, prefix, "_CUB_obs_bin_", suffix, ".pdf", sep = ""), width = 12, height = 11)
#plotCUB(data$reu13.df, chain$b.Mat, emp, estm.phi[estm.phi[, 1] %in% names(emp), 2],# rescale=T,
plotCUB.NSE(data$reu13.df, bMat=chain$b.Mat, phi.bin=emp,
        model.label="MCMC Posterior", main="CUB binning of observed phi", delta_a12=delta.a_12, a_2=a_2,
        weightedCenters=TRUE, logBins=TRUE)
dev.off()

bin.phis <- estm.phi[estm.phi[, 1] %in% names(emp), 2]
names(bin.phis) <- names(emp)
pdf(paste(result.folder, prefix, "_CUB_est_bin_", suffix, ".pdf", sep = ""), width = 12, height = 11)
#plotCUB(data$reu13.df, chain$b.Mat, bin.phis, estm.phi[estm.phi[, 1] %in% names(emp), 2], 
plotCUB.NSE(data$reu13.df, bMat=chain$b.Mat, phi.bin=bin.phis,
        model.label="MCMC Posterior", main="CUB binning of estimated phi"
        ,delta_a12=delta.a_12, a_2=a_2
        ,weightedCenters=TRUE, logBins=TRUE
        )
dev.off()

pdf(paste(result.folder, prefix, "_vs_obs_phi_", suffix , ".pdf", sep = ""), width = 10, height = 10)
reg <- lm(log10(estm.phi[estm.phi[, 1] %in% names(emp), 2])~log10(emp))
plot(log10(emp), log10(estm.phi[estm.phi[, 1] %in% names(emp), 2]), 
     main=expression(paste("Estimated ", phi, " vs. Observed", phi)),
     sub=paste("R^2=", format(summary(reg)$r.squared, digits = 3), sep=""), 
     xlab=expression(paste("Observed ", phi, " (log10)")), ylab=expression(paste("Est. ", phi, " (log10)")))
abline(reg, col="red", lwd=2)
dev.off()

pdf(paste(result.folder, prefix, "_histogram_", suffix , ".pdf", sep = ""), width = 10, height = 10)
layout(c(1,2))
hist(log10(estm.phi[estm.phi[, 1] %in% names(emp), 2]), nclass=50, main=paste(prefix, " X obs.", sep=""), xlab="estim. phi (log10)",
     sub=paste("mean=", format(mean(log10(estm.phi[estm.phi[, 1] %in% names(emp), 2])), digits = 3), sep=""))
hist(log10(emp), nclass=50, main="excluding unreliable X", xlab="emp. phi (log10)", 
     sub=paste("mean=", format(mean(log10(emp)), digits = 3), sep=""))
dev.off()
}

if(length(convergence)>2)
{
  pdf(paste(result.folder, prefix, "_convergence_trace_", suffix, ".pdf", sep = ""), width = 6, height = 4)
  plot(convergence, xlab="Iteration", ylab="Gelman Score")
  dev.off()
}
