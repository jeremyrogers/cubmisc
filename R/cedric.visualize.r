rm(list=ls())
suppressMessages(library(cubfits, quietly = TRUE))
#setwd("~/CodonUsageBias/R_roc")
source("visualize_utility.r")
source("run_utility.r")
source("config.r")

my.function <- init.function(model = 'roc', model.Phi = 'lognormal',
                                adaptive = 'simple')

has.emp.phi <- F

#organism <- ""
chainnumber <- 4
#result.folder <- paste("../", organism, "/results/withoutXobs_noprior_purpureum_5k_randgenesamples/", sep="")
result.folder <- "../results/brewersyeast/"
#result.folder <- paste("../organisms/", organism, "/results/CodonParamEvolution/Skluyveri_normPrior_preDE/chain_", chainnumber, "/", sep="")
#suffix <- paste("multichain_", chainnumber, sep="")
suffix <- "singlechain"
#result.folder <- paste("../organisms/", organism, "/results/Cbescii_excludeShortGenes150aa_18k/", sep="")
prefix <- "w_o"

#### read data
#data.set <- load.data.file(paste(result.folder, prefix, "_xobs_", suffix, ".dat", sep=""))
#estm.phi <- read.csv(paste(result.folder, prefix, "_xobs_", suffix, ".phi", sep=""))
#data.set <- load.data.file("../results/rm_first_50_aa/without_xobs_multichain.dat")
#load(paste(result.folder, "without_xobs_singlechain.dat", sep=""))
load("../results/brewersyeast/brewersYeastwithPhi.dat")
estm.phi <- read.csv(paste(result.folder, "brewersYeastwithPhi.phi", sep=""), as.is=T)

estm.phi <- data.frame(ORF_Info = as.character(estm.phi[, 1]), Phi_Post_Arith_Mean=as.numeric(estm.phi[, 2]), Phi_Post_Median=as.numeric(estm.phi[, 3]), 
           Phi_Post_Geom_Mean=as.numeric(estm.phi[, 5]), Phi_Post_Harm_Mean=as.numeric(estm.phi[, 5]), Phi_Post_SD=as.numeric(estm.phi[, 6]))


if(has.emp.phi)
{
  seq.string.names <- estm.phi[,1]
  names(seq.string.names) <- seq.string.names
  emp <- read.empirical.data("../brewersYeast/yassour2009.phi.tsv", seq.string.names, config$selected.env, th=0)
  emp <- emp$empirical[names(emp$empirical) %in% estm.phi[, 1]]
  emp <- emp[order(names(emp))]
}


#### PLOT TRACES
pdf(paste(result.folder, prefix, "_logmu_", suffix, ".pdf", sep = ""), width = 11, height = 12)
#plotTraces(data.set$res[[1]]$b.Mat, config$aa, param="logmu", main=paste("AA parameter trace ", prefix, sep=""))
plotTraces(chain$b.Mat, config$aa, param="logmu", main=paste("AA parameter trace ", prefix, sep=""))
dev.off()
 
pdf(paste(result.folder, prefix, "_deltaeta_", suffix, ".pdf", sep = ""), width = 11, height = 12)
#plotTraces(data.set$res[[1]]$b.Mat, config$aa, param="deltat", main=paste("AA parameter trace ", prefix, sep=""))
plotTraces(chain$b.Mat, config$aa, param="deltaeta", main=paste("AA parameter trace ", prefix, sep=""))
dev.off()

interval <- (length(chain$b.Mat)-config$use.n.samples):length(chain$b.Mat)
pdf(paste(result.folder, prefix, "_deltaeta_histogram_", suffix, ".pdf", sep = ""), width = 11, height = 12)
plotBMatrixPosterior(chain$b.Mat, config$aa, interval, param="deltaeta", main=paste("AA parameter histogram", prefix, sep=""), nclass=100, center=F)
dev.off()

pdf(paste(result.folder, prefix, "_logmu_histogram_", suffix, ".pdf", sep = ""), width = 11, height = 12)
plotBMatrixPosterior(chain$b.Mat, config$aa, interval, param="logmu", main=paste("AA parameter histogram ", prefix, sep=""), nclass=100, center=F)
dev.off()

pdf(paste(result.folder, prefix, "_pMat_trace_", suffix, ".pdf", sep = ""), width = 6, height = 4)
if(prefix == "w_o"){
  ylabs <- c(expression("M"[phi]), expression("log("~sigma[phi]~")"))
}else{
  ylabs <- c(expression("log("~sigma[epsilon]~")"), expression("M"[phi]), expression("log("~sigma[phi]~")"), expression("A"[phi]))
}
plotPTraces(chain$p.Mat, log = "y", ylab = ylabs)
dev.off()

pdf(paste(result.folder, prefix, "_expPhi_trace_", suffix, ".pdf", sep = ""), width = 6, height = 4)
if(prefix == "w_o")
{
  #plotExpectedPhiTrace(data.set$res[[1]]$phi.pred.Mat)
  plotExpectedPhiTrace(chain$phi.pred.Mat)
}else{
  #plotExpectedPhiTrace(data.set$res[[1]]$phi.Mat)
  plotExpectedPhiTrace(chain$phi.Mat)
}
dev.off()

pdf(paste(result.folder, prefix, "_logL_trace_", suffix, ".pdf", sep = ""), width = 8, height = 6)
chainLength <- length(chain$logL.Mat)
zoomStart <- round(0.9*length(chain$logL.Mat))
ll <- unlist(chain$logL.Mat)
logL <- mean(ll[zoomStart:chainLength])
plot(ll, type="l", main=paste("logL:", logL), xlab="Sample", ylab="logL")
grid (NULL,NULL, lty = 6, col = "cornsilk2") 
Hmisc::subplot(
  plot(zoomStart:chainLength, ll[zoomStart:chainLength], type="l", xlab=NA, ylab=NA, las=2, cex.axis=0.55), 
  0.8*zoomStart, (min(ll)+max(ll))/2, size=c(3,2))
#ll <- plot.likelihood.trace(chain, data, config$use.n.samples)
dev.off()


#### PLOT CUB
if(has.emp.phi)
{
  scaled.emp <- 10^( log10(emp) - log10(mean(emp)) )
  pdf(paste(result.folder, prefix, "_CUB_obs_bin_", suffix, ".pdf", sep = ""), width = 11, height = 12)
  # plotCUB(data.set$data$reu13.df, data.set$res[[1]]$b.Mat, emp, estm.phi[estm.phi[, 1] %in% names(emp), 2], rescale=T,
  #          model.label="MCMC Posterior", main="CUB binning of observed phi")
  plotCUB(data$reu13.df, bMat=chain$b.Mat, phi.bin=scaled.emp,
          model.label="MCMC Posterior", main="CUB binning of observed phi")
  dev.off()
}

bin.phis <- estm.phi[, 2]
names(bin.phis) <- unique(data$reu13.df$A$ORF) estm.phi[, 1]
pdf(paste(result.folder, prefix, "_CUB_est_bin_", suffix, ".pdf", sep = ""), width = 11, height = 12)
# plotCUB(data.set$data$reu13.df, data.set$res[[1]]$b.Mat, bin.phis, estm.phi[estm.phi[, 1] %in% names(emp), 2], 
#          model.label="MCMC Posterior", main="CUB binning of estimated phi")
plotCUB(data$reu13.df, bMat=chain$b.Mat, phi.bin=bin.phis, 
        model.label="MCMC Posterior", main="CUB binning of estimated phi")
dev.off()

if(has.emp.phi){
  pdf(paste(result.folder, prefix, "_vs_obs_phi_", suffix , ".pdf", sep = ""), width = 10, height = 10)
  reg <- lm(log10(estm.phi[estm.phi[, 1] %in% names(emp), 2])~log10(emp))
  plot(log10(emp), log10(estm.phi[estm.phi[, 1] %in% names(emp), 2]), 
       main=expression(paste("Estim. ", phi, " vs. Obs. ", phi)),
       sub=paste("R^2=", format(summary(reg)$r.squared, digits = 3), sep=""), 
       xlab=expression(paste("Obs. ", phi, " (log10)")), ylab=expression(paste("Est. ", phi, " (log10)")))
  grid (NULL,NULL, lty = 6, col = "cornsilk2") 
  #points(log10(emp[idxemp]), log10(phi.v), col="red")
  abline(reg, col="red", lwd=2)
  dev.off()
}
# reg <- lm(log10(scuo)~log10(phi.obs))
# plot(log10(phi.obs), log10(scuo), 
#      xlab=expression(paste("Emp. ", phi, " (log10)")), ylab="SCUO (log10)", 
#      main=expression(paste("SCUO vs. Emp. ", phi, " (GSM552568)")),
#      sub=paste("R^2=", format(summary(reg)$r.squared, digits = 3), sep=""))
# abline(reg, col="red", lwd=2)

pdf(paste(result.folder, prefix, "_histogram_", suffix , ".pdf", sep = ""), width = 10, height = 10)
if(has.emp.phi){
  layout(c(1,2))
  orfs <- estm.phi[, 1] %in% names(emp)
}else{
  orfs <- estm.phi[, 1]
}
hist(log10(estm.phi[orfs, 2]), nclass=100, main=paste("estm phi ", prefix, " X obs.", sep=""), xlab="estim. phi (log10)",
     sub=paste("mean=", format(mean(log10(estm.phi[orfs, 2])), digits = 3), sep=""))
if(has.emp.phi){
  hist(log10(emp), nclass=100, main="Obs. X", xlab="emp. phi (log10)", 
       sub=paste("mean=", format(mean(log10(emp)), digits = 3), sep=""))
}
dev.off()

pdf(paste(result.folder, prefix, "_convergence_trace_", suffix, ".pdf", sep = ""), width = 6, height = 4)
plot(convergence, xlab="Iteration", ylab="Gelman Score")
grid (NULL,NULL, lty = 6, col = "cornsilk2") 
dev.off()

phis <- do.call("cbind", chain$phi.pred.Mat)
pdf(paste(result.folder, prefix, "_phi_trace_", suffix, ".pdf", sep = ""), width = 4, height = 4)
candidates <- which(log10(estm.phi[, 2]) < -0.25)
samples <- sample(x = candidates, size = 25)
for(idx in candidates){
  plot(log10(phis[idx,2000:10000]), type="l", main=names(phis[idx,1]), xlab="Sample", ylab=expression("log10("~phi~")"))
}
dev.off()
