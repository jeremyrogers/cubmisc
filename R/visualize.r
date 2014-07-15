suppressMessages(library(cubfits, quietly = TRUE))
setwd("~/CodonUsageBias/ecoli/R")
source("visualize_utility.r")
source("run_utility.r")
source("config.r")

my.function <- init.function(model = 'roc', model.Phi = 'lognormal',
                                adaptive = 'simple')

result.folder <- "../results/conv_test/chain_1/"
niter <- "multichain_1"
prefix <- "w_o"

#### read data
#data.set <- load.data.file(paste(result.folder, prefix, "_xobs_", niter, ".dat", sep=""))
#estm.phi <- read.csv(paste(result.folder, prefix, "_xobs_", niter, ".phi", sep=""))
#data.set <- load.data.file("../results/rm_first_50_aa/without_xobs_multichain.dat")
load(paste(result.folder, "without_xobs_multichain_lapply_4cores.dat", sep=""))
estm.phi <- read.csv(paste(result.folder, "without_xobs_multichain_lapply_4cores.phi", sep=""))


seq.string <- readGenome("../data/ecoli_K12_MG1655_genome_filtered.fasta", config$rm.short, config$rm.first.aa)
emp <- read.empirical.data("../data/ecoli_X_obs.csv", seq.string, config$selected.env, th=0)
#seq.string <- emp$genome
emp <- emp$empirical[names(emp$empirical) %in% estm.phi[, 1]]


#### PLOT TRACES
pdf(paste(result.folder, prefix, "_logmu_", niter, ".pdf", sep = ""), width = 12, height = 11)
#plotTraces(data.set$res[[1]]$b.Mat, config$aa, param="logmu", main=paste("AA parameter trace ", prefix, sep=""))
plotTraces(chain$b.Mat, config$aa, param="logmu", main=paste("AA parameter trace ", prefix, sep=""))
dev.off()
 
pdf(paste(result.folder, prefix, "_deltat_", niter, ".pdf", sep = ""), width = 12, height = 11)
#plotTraces(data.set$res[[1]]$b.Mat, config$aa, param="deltat", main=paste("AA parameter trace ", prefix, sep=""))
plotTraces(chain$b.Mat, config$aa, param="deltat", main=paste("AA parameter trace ", prefix, sep=""))
dev.off()

pdf(paste(result.folder, prefix, "_pMat_trace_", niter, ".pdf", sep = ""), width = 6, height = 4)
#plotPTraces(data.set$res[[1]]$p.Mat)
plotPTraces(chain$p.Mat)
dev.off()

pdf(paste(result.folder, prefix, "_expPhi_trace_", niter, ".pdf", sep = ""), width = 6, height = 4)
if(prefix == "w_o")
{
  #plotExpectedPhiTrace(data.set$res[[1]]$phi.pred.Mat)
  plotExpectedPhiTrace(chain$phi.pred.Mat)
}else{
  #plotExpectedPhiTrace(data.set$res[[1]]$phi.Mat)
  plotExpectedPhiTrace(chain$phi.Mat)
}
dev.off()



ll <- plot.likelihood.trace(data.set)
mean(ll)



#### PLOT CUB
pdf(paste(result.folder, prefix, "_CUB_obs_bin_", niter, ".pdf", sep = ""), width = 12, height = 11)
# plotCUB(data.set$data$reu13.df, data.set$res[[1]]$b.Mat, emp, estm.phi[estm.phi[, 1] %in% names(emp), 2], rescale=T,
#          model.label="MCMC Posterior", main="CUB binning of observed phi")
plotCUB(data$reu13.df, chain$b.Mat, emp, estm.phi[estm.phi[, 1] %in% names(emp), 2], rescale=T,
        model.label="MCMC Posterior", main="CUB binning of observed phi")
dev.off()

bin.phis <- estm.phi[estm.phi[, 1] %in% names(emp), 2]
names(bin.phis) <- names(emp)
pdf(paste(result.folder, prefix, "_CUB_est_bin_", niter, ".pdf", sep = ""), width = 12, height = 11)
# plotCUB(data.set$data$reu13.df, data.set$res[[1]]$b.Mat, bin.phis, estm.phi[estm.phi[, 1] %in% names(emp), 2], 
#          model.label="MCMC Posterior", main="CUB binning of estimated phi")
plotCUB(data$reu13.df, chain$b.Mat, bin.phis, estm.phi[estm.phi[, 1] %in% names(emp), 2], 
        model.label="MCMC Posterior", main="CUB binning of estimated phi")
dev.off()





### COMPARE TWO RUNS
w.phi <- read.csv(paste(result.folder, "wphi_no_low_conf_80k.txt", sep=""))
wo.phi <- read.csv(paste(result.folder, "wophi_no_low_conf_80k.txt", sep=""))

pdf(paste(result.folder, "withphi_vs_fitIC_sdlog05", ".pdf", sep = ""), width = 10, height = 10)
xy.max <- max(cbind(w.phi[,2], wo.phi[,2]))
reg <- lm(wo.phi[, 2]~w.phi[, 2])
plot(w.phi[, 2], wo.phi[,2], xlim=c(0, xy.max), ylim=c(0,xy.max),
     main=expression(paste("with ", phi, " fit vs. fitted IC")),
     sub=paste("R^2=", format(summary(reg)$r.squared, digits = 3), sep=""),
     xlab=expression(paste("with obs. ", phi)), ylab=expression(paste("fit. I.C. ", phi)))
abline(0,1, lty=2)
abline(reg, col="red", lwd=2)
dev.off()




pdf(paste(result.folder, prefix, "_vs_obs_phi_", niter , ".pdf", sep = ""), width = 10, height = 10)
reg <- lm(log10(estm.phi[estm.phi[, 1] %in% names(emp), 2])~log10(emp))
plot(log10(emp), log10(estm.phi[estm.phi[, 1] %in% names(emp), 2]), 
     main=expression(paste("Estim. ", phi, " vs. Obs. ", phi)),
     sub=paste("R^2=", format(summary(reg)$r.squared, digits = 3), sep=""), 
     xlab=expression(paste("Obs. ", phi, " (log10)")), ylab=expression(paste("Est. ", phi, " (log10)")))
abline(reg, col="red", lwd=2)
dev.off()


pdf(paste(result.folder, prefix, "_histogram_", niter , ".pdf", sep = ""), width = 10, height = 10)
layout(c(1,2))
hist(log10(estm.phi[estm.phi[, 1] %in% names(emp), 2]), nclass=50, main=paste(prefix, " X obs.", sep=""), xlab="estim. phi (log10)",
     sub=paste("mean=", format(mean(log10(estm.phi[estm.phi[, 1] %in% names(emp), 2])), digits = 3), sep=""))
hist(log10(emp), nclass=50, main="excluding unreliable X", xlab="emp. phi (log10)", 
     sub=paste("mean=", format(mean(log10(emp)), digits = 3), sep=""))
dev.off()




plot(convergence, xlab="Iteration", ylab="Gelman Score")
data <- as.data.frame(convergence)
data <- list(x=data[,1], y=data[,2])
model <- nls(y ~ 1.213+exp(c + d * x), data = data, start = list(c = 0, d = 0))


iter <- seq(1, max(convergence[,1]))
pred.score <- (predict(model, list(x=iter)))
plot(data$x, data$y, xlab="Iteration", ylab="Gelman Score")
lines(iter, pred.score, type="l")


plot.w.phi.vs.wo.phi(res.w.phi, res.wo.phi))
plot.pred.vs.emp(data.set$mean.phis)



estim.phi <- read.csv(paste(result.folder, "sdlog05_fittedIC.txt", sep=""))
pdf(paste(result.folder, "metrics_sdlog05_fittedIC.pdf", sep=""), width=16, height=10)
plot(estim.phi[, 2:5], main=expression(paste("est. ", bar(phi))))
dev.off()







#ecoli_sdlog1_vs_mRNA.png

## 3D plot all empirical data in 3D
#library(scatterplot3d)
## delete genes with 0 synthesis rate
#emp <- empirical[empirical[,1] != 0, 1:3]
#emp <- emp[emp[,2] != 0, 1:3]
#emp <- emp[emp[,3] != 0, 1:3]

#reg2 <- lm(log10(emp[,3])~log10(emp[,1])+log10(emp[,2]))
#scatterplot3d(log10(emp), highlight.3d=T, main="log10(empirical phi)", 
#              sub=paste("R^2=", format(summary(reg2)$r.squared, digits = 3), sep="") )
#reg2 <- lm(log10(emp[,2])~log10(emp[,1]))
#plot(log10(emp[,1]), log10(emp[,2]))



mean(c(-1325842.757, 
  -1325578.915, 
  -1325357.836, 
  -1325792.283, 
  -1325699.128, 
  -1325463.118, 
  -1325262.055, 
  -1325442.114))