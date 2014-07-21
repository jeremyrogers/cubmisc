suppressMessages(library(cubfits, quietly = TRUE))
##setwd("~/CodonUsageBias/ecoli/R")
source("visualize_utility.r")
source("run_utility.r")
source("config.r")

##Local variable definitions
##Genome to be analyzed
genome.name = "EE36"
fasta.name.suffix = "_genome_curated.fasta"
genome.folder.location="../../"
genome.file = paste(genome.folder.location, genome.name, fasta.name.suffix, sep="")
run.suffix = "" #any suffix in results folder name
chain.id = 4;
results.folder <- paste("../../results/", genome.name, run.suffix, "/chain_", chain.id, "/", sep="")
plot.label <- paste("multichain", chain.id, sep="") ##analysis label
read.xobs.data=FALSE;
model.fit = "without";

my.function <- init.function(model = 'roc', model.Phi = 'lognormal',
                                adaptive = 'simple')

##Redefine plotPTraces so you can request logrithmic axes
plotPTraces <- function (ret.p.Mat, main = "Hyperparameter Traces", lty = 1, ...) 
{
    ret.p.Mat <- do.call("rbind", ret.p.Mat)
    n.traces <- dim(ret.p.Mat)[2]
    if (n.traces == 2) {
        ylabs <- c("M", expression("s"[phi]))
    }
    else {
        ylabs <- c(expression("s"[epsilon]), "M", expression("s"[phi]), 
            "K")
    }
    par(oma = c(1, 1, 2, 1), mgp = c(2, 1, 0), mar = c(3, 4, 
        2, 1), mfrow = c(2, ceiling(n.traces/2)))
    for (i in 1:n.traces) {
        plot(ret.p.Mat[, i], xlab = "Iteration", ylab = ylabs[i], 
            type = "l", lty = 1, ...)
    }
    title(main = main, outer = T)
}


#### read data
#data.set <- load.data.file(paste(results.folder, model.fit, "_xobs_", plot.label, ".dat", sep=""))
#estm.phi <- read.csv(paste(results.folder, model.fit, "_xobs_", plot.label, ".phi", sep=""))
#data.set <- load.data.file("../results/rm_first_50_aa/without_xobs_multichain.dat")
data.set <- load(paste(results.folder, "without_xobs_multichain_lapply_4cores.dat", sep=""))
estm.phi <- read.csv(paste(results.folder, "without_xobs_multichain_lapply_4cores.phi", sep=""))


seq.string <- readGenome(genome.file, config$rm.short, config$rm.first.aa)

if(read.xobs.data){
  emp <- read.empirical.data("../data/ecoli_X_obs.csv", seq.string, config$selected.env, th=0) #should be called 'read.Xobs.data'
  ##seq.string <- emp$genome
  emp <- emp$empirical[names(emp$empirical) %in% estm.phi[, 1]]
}

#### PLOT TRACES
pdf(paste(results.folder, model.fit, "_logmu_", plot.label, ".pdf", sep = ""), width = 12, height = 11)
#plotTraces(data.set$res[[1]]$b.Mat, config$aa, param="logmu", main=paste("AA parameter trace ", model.fit, sep=""))
plotTraces(chain$b.Mat, config$aa, param="logmu", main=paste("AA parameter trace ", model.fit, sep=""))
dev.off()
 
pdf(paste(results.folder, model.fit, "_deltat_", plot.label, ".pdf", sep = ""), width = 12, height = 11)
#plotTraces(data.set$res[[1]]$b.Mat, config$aa, param="deltat", main=paste("AA parameter trace ", model.fit, sep=""))
plotTraces(chain$b.Mat, config$aa, param="deltat", main=paste("AA parameter trace ", model.fit, sep=""))
dev.off()

pdf(paste(results.folder, model.fit, "_pMat_trace_", plot.label, ".pdf", sep = ""), width = 6, height = 4)
#plotPTraces(data.set$res[[1]]$p.Mat)
plotPTraces(chain$p.Mat, log="y")
dev.off()

pdf(paste(results.folder, model.fit, "_expPhi_trace_", plot.label, ".pdf", sep = ""), width = 6, height = 4)
if(model.fit == "without")
{
  #plotExpectedPhiTrace(data.set$res[[1]]$phi.pred.Mat)
  plotExpectedPhiTrace(chain$phi.pred.Mat)
}else{
  #plotExpectedPhiTrace(data.set$res[[1]]$phi.Mat)
  plotExpectedPhiTrace(chain$phi.Mat)
}
dev.off()


#### PLOT LLIK TRACE
##THIS DOESN"T WORK!! GET: Error in par(oma = c(1, 1, 2, 1), mgp = c(2, 1, 0), mar = c(3, 4, 2, 1),  : 
##  invalid value specified for graphical parameter "mfrow"
##
## This seems to be because input and, therefore, input$res are not defined
ll <- plot.likelihood.trace(data.set)
mean(ll)



#### PLOT CUB
## Plot vs. Xobs data
pdf(paste(results.folder, model.fit, "_CUB_obs_bin_", plot.label, ".pdf", sep = ""), width = 12, height = 11)
# plotCUB(data.set$data$reu13.df, data.set$res[[1]]$b.Mat, emp, estm.phi[estm.phi[, 1] %in% names(emp), 2], rescale=T,
#          model.label="MCMC Posterior", main="CUB binning of observed phi")
plotCUB(data$reu13.df, chain$b.Mat, emp, estm.phi[estm.phi[, 1] %in% names(emp), 2], rescale=T,
        model.label="MCMC Posterior", main="CUB binning of observed phi")
dev.off()

##Plot vs. estimated phi values 
bin.phis <- estm.phi[ ,2] #[estm.phi[, 1] %in% names(emp), 2]
names(bin.phis) <- estm.phi[,1]
pdf(paste(results.folder, model.fit, "_CUB_est_bin_", plot.label, ".pdf", sep = ""), width = 12, height = 11)
# plotCUB(data.set$data$reu13.df, data.set$res[[1]]$b.Mat, bin.phis, estm.phi[estm.phi[, 1] %in% names(emp), 2], 
#          model.label="MCMC Posterior", main="CUB binning of estimated phi")
plotCUB(data$reu13.df, chain$b.Mat, bin.phis, bin.phis, 
        model.label="MCMC Posterior", main="CUB binning of estimated phi")
dev.off()



### COMPARE TWO RUNS
run.one <- read.csv(paste(results.folder, "wphi_no_low_conf_80k.txt", sep=""))
run.two <- read.csv(paste(results.folder, "wophi_no_low_conf_80k.txt", sep=""))

pdf(paste(results.folder, "withphi_vs_fitIC_sdlog05", ".pdf", sep = ""), width = 10, height = 10)
xy.max <- max(cbind(run.one[,2], run.two[,2]))
reg <- lm(run.two[, 2]~run.one[, 2])
plot(run.one[, 2], run.two[,2], xlim=c(0, xy.max), ylim=c(0,xy.max),
     main=expression(paste("with ", phi, " fit vs. fitted IC")),
     sub=paste("R^2=", format(summary(reg)$r.squared, digits = 3), sep=""),
     xlab=expression(paste("with obs. ", phi)), ylab=expression(paste("fit. I.C. ", phi)))
abline(0,1, lty=2)
abline(reg, col="red", lwd=2)
dev.off()




pdf(paste(results.folder, model.fit, "_vs_obs_phi_", plot.label , ".pdf", sep = ""), width = 10, height = 10)
reg <- lm(log10(estm.phi[estm.phi[, 1] %in% names(emp), 2])~log10(emp))
plot(log10(emp), log10(estm.phi[estm.phi[, 1] %in% names(emp), 2]), 
     main=expression(paste("Estim. ", phi, " vs. Obs. ", phi)),
     sub=paste("R^2=", format(summary(reg)$r.squared, digits = 3), sep=""), 
     xlab=expression(paste("Obs. ", phi, " (log10)")), ylab=expression(paste("Est. ", phi, " (log10)")))
abline(reg, col="red", lwd=2)
dev.off()


pdf(paste(results.folder, model.fit, "_histogram_", plot.label , ".pdf", sep = ""), width = 10, height = 10)
layout(c(1,2))
hist(log10(estm.phi[estm.phi[, 1] %in% names(emp), 2]), nclass=50, main=paste(model.fit, " X obs.", sep=""), xlab="estim. phi (log10)",
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


plot.w.phi.vs.wo.phi(res.run.one, res.run.two))
plot.pred.vs.emp(data.set$mean.phis)



estim.phi <- read.csv(paste(results.folder, "sdlog05_fittedIC.txt", sep=""))
pdf(paste(results.folder, "metrics_sdlog05_fittedIC.pdf", sep=""), width=16, height=10)
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
