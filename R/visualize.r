library(cubfits, quietly = TRUE, lib.loc = "~/cubfitsNSEdebug/")
setwd("~/cubfits/misc/R/")
source("visualize_utility.r")
source("run_utility.r")
source("config.r")

model <- 'nse';              #nse \ roc

genome <- 'brewYeast';       #ecoli \ REUyeast \ pYeast \ brewYeast
prefix <- "11-17"            #generally, date the run started in "MM-DD" format
suffix <- "FULL"             #What was special about this run?

delta_a12 <- 0
a_2 <- 4
  

if(tolower(genome)=="ecoli"){
  simulated_data <- FALSE;

#### read data
  result.folder <- paste("../results/test/", substr(tolower(model),1,1), substr(tolower(genome),1,1), "/", sep="")
  load(paste(result.folder, paste("debug", model, genome, ".dat", sep="")))
  estm.phi <- read.csv(paste(result.folder, genome, ".", model, ".", suffix, ".phi" , sep=""))
  seq.string <- readGenome("../data/ecoli_K12_MG1655_genome_filtered.fasta", config$rm.short, config$rm.first.aa)
  emp <- read.empirical.data("../data/ecoli_X_obs.csv", seq.string, config$selected.env, th=0)
} else
if(tolower(genome)=="reuyeast"){
  simulated_data <- TRUE;

#### read data
  result.folder <- paste("../results/ny/", prefix, "/", sep="")
  load(paste(result.folder, "REU.nse.", suffix, ".dat", sep=""))
  estm.phi <- read.csv(paste(result.folder, "REU.nse.", suffix, ".phi", sep=""))
  seq.string <- readGenome("../S.cervisiae.REU13/section1.fasta", config$rm.short, config$rm.first.aa)
  emp <- read.empirical.data("../S.cervisiae.REU13/section1sorted.csv", seq.string, config$selected.env, th=0)


  omega_true <- read.table("../S.cervisiae.REU13/scaled_omega.csv", header=TRUE)
  mu_true <- read.table("../S.cervisiae.REU13/scaled_logMu.csv", header=TRUE)
  names(chain$b.Mat[[1]]) <- read.table("names.txt")[[1]]
    ###this needs to be fixed!!!
} else
if(tolower(genome)=="pyeast"){
  simulated_data <- TRUE;
  
#### read data
  result.folder <- paste("../results/", substr(tolower(model),1,1), substr(tolower(genome),1,1), "/", prefix, "/", sep="")
  load(paste(result.folder, genome, ".", model, ".", suffix, ".dat" , sep=""))
  estm.phi <- read.csv(paste(result.folder, genome, ".", model, ".", suffix, ".phi" , sep=""))

#  seq.string <- readGenome("../prestonYeast/section1of5.fasta", config$rm.short, config$rm.first.aa)
#  emp <- read.empirical.data("../prestonYeast/section1of5_sorted.csv", seq.string, config$selected.env, th=0)
  seq.string <- readGenome("../prestonYeast/S.cerevisiae.S288c.fasta", config$rm.short, config$rm.first.aa)
  emp <- read.empirical.data("../prestonYeast/pyeast.phi.tsv", seq.string, config$selected.env, th=0)
  omega_true <- read.table("../prestonYeast/scaled_omega.csv", header = TRUE)
  mu_true <- read.table("../prestonYeast/scaled_logMu.csv", header = TRUE)
  names(chain$b.Mat[[1]]) <- read.table("names.txt")[[1]]
    ###this needs to be fixed!!!
} else
if(tolower(genome)=="brewyeast"){
  simulated_data <- FALSE;
  
  #### read data
  result.folder <- paste("../results/", substr(tolower(model),1,1), substr(tolower(genome),1,1), "/", prefix, "/", sep="")
  load(paste(result.folder, genome, ".", model, ".", suffix, ".dat" , sep=""))
  estm.phi <- read.csv(paste(result.folder, genome, ".", model, ".", suffix, ".phi" , sep=""))
  
  seq.string <- readGenome("../brewersYeast/s288c.genome.fasta", config$rm.short, config$rm.first.aa)
  emp <- read.empirical.data("../brewersYeast/yassour2009.phi.tsv", seq.string, config$selected.env, th=0)
#  omega_true <- read.table("../brewYeast/scaled_omega.csv", header = TRUE)
#  mu_true <- read.table("../brewYeast/scaled_logMu.csv", header = TRUE)
  names(chain$b.Mat[[1]]) <- read.table("names.txt")[[1]]
  ###this needs to be fixed!!!
}


emp <- emp$empirical[names(emp$empirical) %in% estm.phi[, 1]]

if(simulated_data){
  require(package = plotrix, lib.loc = "/home/lbrown/cubfits/Dependencies/")
  chunksize <- floor(length(chain$b.Mat)/100)*10
  ci <- .95 #Confidence interval, on a scale of (0,1). e.g. 0.95 means 95%CI
  
  
  #average converged omega
{
  omega <- matrix(0, nrow=length(omega_true[[3]]), ncol=4)
  colnames(omega) <- c("'true' deltaomega", "estimated deltaomega"
                       ,paste("Lower ", ci, "%ci", sep="") , paste("Upper ", ci, "%ci", sep="")
  )
  omega[,1] <- omega_true[[3]];
  vec <- grep("eta", names(chain$b.Mat[[1]]))
  tmpomega <- matrix(ncol=chunksize, nrow=length(vec));
  
  for(ii in 1:length(vec)){
    for(jj in 1:chunksize){
      tmpomega[ii,jj] <- chain$b.Mat[[jj+length(chain$b.Mat)-chunksize]][vec[ii]]
    }##finish current codon
    omega[ii,2] <- mean(tmpomega[ii,])
    omega[ii,3:4] <- quantile(tmpomega[ii,], probs = c( (1-ci)/2 , (1+ci)/2))
    #omega[ii,3] <- omega[ii,2] - qnorm((1+ci)/2) * sd(tmpomega[ii,])/sqrt(chunksize)
    #omega[ii,4] <- omega[ii,2] + qnorm((1+ci)/2) * sd(tmpomega[ii,])/sqrt(chunksize)
  }##finish all codons
  
}
#omega[,1] is the 'true' omega and omega[,2] is the estimated omega
#omega[,3] and [,4] are the lower and upper 95% confidence interval

#average converged mu
{
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
    #mu[ii,3] <- mu[ii,2] - qnorm((1+ci)/2) * sd(tmpmu[ii,])/sqrt(chunksize)
    #mu[ii,4] <- mu[ii,2] + qnorm((1+ci)/2) * sd(tmpmu[ii,])/sqrt(chunksize)
  }##finish all codons
  
}

#### PLOT RESULTS VS 'TRUE' VALUES
# ?
#Plot Omega
pdf(paste(result.folder, prefix, "_codonParameters_", suffix, ".pdf", sep = ""), width = 8, height = 10)
par(mfrow=c(2,1))
reg <- lm(omega[,2]~omega[,1]);
plot(omega, main=expression(paste("'true' ",Delta,omega," vs Estimated ",Delta,omega, "without problem entries"))
     ,ylim=range(omega[,3:4])
     ,xlab=expression(paste("'true' ", Delta, omega))
     ,ylab=expression(paste("estimated ", Delta, omega))
     ,sub=paste("R^2=", format(summary(reg)$r.squared, digits = 3), sep="") 
)
abline(reg, lwd=2, col="red")
abline(0,1, lwd=2, col="blue")
plotCI(x = omega[,1], y=omega[,2], ui=omega[,4], li=omega[,3], add=TRUE, scol = "darkslategrey");
legend("bottomright", legend = c("1-1 line", "Linear Fit Line", paste(100*ci, "% confidence", sep=""))
       ,fill=c("blue", "red", "darkslategrey")
)

###Plot Mu
reg <- lm(mu[,2]~mu[,1]);
plot(mu, main=expression(paste("'true' ",Delta,"M vs Estimated ",Delta,"M"))
     ,ylim=range(mu[,3:4])
     ,xlab=expression(paste("'true' ", Delta, "M"))
     ,ylab=expression(paste("estimated ", Delta,"M"))
     ,sub=paste("R^2=", format(summary(reg)$r.squared, digits = 3), sep="")
)
abline(0,1, lwd=2, col="blue")
abline(reg, lwd=2, col="red")
plotCI(x = mu[,1], y=mu[,2], ui=mu[,4], li=mu[,3], add=TRUE, scol = "darkslategrey");
legend("bottomright", legend = c("1-1 line", "Linear Fit Line", paste(100*ci, "% confidence", sep=""))
       ,fill=c("blue", "red", "darkslategrey")
)

dev.off()
}

#### PLOT TRACES
{
pdf(paste(result.folder, prefix, "_logmu_", suffix, ".pdf", sep = ""), width = 12, height = 11)
#plotTraces(data.set$res[[1]]$b.Mat, config$aa, param="logmu", main=paste("AA parameter trace ", prefix, sep=""))
plotTraces(chain$b.Mat, config$aa, param="logmu", main=paste("AA parameter trace ", prefix, sep=""))
dev.off()
 
pdf(paste(result.folder, prefix, "_deltaeta_", suffix, ".pdf", sep = ""), width = 12, height = 11)
#plotTraces(data.set$res[[1]]$b.Mat, config$aa, param="deltaeta", main=paste("AA parameter trace ", prefix, sep=""))
plotTraces(chain$b.Mat, config$aa, param="deltaeta", main=paste("AA parameter trace ", prefix, sep=""))
dev.off()

pdf(paste(result.folder, prefix, "_pMat_trace_", suffix, ".pdf", sep = ""), width = 6, height = 4)
#plotPTraces(data.set$res[[1]]$p.Mat)
plotPTraces(chain$p.Mat)
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

### Log Likelihood 
  pdf(paste(result.folder, prefix, "_logL_trace_", suffix, ".pdf", sep = ""), width = 6, height = 4)
  ll <- plot.likelihood.trace(chain, data, config$use.n.samples)
  dev.off()
}

#### PLOT CUB
{
pdf(paste(result.folder, prefix, "_CUB_obs_bin_", suffix, ".pdf", sep = ""), width = 12, height = 11)
plotCUB(data$reu13.df, chain$b.Mat, emp, estm.phi[estm.phi[, 1] %in% names(emp), 2], rescale=T,
        model.label="MCMC Posterior", main="CUB binning of observed phi"
        #,model='nse'
        ,delta_a12=delta_a12, a_2=a_2
        )
dev.off()

bin.phis <- estm.phi[estm.phi[, 1] %in% names(emp), 2]
names(bin.phis) <- names(emp)
pdf(paste(result.folder, prefix, "_CUB_est_bin_", suffix, ".pdf", sep = ""), width = 12, height = 11)
plotCUB(data$reu13.df, chain$b.Mat, bin.phis, estm.phi[estm.phi[, 1] %in% names(emp), 2], 
        model.label="MCMC Posterior", main="CUB binning of estimated phi"
        #,model=model
        ,delta_a12=delta_a12, a_2=a_2
        )
dev.off()

pdf(paste(result.folder, prefix, "_vs_obs_phi_", suffix , ".pdf", sep = ""), width = 10, height = 10)
reg <- lm(log10(estm.phi[estm.phi[, 1] %in% names(emp), 2])~log10(emp))
plot(log10(emp), log10(estm.phi[estm.phi[, 1] %in% names(emp), 2]), 
     main=expression(paste("Estim. ", phi, " vs. 'true'", phi)),
     sub=paste("R^2=", format(summary(reg)$r.squared, digits = 3), sep=""), 
     xlab=expression(paste("'true' ", phi, " (log10)")), ylab=expression(paste("Est. ", phi, " (log10)")))
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

pdf(paste(result.folder, prefix, "_convergence_trace_", suffix, ".pdf", sep = ""), width = 6, height = 4)
plot(convergence, xlab="Iteration", ylab="Gelman Score")
dev.off()

# data <- as.data.frame(convergence)
# data <- list(x=data[,1], y=data[,2])
# model <- nls(y ~ 1.2+exp(c + d * x), data = data, start = list(c = 0, d = 0))
# 
# 
# iter <- seq(1, 2*max(convergence[,1]))
# pred.score <- (predict(model, list(x=iter)))
# plot(iter, pred.score, type="l", xlab="Iteration", ylab="Gelman Score", main="Convergence trend")
# points(data$x, data$y)
# 
# 
# plot.w.phi.vs.wo.phi(res.w.phi, res.wo.phi))
# plot.pred.vs.emp(data.set$mean.phis)



# estim.phi <- read.csv(paste(result.folder, "sdlog05_fittedIC.txt", sep=""))
# pdf(paste(result.folder, "metrics_sdlog05_fittedIC.pdf", sep=""), width=16, height=10)
# plot(estim.phi[, 2:5], main=expression(paste("est. ", bar(phi))))
# dev.off()




# ### COMPARE TWO RUNS
# w.phi <- read.csv(paste(result.folder, "wphi_no_low_conf_80k.txt", sep=""))
# wo.phi <- read.csv(paste(result.folder, "wophi_no_low_conf_80k.txt", sep=""))
# 
# pdf(paste(result.folder, "withphi_vs_fitIC_sdlog05", ".pdf", sep = ""), width = 10, height = 10)
# xy.max <- max(cbind(w.phi[,2], wo.phi[,2]))
# reg <- lm(wo.phi[, 2]~w.phi[, 2])
# plot(w.phi[, 2], wo.phi[,2], xlim=c(0, xy.max), ylim=c(0,xy.max),
#      main=expression(paste("with ", phi, " fit vs. fitted IC")),
#      sub=paste("R^2=", format(summary(reg)$r.squared, digits = 3), sep=""),
#      xlab=expression(paste("with obs. ", phi)), ylab=expression(paste("fit. I.C. ", phi)))
# abline(0,1, lty=2)
# abline(reg, col="red", lwd=2)
# dev.off()



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