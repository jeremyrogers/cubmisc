suppressMessages(library(cubfits, quietly = TRUE))
suppressMessages(library(psych, quietly = TRUE))
##suppressMessages(library(Rmisc, quietly = TRUE))
suppressMessages(library(getopt, quietly = TRUE))

##load helper functions used to process current project.
## not general enough to include in package.
source("run_utility.r")
source("config.r")


##Process input from command line
# using getopt library

##define matrix for options and flags that are passed to "Rscript run_roc.r ..." command 
spec <- matrix(c(
  'cubmethod', 'c', 2, "character",
  'fnin' , 'f', 1, "character",
  'fnpin' , 'p', 2, "character",
  'fnout' , 'o', 1, "character",
  'sdlog' , 's', 2, "character"
), byrow=TRUE, ncol=4);


args <- commandArgs(trailingOnly = TRUE)
opt <- getopt(spec, args)
## print config file to log
print.config(config)


##Assigning options to variables
##Question: why do they have different names?

if(debug.code){
  ## hard coded variables to be used when debugging.
  ## to be removed later?
  
  ## I/O variables
  sdlog.phi.init <- 5 # has to be non 0 for cubappr
  fn.in <- "../data/ecoli_K12_MG1655_genome_filtered.fasta"
  fn.phi.in <- "../data/ecoli_X_obs.csv"
  fn.phi.out <- "../results/test/test.phi"
  fn.out <- "../results/test/test.dat"
  cubmethods <- "cubappr"
}else{    
  sdlog.phi.init <- opt$sdlog
  sdlog.phi.init <- as.double(unlist(strsplit(sdlog.phi.init, " ")))
  ## set method for either with Xobs, cubfits, or without Xobs, cubappr,
  ## NOTE: cross validation approach, cubpred, is not supported in this script
  cubmethods <- opt$cubmethod
  fn.in <- opt$fnin
  fn.phi.in <- opt$fnpin
  fn.phi.out <- paste(opt$fnout, ".phi", sep="")
  fn.out <- paste(opt$fnout, ".dat", sep="")
}

cat(paste("starting time:", Sys.time(), "\n"))
cat(paste("using", cubmethods, "\n"))



## cat vs. print?  Print to file?
cat(paste("reading sequences from file", fn.in, "\n"))
## readGenome function in cubfits package
seq.string <- readGenome(fn.in, config$rm.short, config$rm.first.aa)




## if using Xobs data 
if(cubmethods == "cubfits")
{
  cat(paste("reading gene expression measurements (Xobs) from file\n", fn.phi.in, "\nand compare to ORF list from FASTA file\n", fn.phi, "\n"))
  ret.list <- read.empirical.data(fn.phi.in, seq.string, config$selected.env)
  phi.obs <- ret.list$empirical
  phi.obs <- phi.obs[order(names(phi.obs))]
  seq.string <- ret.list$genome
  rm("ret.list")
  
  ## setting the sdlog.phi.init to sd(log(Xobs)) if sdlog.phi.init == 0. 
  for(i in 1:config$n.chains)
  {  
    sdlog.phi.init[i] <- ifelse(sdlog.phi.init[i] == 0, sdlog.phi.init[i] <- sd(log(phi.obs)), sdlog.phi.init[i]) 
  }
}

## MORE DETAILS!!!
cat("generating other data...\n")
data <- generate.data(seq.string, config$aa)


## 
cat("generate initial phi using ")
init.phi <- list()
length(init.phi) <- config$n.chains
if(config$use.scuo)
{
  cat("SCUO with sd(ln(phi)) values\n")
  scuo <- gen.scuo(seq.string, config$aa)
  scuo <- calc_scuo_values(scuo)$SCUO
  for(i in 1:config$n.chains)
  {
    cat( paste("\t", sdlog.phi.init[i],"\n") )
    randscuo <- scuo.random(scuo, meanlog = -sdlog.phi.init[i]^2 / 2, sdlog = sdlog.phi.init[i])
    randscuo <- randscuo / mean(randscuo)
    names(randscuo) <- names(seq.string)
    init.phi[[i]] <- randscuo
  }
}else{
  cat("Xobs data. Note: All chains starting with same initial phi values\n")
  if(cubmethods == "cubfits")
  {
    for(i in 1:config$n.chains){ init.phi[[i]] <- phi.obs }
  }else{
    ret.list <- read.empirical.data(fn.phi.in, seq.string, config$selected.env)
    phi.obs <- ret.list$empirical
    phi.obs <- phi.obs[order(names(phi.obs))]
    for(i in 1:config$n.chains){ init.phi[[i]] <- phi.obs }
  }
}

##define objects to hold initial values and output from MCMC chains
results <- list()
length(results) <- config$n.chains
p.init <- list(NULL)
length(p.init) <- config$n.chains

##Hard coding of initial values for pMat of each chain
## c(sd_epsilon, meanlog(phi), sdlog(phi), K_x scaling constant)
if(cubmethods == "cubfits"){
  p.init[[1]] <- c(0.5, -0.5, 1, 7)
  p.init[[2]] <- c(0.1, -0.03125, 0.25, 1)
  p.init[[3]] <- c(1.5, -2.0, 2, 5)
  p.init[[4]] <- c(1.0, -4.5, 3, 10)
} else if(cubmethods == "cubappr") {
  p.init[[1]] <- c(-0.5, 1)
  p.init[[2]] <- c(-0.03125, 0.25)
  p.init[[3]] <- c(-2.0, 2)
  p.init[[4]] <- c(-4.5, 3)
}

##Change to reflect single chain possibility
cat(paste("running ", cubmethods, " using multiple chains \n"))
seeds <- round(runif(config$n.chains, 1, 100000))

##RENAME!
s <- system.time(
  {
    .CF.CT$parallel <- "mclapply"
    if(cubmethods == "cubfits"){
      .CF.CT$type.p <- "lognormal_bias"
      .CF.CONF$scale.phi.Obs <- F
      .CF.CONF$estimate.bias.Phi <- T
      if(config$n.chains < 2)
      {
        results <- cubsinglechain(cubmethods, niter=config$use.n.iter, frac1=0.1, frac2=0.5, reset.qr=config$reset.qr, seed=seeds, teston="sphi",
                               min=config$min.iter, max=config$max.iter, thin=config$thin, eps=1, 
                               reu13.df.obs=data$reu13.df, phi.Obs=phi.obs, y=data$y, n=data$n, phi.Init=init.phi[[1]],
                               nIter=config$n.iter, burnin=config$burn.in, p.Init=p.init[[1]],
                               model="roc", adaptive="simple", .CF.CT=.CF.CT, .CF.CONF=.CF.CONF)
      }else{
        results <- cubmultichain(cubmethods, niter=config$use.n.iter, reset.qr=config$reset.qr, seeds=seeds, teston="sphi",
                               min=config$min.iter, max=config$max.iter, nchains=config$n.chains, thin=config$thin, eps=0.05, ncores=config$n.cores,
                               reu13.df.obs=data$reu13.df, phi.Obs=phi.obs, y=data$y, n=data$n, phi.Init=init.phi,
                               nIter=config$n.iter, burnin=config$burn.in, p.Init=p.init,
                               model="roc", adaptive="simple", .CF.CT=.CF.CT, .CF.CONF=.CF.CONF)
      }
    }else if(cubmethods == "cubappr") {
      if(config$n.chains < 2)
      {
        results <- cubsinglechain(cubmethods, niter=config$use.n.iter, frac1=0.1, frac2=0.5, reset.qr=config$reset.qr, seed=seeds, teston="sphi",
                               min=config$min.iter, max=config$max.iter, thin=config$thin, eps=1, 
                               reu13.df.obs=data$reu13.df, y=data$y, n=data$n, phi.pred.Init=init.phi[[1]],
                               nIter=config$n.iter, burnin=config$burn.in, p.Init=p.init[[1]],
                               model="roc", adaptive="simple", .CF.CT=.CF.CT, .CF.CONF=.CF.CONF)
      }else{
        results <- cubmultichain(cubmethods, niter=config$use.n.iter, reset.qr=config$reset.qr, seeds=seeds, teston="sphi",
                               min=config$min.iter, max=config$max.iter, nchains=config$n.chains, thin=config$thin, eps=0.05, 
                               ncores=config$n.cores, reu13.df.obs=data$reu13.df, y=data$y, n=data$n, phi.pred.Init=init.phi,
                               nIter=config$n.iter, burnin=config$burn.in, p.Init=p.init,
                               model="roc", adaptive="simple", .CF.CT=.CF.CT, .CF.CONF=.CF.CONF)        
      }
    }
  }
) 

cat(paste("Elapsed time for", config$n.chains, "runs on", config$n.cores, "cores was", round(s["elapsed"]/60, digits=2), "min\n" ))
seq.string.names <- names(seq.string)
rm("seq.string")



## process results
cat("process results...\n")
phi.pred <- list()
for(i in 1:length(results$chains))
{
  if(cubmethods == "cubfits"){
    phi.pred[[i]] <- do.call("cbind", results$chains[[i]]$phi.Mat) ### matrix n.G by nIter+1
  }
  if(cubmethods == "cubappr") {
    phi.pred[[i]] <- do.call("cbind", results$chains[[i]]$phi.pred.Mat) ### matrix n.G by nIter+1
  }
}

mean.phis <- list()
median.phis <- list()
geo.mean.phis <- list()
harm.mean.phis <- list()
sd.phi <- list()
interval <- (dim(phi.pred[[i]])[2]-config$use.n.iter):dim(phi.pred[[i]])[2]
for(i in 1:length(results$chains))
{
  mean.phis[[i]] <- rowMeans(phi.pred[[i]][, interval]) ### mean of the last 2000 iterations
  median.phis[[i]] <- apply(phi.pred[[i]][, interval], 1, median)
  geo.mean.phis[[i]] <- apply(phi.pred[[i]][, interval], 1, geometric.mean)
  harm.mean.phis[[i]] <- apply(phi.pred[[i]][, interval], 1, harmonic.mean)
  sd.phi[[i]] <- apply(phi.pred[[i]][, seq(interval[1], interval[length(interval)], by=config$use.n.iter/10)], 1, sd)
}
## save results and input
cat("saving results...\n")
mean.phis <- do.call("cbind", mean.phis)
median.phis <- do.call("cbind", median.phis)
geo.mean.phis <- do.call("cbind", geo.mean.phis)
harm.mean.phis <- do.call("cbind", harm.mean.phis)
sd.phi <- do.call("cbind", sd.phi)

if(config$n.chains > 1) # only need that if there is more than one run
{
  mean.phis <- rowMeans(mean.phis)
  median.phis <- apply(median.phis, 1, median)
  geo.mean.phis <- apply(geo.mean.phis, 1, geometric.mean)
  harm.mean.phis <- apply(harm.mean.phis, 1, harmonic.mean)
  sd.phi <- rowMeans(sd.phi)
}
means <- cbind(seq.string.names, mean.phis, median.phis, geo.mean.phis, harm.mean.phis, sd.phi)
colnames(means) <- c("ORF_Info", "Phi_Post_Arith_Mean", "Phi_Post_Median", "Phi_Post_Geom_Mean", "Phi_Post_Harm_Mean", "Phi_Post_SD") 
write.csv(means, file = fn.phi.out, row.names=F)
#write.table(means, file = fn.phi.out)

for(i in 1:config$n.chains)
{
  chain <- results$chains[[i]]
  convergence <- results$convergence
  list.save <- c("data", "mean.phis", "seeds", "chain", "seq.string.names", "config", "means", "sdlog", "convergence")
  save(list = list.save, file = paste(fn.out, "_chain_", i, sep=""))
  
}
rm("phi.pred")
