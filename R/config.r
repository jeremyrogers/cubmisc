print.config <- function(config)
{
  cat("MCMC parameters:\n")
  cat(paste("Number of Iterations:", config$n.iter, "\n"))
  cat(paste("Min iterations:", config$min.iter, "\n"))
  cat(paste("Max iterations:", config$max.iter, "\n"))
  cat(paste("Reset QR until:", config$reset.qr, "Iterations\n"))
  cat(paste("Thining: store every", config$chain.thin, "iteration\n"))
  
  cat("\n")
  
  cat("Simulation parameters\n")
  cat(paste("Number of Cores", config$n.cores, "\n"))
  cat(paste("Number of Chains:", config$n.chains, "\n"))
  cat(paste("Parallel mode within chain:", config$parallel, "\n"))
  cat(paste("Iterations used:", config$use.n.iter, "\n"))
  cat(paste("First", config$rm.first.aa, " AAs removed from sequence\n"))
  cat(paste("Sequences with less than", config$rm.short, "AAs ignored\n"))
  cat("List of AAs taken into account:\n\t");cat(config$aa);cat("\n")
  
  cat("\n") 
  
  cat("Convergence criteria\n")
  if(config$n.chains < 2)
  {
    cat("Convergence test: Geweke\n")
    cat(paste("Convergence criterium: Z Score <", config$eps, "\n"))
    cat(paste("% chain used at the begining:", config$frac1, "\n"))
    cat(paste("% chain used at the end:", config$frac2, "\n"))
  }else{
    cat("Convergence test: Gelman & Rubin\n")
    cat(paste("Convergence criterium: Gelman Score <", config$eps, "\n"))
  }
  cat(paste("Use every", config$conv.thin, "iteration for convergence test\n"))
  
  cat("\n")
}


config <- list(
  n.iter  = 5000,  # steps after burn.in between convergence checks.
  use.n.iter = 1000, #window for convergence test (multichain) and window to calculate mean phi. IMPORTANT: Make sure this one is smaller than the thined chain, otherwise saving will crash!
  n.chains = 4, # num chains
  n.cores = 4, # total num of cpu cores (should be about 5*n.chains when using parallel method other then "lapply")
  selected.env = 1, # deprecated (us if more than one dataset is stored in csv, e.g. different conditions)
  min.iter=200000, # minimum steps each chain has to do (without thining), convergence criterium is ignored until min.iter is reached 
  max.iter=220000, # maximum steps for each chain (without thining). MCMC will be stoped if the chain reached max.iter iterations (convergence criterium is ignored)
  reset.qr=400, # stop resetting qr matrix when checking for convergence after this iteration (after thining)
  conv.thin=1,  # thining for convergence test (recommend 1 if chain.thin != 1, otherwise double thining)
  chain.thin=50, # thining of the chain during runtime
  rm.first.aa=0, # remove first rm.first.aa AAs (after the first codon which is expected to be the start codon)
  rm.short=0, # # ignore sequences with length < rm.short AAs after the first rm.first.aa AAs are removed
  parallel="lapply", # parallel method within chain (lapply = no parallelization within chain)
  eps=0.05, # I usually use 0.05 when n.chain > 1 and XXXX otherwise
  frac1=0.1, # % of chain at the beginning of the chain (use.n.iter window) used to calculate mean for convergence test
  frac2=0.5, # % of chain at the end of the chain (use.n.iter window) used to calculate mean for convergence test
  
  #aa = c("D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "W", "Y", "Z"),
  #aa = c("A", "V"), 
  aa = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "Z"), # AAs to take into account
  
  use.scuo = T # false means empirical data is used as initial conditions
)
