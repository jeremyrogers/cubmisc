print.config <- function(config)
{
  cat("MCMC parameters:\n")
  cat(paste("Number of Iterations:", config$n.iter, "\n"))
  cat(paste("Burn In steps:", config$burn.in, "\n"))
  cat(paste("Min iterations:", config$min.iter, "\n"))
  cat(paste("Max iterations:", config$max.iter, "\n"))
  cat(paste("Thin every:", config$thin, "Iterations\n"))
  cat(paste("Reset QR until:", config$reset.qr, "Iterations\n"))
  
  cat("\n")
  
  cat("Simulation parameters\n")
  cat(paste("Number of Cores", config$n.cores, "\n"))
  cat(paste("Number of Chains:", config$n.chains, "\n"))
  cat(paste("Iterations used:", config$use.n.iter, "\n"))
  cat(paste("First", config$rm.first.aa, " removed\n"))
  cat(paste("Sequences with less than", config$rm.short, " AA ignored\n"))
  
  cat("\n") 
}


config <- list(
  n.iter  = 500,  # steps after burn.in between convergence checks.
  burn.in = 0, # burn.in after restarting
  use.n.iter = 2000, #window for convergence test
  n.chains = 4, # num chains
  n.cores = 4, # total num of cpu cores
  selected.env = 1, #deprecated
  min.iter=1000, 
  max.iter=50000,
  reset.qr=2200, #stop resetting qr matrix when checking for convergence after this iteration  
  thin=20,
  rm.first.aa=0,
  rm.short=0,
  
  #aa = c("D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "W", "Y", "Z"),
  #aa = c("A", "V"), 
  aa = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "Z"),
  
  use.scuo = T # false means empirical data is used as initial conditions
)
