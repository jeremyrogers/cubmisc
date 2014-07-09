print.config <- function(config)
{
  cat("MCMC parameters:\n")
  cat(paste("Number of Iterations:", config$n.iter, "\n"))
  cat(paste("Burn In steps:", config$burn.in, "\n"))
  cat(paste("Iterations used:", config$use.n.iter, "\n"))
  
  cat("\n")
  
  cat("Simulation parameters\n")
  cat(paste("Number of Cores", config$n.cores, "\n"))
  cat(paste("Number of Runs:", config$n.runs, "\n"))
  cat(paste("Data coloumn selected:", config$selected.env, "\n"))
}


config <- list(
  n.iter = 200,
  burn.in = 100,
  use.n.iter = 1000,
  n.runs = 4,
  n.cores = 4,
  selected.env = 1,
  min.iter=1000,
  max.iter=10000,
  reset.qr=1000,
  thin=100,
  
  #aa = c("D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "W", "Y", "Z"),
  #aa = c("A", "V"), 
  aa = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "Z"),
  
  use.scuo = T, # false means empirical data is used as initial conditions
  randomize.scuo = T
)
