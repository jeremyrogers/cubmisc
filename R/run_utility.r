## set mean of dataset to 1
normalize.data.set <- function(data)
{
  data <- data/mean(data)
}


# get.init.phi <- function(use.scuo=T, fn.phi.in, seq.string, aa.names, sdlog=0.5, randomize=T)
# {
#   scuo <- 0
#   if(use.scuo)
#   {
#     scuo <- get.scuo(genome=seq.string, aa.names=aa.names, randomize=T, sdlog=sdlog) 
#   }else{
#     estim.phi <- read.csv(fn.phi.in)
#     names <- estim.phi[, 1]
#     estim.phi <- estim.phi[, config$selected.env + 1]
#     names(estim.phi) <- names
#     scuo <- normalize.data.set(estim.phi)
#   }
#   return(scuo)
# }


get.logL <- function(ret, data, model="roc")
{
  require(multicore)
  if(is.null(ret[["phi.Mat"]])){
    ret$phi.Mat <- ret$phi.pred.Mat
  }
  
  if("sigmaW" %in% rownames(ret$p.Mat[[1]])){
    logL <- mcllapply(1:length(ret$phi.Mat),
                      function(i.iter){
                        xx <- ret$phi.Mat[[i.iter]]
                        bInit <- convert.bVec.to.b(ret$b.Mat[[i.iter]],
                                                   names(data$reu13.df),
                                                   model = model)
                        bInit <- lapply(bInit, function(B) B$coefficients)
                        sigmaWsq <- ret$p.Mat[[i.iter]][1]^2
                        tmp <- .cubfitsEnv$my.logLAll(xx, phi.Obs, data$y, data$n, bInit,
                                                      sigmaWsq,
                                                      reu13.df = data$reu13.df)
                        sum(tmp)
                      }, mc.cores=2)
  }
  else
  {
    logL <- mclapply(1:length(ret$phi.Mat),
                     function(i.iter){
                       xx <- ret$phi.Mat[[i.iter]]
                       bInit <- convert.bVec.to.b(ret$b.Mat[[i.iter]],
                                                  names(data$reu13.df),
                                                  model = model)
                       bInit <- lapply(bInit, function(B) B$coefficients)
                       tmp <- .cubfitsEnv$my.logLAllPred(phi=xx, y=data$y, n=data$n, b=bInit,
                                                         reu13.df = data$reu13.df)
                       sum(tmp)
                     }, mc.cores=2)
  }
}

# reading.genome <- function(fn.genome, ex.sh.aa=-1)
# {
#   cat("reading sequences...\n")
#   seq.data <- read.seq(fn.genome)
#   ind <- unlist(lapply(1:length(seq.data), function(i)
#           {
#             return(length(seq.data[[i]]) > (ex.sh.aa*3))
#           }))
#   seq.data <- seq.data[ind]
#   
#   cat("converting sequences...\n")
#   seq.string <- convert.seq.data.to.string(seq.data)
#   seq.string <- seq.string[order(names(seq.string))]
#   
#   return(seq.string)
# }

read.empirical.data <- function(fn, genome, env, th=-1)
{
  ## load empirical data
  empirical <- read.table(file=fn, header=T, sep=",", fill=T)
  gene.names <- empirical[,1]
  rownames(empirical) <- gene.names
  
  empirical <- empirical[gene.names %in% names(genome), env + 1]
  gene.names <- gene.names[gene.names %in% names(genome)]
  names(empirical) <- gene.names
  empirical <- empirical[empirical > th]
  empirical <- empirical[!is.na(empirical)]
  
  genome <- genome[names(genome) %in% names(empirical)]
  return(list(empirical=empirical, genome=genome))
}

generate.data <- function(genome, aa.names)
{
  ## generate data set
  phi <- data.frame(names(genome), rlnorm(length(genome), 0, 1))
  reu13.df <- gen.reu13.df(genome, phi, aa.names)
  
  y <- gen.y(genome, aa.names)
  n <- gen.n(genome, aa.names)
  return(list(reu13.df=reu13.df, y=y, n=n))
}

# get.scuo <- function(genome, aa.names, randomize=F, sdlog=1.5)
# {
#   scuo <- gen.scuo(genome, aa.names)
#   scuo <- calc_scuo_values(scuo)$SCUO
#   if(randomize)
#   {
#     scuo <- scuo.random(scuo, meanlog = -sdlog^2 / 2,
#                         sdlog = sdlog)  
#   }
#   scuo <- scuo / mean(scuo)
#   names(scuo) <- names(genome)
#   return(scuo)
# }