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
