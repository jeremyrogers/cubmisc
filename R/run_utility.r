
mapBMatNames <- function(in.names, aa.names)
{
  codon.count <- list(A=3, R=5, N=1, D=1, C=1, Q=1, E=1, G=3, H=1, I=2, L=5, K=1, F=1, P=3, S=3, Z=1, T=3, Y=1, V=3, M=0, W=0)
  id.intercept <- grep("Intercept", in.names)
  id.slope <- 1:length(in.names)
  id.slope <- id.slope[-id.intercept]
  
  out.names <- vector(mode = "character", length=length(in.names))
  start <- 1
  for(aa in aa.names)
  {
    ncodons <- codon.count[[aa]]*2
    if(ncodons == 0) next # for M and W
    out.names[start:(start+ncodons-1)] <- rep(aa, ncodons)
    start <- start + ncodons
  }
  out.names[id.intercept] <- paste(out.names[id.intercept], "log(mu)")
  out.names[id.slope] <- paste(out.names[id.slope], "delta_t")
  return(out.names)
}



get.logL <- function(ret, data, model="nsef")
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
                        tmp <- .cubfitsEnv$my.logLAll(xx, phi.Obs, data$y, data$n, bInit,
                                                      ret$p.Mat[[i.iter]],
                                                      reu13.df = data$reu13.df)
                        sum(tmp)
                      }, mc.cores=4)
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
                     }, mc.cores=4)
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
