plot.likelihood.trace <- function(chain, data, windowsize)
{ 
  
  #logTraces <- get.logL(chain, data)
  logTraces <- unlist(chain$logL.Mat)
  ll <- format(mean(as.numeric(logTraces[(length(logTraces)-windowsize):length(logTraces)])), digits = 10)
  plot(1:length(logTraces), logTraces, 
       main=ll, xlab="Iteration", ylab="LogL", type="l")
  cat(paste(ll, "\n"))
  
  return(ll)
}


# load.data.file <- function(fn)
# {
#   load(fn)
#   return(list(res=results, data=data, mean.phis=mean.phis))
# }
# 
# find.optimal.layout <- function(n)
# {
#   if(n %% 2 != 0) {n <- n + 1}
#   for(i in 2:(n-1))
#   {
#     if(n %% i == 0) break
#   }
#   x <- n/i
#   if(x %% 2 == 0 & x > 2)
#   {
#     x <- x/2
#   }
#   return(list(x=x, y=n/x))
# }
# 
# plot.w.phi.vs.wo.phi <- function(w.phi, wo.phi, main="with phi vs without phi", wtf=F)
# {
#   if(length(w.phi$mean.phis) != length(wo.phi$mean.phis)) {stop("x and y data set of unequal length")}
#   
#   lay <- find.optimal.layout(length(wo.phi$mean.phis))
#   par(oma=c(1,1,2,1), mgp=c(2,1,0), mar = c(3,4,2,1), mfrow=c(lay$x, lay$y))
#   for(i in 1:length(wo.phi$mean.phis))
#   {
#     reg <- lm(log10(wo.phi$mean.phis[[i]])~log10(w.phi$mean.phis[[i]]))
#     if(wtf) # save plot directly as pdf
#     {
#       
#     }else{
#       plot(log10(w.phi$mean.phis[[i]]), log10(wo.phi$mean.phis[[i]]), xlab="log10(with phi)", ylab="log10(w/o phi)",
#            main=paste("R^2=", format(summary(reg)$r.squared, digits = 3), sep=""))
#       abline(reg, col="red", lwd=2)   
#     }
#   }
#   title(main=main,outer=T)
# }
# 
# plot.pred.vs.emp <- function(mean.phis, main="predicted vs empirical")
# {
#   seq.data <- read.seq(config$genome.data)
#   aa.names <- config$aa
#   
#   accessions <- scan(config$accession.data, what="", sep="\n")
#   seq.data <- seq.data[names(seq.data) %in% accessions]
#   
#   seq.string <- convert.seq.data.to.string(seq.data)
#   seq.string <- seq.string[order(names(seq.string))]
#   rm("seq.data", "accessions")
#   emp <- read.empirical.data(config$empirical.data, seq.string, config$selected.env)
#   emp <- normalize.data.set(emp$empirical)
#   
# 
#   lay <- find.optimal.layout(length(mean.phis))
#   par(oma=c(1,1,2,1), mgp=c(2,1,0), mar = c(3,4,2,1), mfrow=c(lay$x, lay$y))  
#   
#   for(i in 1:length(mean.phis))
#   { 
#     curr.mean.phis <- mean.phis[[i]][names(mean.phis[[i]]) %in% names(emp)]
#     reg <- lm(log10(emp)~log10(curr.mean.phis))
#     plot(log10(curr.mean.phis), log(emp), main=paste("R^2=", format(summary(reg)$r.squared, digits = 3), sep=""),
#          xlab="predicted phi", ylab="empirical phi")
#     abline(reg, col="red", lwd=2)  
#   }
#   title(main=main,outer=T)
# }


