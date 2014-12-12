### SEE ALSO: errors.r ###

OmegaStudy <- function(omega, cutoff=.001, Title="Omega Study", Subtitle=NULL, details="raw error")
{
  
if(is.null(Subtitle)){Subtitle <- paste("Error Range is +/-", cutoff)}
  
m <- rbind(c(1, 1, 2, 2), c(1, 1, 2, 2), c(3, 3, 4, 4), c(3, 3, 5, 5))
layout(m)
par(oma = c(2, 0, 3, 0))

ci <- .95

distances <- omega[,2]-omega[,1]
#distances <- 100*(omega[,2]-omega[,1])/abs(omega[,1]); cutoff=400

bad <- which(abs(distances) > cutoff)
omega.bad <- omega[bad,]
omega.good <- omega[-bad,]
reg <- lm(omega[,2]~omega[,1]);
reg.good <- lm(omega.good[,2]~omega.good[,1]);
reg.bad <- lm(omega.bad[,2]~omega.bad[,1]);

#plot regular
{
  plot(omega, main=bquote(.(details) ~ "'true'" ~ Delta * omega ~ "vs estimated" ~ Delta * omega)
       ,ylim=range(omega[,3:4])
       ,xlab=expression(paste("'true' ", Delta, omega))
       ,ylab=expression(paste("estimated ", Delta, omega))
       ,sub=paste("R^2=", format(summary(reg)$r.squared, digits = 3), sep="") 
  )
  points(omega[bad,], bg="black", pch=21)
  abline(reg, lwd=2, col="red")
  abline(0,1, lwd=1, col="lightslategrey")
  abline(cutoff,1, lwd=1, lty="dotted", col="lightslategrey")
  abline(-cutoff,1, lwd=1, lty="dotted", col="lightslategrey")

  plotCI(x = omega[,1], y=omega[,2], ui=omega[,4], li=omega[,3], add=TRUE, scol = "darkslategrey");
  
  legend("topleft", legend = c("1-1 line", "Error Range", "Linear Fit", paste(100*ci, "% CI", sep=""))
         ,lty=c("solid", "dotted", "solid", "solid")
         ,col=c("lightslategrey", "lightslategrey", "red", "darkslategrey")
         #,bty="n"
  )
}

#plot good
{
  plot(omega.good, main=bquote(.(details) ~ Delta * omega ~ "without problem values")
#       ,xlim=range(omega[,1])
#       ,ylim=range(omega[,3:4])
       ,ylim=range(omega.good[,3:4])
       ,xlab=expression(paste("'true' ", Delta, omega))
       ,ylab=expression(paste("estimated ", Delta, omega))
       ,sub=paste("R^2=", format(summary(reg.good)$r.squared, digits = 3), sep="") 
  )
  abline(reg, lwd=2, col="red")
  abline(reg.good, lwd=1, col="tomato", lty="dashed")
  abline(0,1, lwd=1, col="lightslategrey")
  plotCI(x = omega.good[,1], y=omega.good[,2], ui=omega.good[,4], li=omega.good[,3], add=TRUE, scol = "darkslategrey");
  
  legend("topleft", legend = c("1-1 line", "Fit Good Values", "Fit All Values", paste(100*ci, "% CI", sep=""))
         ,lty=c("solid", "solid", "dashed", "solid")
         ,col=c("lightslategrey", "red", "tomato", "darkslategrey")
  )
}

#plot bad
{
  plot(omega.bad, main=bquote(.(details) ~ Delta * omega ~ " problem values")
       ,xlim=range(omega[,1])
       ,ylim=range(omega[,3:4])
       ,xlab=expression(paste("'true' ", Delta, omega))
       ,ylab=expression(paste("estimated ", Delta, omega))
       ,sub=paste("R^2=", format(summary(reg.bad)$r.squared, digits = 3), sep="") 
  )
  abline(reg.bad, lwd=2, col="red")
  abline(reg, lwd=1, col="tomato", lty="dashed")
  abline(0,1, lwd=1, col="lightslategrey")
  plotCI(x = omega.bad[,1], y=omega.bad[,2], ui=omega.bad[,4], li=omega.bad[,3], add=TRUE, scol = "darkslategrey");
  
  legend("topleft", legend = c("1-1 line", "Fit Bad Values", "Fit All Values", paste(100*ci, "% CI", sep=""))
         ,lty=c("solid", "solid", "dashed", "solid")
         ,col=c("lightslategrey", "red", "tomato", "darkslategrey")
  )
}

off <- par()$mar[3]-2

breaks=seq(from=min(distances), to=max(distances), length.out=15)
info <- hist(distances, ,breaks=breaks, plot=FALSE)
par(mar=par()$mar+c(-1,0,-off,0)) #make extra room on top and bottom for the lack of a title
hist(distances
     ,main=NULL
     ,xlab=details
     ,ylim=c(0, max(info$counts)*6/5)
     ,labels=TRUE
     ,breaks=breaks
)

par(mar=par()$mar+c(1,0,0,0))

plot(x=omega[,1], y=distances, ylab=details, xlab=expression(paste("'true' ", Delta, omega)))
lines(x = range(omega[,1]),y = c(cutoff,cutoff))
lines(x = range(omega[,1]),y = c(-cutoff,-cutoff))

par(mar=par()$mar+c(0,0,off,0))

mtext(text = Title, outer=TRUE, side = 3, cex = 1.5)
mtext(text = Subtitle, outer = TRUE, side = 1)

}

init <- function()
{
  library(cubfits, quietly = TRUE, lib.loc = "~/cubfitsNSEdebug/")
  setwd("~/cubfits/misc/R/")
  source("visualize_utility.r")
  source("run_utility.r")
  source("config.r")
  
  model <- 'nse';              #nse \ roc
  genome <- 'pYeast';       #ecoli \ REUyeast \ pYeast \ brewYeast
  #prefix <- "11-21"            #generally, date the run started in "MM-DD" format
  #suffix <- "6000.1a2.Mflip"             #What was special about this run?
  prefix <- "11-24"            #generally, date the run started in "MM-DD" format
  suffix <- "6000.4a2.noMflip"             #What was special about this run?
  
  simulated_data <- TRUE;
  
  #### read data
  result.folder <- paste("../results/", substr(tolower(model),1,1), substr(tolower(genome),1,1), "/", prefix, "/", sep="")
  load(paste(result.folder, genome, ".", model, ".", suffix, ".dat" , sep=""))
  estm.phi <- read.csv(paste(result.folder, genome, ".", model, ".", suffix, ".phi" , sep=""))
  
  seq.string <- readGenome("../prestonYeast/S.cerevisiae.S288c.fasta", config$rm.short, config$rm.first.aa)
  emp <- read.empirical.data("../prestonYeast/pyeast.phi.tsv", seq.string, config$selected.env, th=0)
  omega_true <- read.table("../prestonYeast/scaled_omega.csv", header = TRUE)
  names(chain$b.Mat[[1]]) <- read.table("names.txt")[[1]]
  ###this needs to be fixed!!!
  
  require(package = plotrix, lib.loc = "/home/lbrown/cubfits/Dependencies/")
  chunksize <- floor(length(chain$b.Mat)/100)*10
  ci <- .95 #Confidence interval, on a scale of (0,1). e.g. 0.95 means 95%CI
  
  #Generate Omegas
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
  }##finish all codons
  
  omega
}