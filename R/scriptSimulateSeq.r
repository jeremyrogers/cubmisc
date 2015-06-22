setwd("~/cubmisc/jeremyYeast/")
rm(list=ls())
require(plyr, lib.loc = "~/cubfits/Dependencies/")
require(parallel)
require(NSEexchange) 
#suppressMessages(library(NSEexchange, lib.loc="~/NSEexchangeBuild", quietly=TRUE))

SaveParams <- TRUE
GenOptPes <- FALSE

data.folder <- "~/cubmisc/jeremyYeast/"
dna.file <- "~/cubmisc/jeremyYeast/Genome/jeremyYeast.fasta";
phi.file <- paste(data.folder, "jyeast.phi.csv", sep="")
#nse.pr.file <- paste(data.folder, "codon/NSE_pr/S.cerevisiae.2007.NSE.csv", sep="")
nse.pr.file <- paste(data.folder, "codon/NSE_pr/modded_NsePr.tsv", sep="")
mut.file <- paste(data.folder, "codon/Mut_rate/PNAS2011_S.cer.mut.tsv", sep="")
out.fasta <- paste(data.folder, "Genome/jeremyYeast2.fasta", sep="")

codon.params <- load.codon.parms(init.nse.pr.file = nse.pr.file, init.mut.file = mut.file, obs.nse.pr.file = nse.pr.file, obs.mut.file = mut.file)
phi <- load.phi(init.phi.file = phi.file, obs.phi.file = phi.file)
obs.genome <- load.genome(genome.fasta = dna.file, codon.parms = codon.params, phi = phi$init)

if(SaveParams){
  save(codon.params, file="S.Cer.codon.params_test.rda")
  save(obs.genome, phi, file="S.Cer.genome_test.rda")
}

## clean up workspace for better overview during code development 
rm(list=c("dna.file", "mut.file", "nse.pr.file", "phi.file"))
seq.ids <- mclapply(1:length(obs.genome), function(i) {paste(obs.genome[[i]]$name, collapse = "")}, mc.cores=4)


if(GenOptPes){
  ## generate optimal an pesimal sequence
  opt.seqs <- mclapply(1:length(obs.genome), function(i) {paste(codon.params$codon[opt.seq(obs.genome[[i]], codon.params)], collapse = "") }, mc.cores=4 )
  write.fasta(sequences = opt.seqs, names = seq.ids, file.out = "yeast_opt_seqs.fasta")
  pes.seqs <- mclapply(1:length(obs.genome), function(i) {paste(codon.params$codon[pess.seq(obs.genome[[i]], codon.params)], collapse = "") }, mc.cores=4 )
  write.fasta(sequences = pes.seqs, names = seq.ids, file.out = "yeast_pes_seqs.fasta")
}

## simulate sequence
obs.genome.parms <- list(A.1 = 4, A.2 = 4, B = 0.0025, Ne = 1, Q = 1)
#sim.genome <- simulate.data.all.genes(obs.genome = obs.genome, obs.phi = phi, obs.codon.parms = codon.params$obs, obs.genome.parms = obs.genome.parms, BIS=4, GMT=0, MES=0, n.cores=4, aux.simulation.type = 'M')
sim.genome <- simulate.data.all.genes(parallel='lapply', obs.genome = obs.genome, obs.phi = phi$obs, obs.codon.parms = codon.params$obs, obs.genome.parms = obs.genome.parms, BIS=4, GMT=0, MES=0, n.cores=4, aux.simulation.type = 'M')
sim.seqs <- mclapply(1:length(sim.genome), function(i) {paste(sim.genome[[i]]$gene.dat$codon, collapse = "")}, mc.cores=4)
write.fasta(sequences = sim.seqs, names = seq.ids, file.out = out.fasta)
