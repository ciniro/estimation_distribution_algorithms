#-----------------------------------------------------------------------------#
#                       COMPACT GENETIC ALGORITHM                             #
#                             REAL PROBLEMS                                   #
#                                                                             #
# Ciniro Ap. Leite Nametala                                                   #
# ciniro@gmail.com                                                            #
# github.com/ciniro                                                           #
# V0.02 em 21/08/2019                                                         #
#-----------------------------------------------------------------------------#

#clean workspace
rm(list=ls())

#clean screen
cat('\014')

#libraries
library("plot3D")
library("stringr")
library("rgl")
library("magick")
library("data.table")
library("pheatmap")

#-----------------------------------------------------------------------------#
#                               FUNCTIONS                                     #
#                                                                             #
#                     CINIRO APARECIDO LEITE NAMETALA                         #
#-----------------------------------------------------------------------------#

#fitness functions
#sphere
#optim: 0..
f.sphere <- function(x) 
{
  return(sum(x^2))
}

#rastrigin
#optim: 0..
f.rastrigin <- function(x)
{
  d <- length(x)
  sum <- sum(x^2 - 10*cos(2*pi*x))
  y <- 10*d + sum
  return(y)
}

#rosenbrock
#optim: 1..
f.rosenbrock <- function(x)
{
  d <- length(x)
  xi <- x[1:(d-1)]
  xnext <- x[2:d]
  
  sum <- sum(100*(xnext-xi^2)^2 + (xi-1)^2)
  
  y <- sum
  return(y)
}

#graph functions
generateContours <- function(pop, x = seq(0, 1, length.out = nrow(z)),
                               y = seq(0, 1, length.out = ncol(z)),
                               z, nlevels=30, ...) {
  ma <- max(abs(z))
  lvls <- seq(-ma, ma, length.out = nlevels)
  cols <- colorRampPalette(c("blue","white","red")) (nlevels - 1)
  filled.contour(x, y, z, main = "Optimization point", plot.axes = {axis(1); axis(2); points(x = pop[1], y = pop[2], pch=20, col="blue")},
                 col=cols, levels=lvls, ...)
}

#representation functions
BinToDec <- function(x) {
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
}

DecToBin <- function(x){
  tmp <- rev(as.numeric(intToBits(x)))
  id <- seq_len(match(1,tmp,length(tmp))-1)
  tmp[-id]
}

BinToReal <- function (bit, p, linf) {
  dec <- BinToDec(bit)
  dec <- linf + (p * dec)
  return (dec)
}

RealToBin <- function (dec, p, linf, bitnec) {
  bit <- (dec-linf)/p
  vet <- paste(round(DecToBin(bit)), collapse = "")
  vet <- str_pad(vet, bitnec, pad = "0")
  return(vet)
}

contaCasasDec <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

splitChunk <- function(x, n) {
  sst <- strsplit(x, '')[[1]]
  m <- matrix('', nrow=n, ncol=(length(sst)+n-1)%/%n)
  m[seq_along(sst)] <- sst
  apply(m, 2, paste, collapse='')
}

convertCromoBinToReal <- function(cromo, bitnec, size){
  n <- length(cromo)
  real <- rep(NA, size)
  strcromo <- paste(as.character(cromo), collapse = '')
  strcromo <- splitChunk(strcromo, bitnec)
  
  for (i in 1:size)
    real[i] <- BinToReal(strcromo[i], precisao, linf)
  
  return(real)
}

convertCromoRealToBin <- function(cromo){
  n <- length(cromo)
  bin <- rep(NA, n)
  for (i in 1:n)
    bin[i] <- RealToBin(cromo[i], precisao, linf, bitnec)
  
  return(bin)
}

#algorithm functions
generate <- function(probvector) {
  
  candidate <- rep(0, length(probvector))
  
  i <- 1
  for (p in probvector) {
    if (runif(1) < p)
      candidate[i] <- 1
    else
      candidate[i] <- 0
    
    i <- i + 1
  }
  
  return(candidate)
}

tournament <- function(cromo1, cromo2) {
  fit_cromo1 <- fitfunction(cromo1)
  fit_cromo2 <- fitfunction(cromo2)
  
  #let them compete, so we can know who is the best of the pair
  winner <- NA
  loser <- NA
  
  if (fit_cromo1 < fit_cromo2) {
    winner <- list(cromo1, fit_cromo1)
    loser <- list(cromo2, fit_cromo2)
  } else {
    winner <- list(cromo2, fit_cromo2)
    loser <- list(cromo1, fit_cromo1)
  }
  
  return (list(winner, loser))
}

update <- function(probvector, winner, loser, npop) {
  for (i in 1:length(probvector)) {
    w <- convertCromoRealToBin(winner[[1]])
    w <- as.numeric(strsplit(paste(w,collapse=""),"")[[1]])
    l <- convertCromoRealToBin(loser[[1]])
    l <- as.numeric(strsplit(paste(l,collapse=""),"")[[1]])
    
    if (w[i] != l[i]) {
      if (w[i] == 1) {
        probvector[i] <- probvector[i] + (1/npop)
      } else {
        probvector[i] <- probvector[i] - (1/npop)
      }
    }
  }
  
  return(round(probvector,2))
}

convergence <- function(probvector) {
  conv <- TRUE
  for (p in probvector) {
    if (p > 0 & p < 1) {
      conv <- FALSE
      break
    }
  }
  
  return(conv)
}

#-----------------------------------------------------------------------------#
#                               ALGORITHM                                     #
#                                                                             #
#                     CINIRO APARECIDO LEITE NAMETALA                         #
#-----------------------------------------------------------------------------#
#algorithm
size <- 2
maxgen <- 2000
npop <- 100
fitfunction <- f.rastrigin

#representation
linf <- -2
lsup <- 2
precisao <- 0.1
bitnec <- ceiling(log2((lsup-linf)/precisao))
casas <- contaCasasDec(precisao)
lsupbin <- RealToBin(lsup, precisao, linf, bitnec)
linfbin <- RealToBin(linf, precisao, linf, bitnec)

#graphs
plotGraphs <- TRUE
resPlotFunc <- 0.1
plotFuncRotate <- TRUE

#plot function fitness
if (plotGraphs == TRUE)
{
  #plotando a funcoes
  x <- seq(linf,lsup,resPlotFunc)
  y <- x
  z <- matrix(NA, length(x), length(y))
  for (i in 1:length(x))
  {
    for (j in 1:length(y))
    {
      z[i,j] <- fitfunction(c(x[i],y[j]))
    }
  }
  
  if (plotFuncRotate == TRUE) {
    persp3D(x,y,z,theta=45,phi=30)
    
    rgl::persp3d(x,y,z, col = "lightblue",
                 ticktype="detailed", xlab="", ylab="", zlab="",axes=FALSE)
  }
}

#inicializes a probability vector
probvector <- rep(0.5, bitnec*size)

#generate initial best cromosoma
best <- convertCromoBinToReal(generate(probvector), bitnec, size)
fit_best <- fitfunction(best)
best <- list(best, fit_best)

#storage best fitness by iteration
bestfitness <- rep(NA, maxgen)
matprob <- probvector

for (iter in 1:maxgen) {
  #generates a new candidate solution based on the probability vector
  cromo1 <- convertCromoBinToReal(generate(probvector), bitnec, size)
  cromo2 <- convertCromoBinToReal(generate(probvector), bitnec, size)
  
  #let them compete, so we can know who is the best of the pair
  result <- tournament(cromo1, cromo2)
  winner <- result[[1]]
  loser <- result[[2]]
  
  #evaluate best against winner
  if (winner[[2]] < best[[2]]) {
    best <- winner
  }
  bestfitness[iter] <- best[[2]]
  
  #update probability vector
  probvector <- update(probvector, winner, loser, npop)
  matprob <- rbind(matprob, probvector)
  
  #convergence verification
  if (convergence(probvector) == TRUE)
    break
}

print("STOP ITERATION:")
print(iter)
print("BEST SOLUTION (REAL):")
print(best[[1]])
print("BEST SOLUTION (BINARY):")
print(convertCromoRealToBin(best[[1]]))
print("FIT BEST SOLUTION:")
print(best[[2]])
print("PROBABILITY VECTOR:")
print(probvector)

if (plotGraphs == TRUE) {
  plot(bestfitness, type="l", xlab="Iterations", 
       ylab="Best fitness", 
       main="Convergence evolution", panel.first=grid(), lwd=2)
  
  generateContours(best[[1]],x,y,z)
  
  rownames(matprob) <- as.character(seq(1:dim(matprob)[1]))
  pheatmap(matprob, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", color = colorRampPalette(c("white", "navy", "black"))(100), border_color = "grey60", breaks = seq(0,1,0.01), legend_breaks = seq(0,1,0.1), show_colnames = FALSE, show_rownames = FALSE, main = "Heatmap convergence")
}