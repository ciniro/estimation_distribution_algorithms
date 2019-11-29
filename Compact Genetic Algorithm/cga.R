#-----------------------------------------------------------------------------#
#                       COMPACT GENETIC ALGORITHM                             #
#                            BINARY PROBLEMS                                  #
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
library("data.table")
library("pheatmap")

#-----------------------------------------------------------------------------#
#                               FUNCTIONS                                     #
#                                                                             #
#                     CINIRO APARECIDO LEITE NAMETALA                         #
#-----------------------------------------------------------------------------#

f.onemax <- function(cromo) {
  return (sum(cromo))
}

f.onemin <- function(cromo) {
  return (length(which(cromo==0)))
}

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
  
  if (fit_cromo1 > fit_cromo2) {
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
    if (winner[[1]][i] != loser[[1]][i]) {
      if (winner[[1]][i] == 1) {
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

size <- 10
maxgen <- 50
npop <- 10
fitfunction <- f.onemax

#inicializes a probability vector
probvector <- rep(0.5, size)

#generate initial best cromosoma
best <- generate(probvector)
fit_best <- fitfunction(best)
best <- list(best, fit_best)

#storage best fitness by iteration
bestfitness <- rep(NA, maxgen)
matprob <- probvector

for (iter in 1:maxgen) {
  #generates a new candidate solution based on the probability vector
  cromo1 <- generate(probvector)
  cromo2 <- generate(probvector)
  
  #let them compete, so we can know who is the best of the pair
  result <- tournament(cromo1, cromo2)
  winner <- result[[1]]
  loser <- result[[2]]
  
  #evaluate best against winner
  if (winner[[2]] > best[[2]]) {
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

plot(bestfitness, type="l", xlab="Iterations", 
     ylab="Best fitness", 
     main="Convergence", panel.first=grid(), lwd=2)

rownames(matprob) <- as.character(seq(1:dim(matprob)[1]))
pheatmap(matprob, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", color = colorRampPalette(c("white", "navy", "black"))(100), border_color = "grey60", breaks = seq(0,1,0.01), legend_breaks = seq(0,1,0.1), show_colnames = FALSE, show_rownames = FALSE, main = "Heatmap convergence")

print("STOP ITERATION:")
print(iter)
print("BEST SOLUTION:")
print(best[[1]])
print("FIT BEST SOLUTION:")
print(best[[2]])
print("PROBABILITY VECTOR:")
print(probvector)

