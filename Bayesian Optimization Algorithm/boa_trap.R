#-----------------------------------------------------------------------------#
#                                  BOA                                        #
#                            BINARY PROBLEMS                                  #
#                                                                             #
# Ciniro Ap. Leite Nametala                                                   #
# ciniro@gmail.com                                                            #
# github.com/ciniro                                                           #
# V0.01 em 09/09/2019                                                         #
#-----------------------------------------------------------------------------#

#prepare environment
rm(list=ls())
cat('\014')

#libraries
library(bnlearn)

#functions--------------------------------------------------
f.trapk <- function(cromo) {
  n <- length(cromo)
  sizebb <- n/qtdebb
  
  fitness <- 0

  for(i in seq(1, n, by = sizebb)) {
    bb <- cromo[i:(i + sizebb - 1)]
    
    ones <- sum(bb)
    
    fitbb <- 0
    if (ones == sizebb)
      fitbb <- sizebb
    else
      fitbb <- ((sizebb - 1) - ones)
    
    fitness <- fitness + fitbb
  }
  
  return(fitness)
}

generate_initial_model <- function(size) {
  xi <- paste("x", 1:size, sep="")
  bn <- empty.graph(nodes = xi)
  
  lv <- c("0","1")
  
  cpt <- list()
  
  for (i in 1:size) {
    cpt[[i]] <- array(c(0.5, 0.5), 
                    dim = 2, 
                    dimnames = list(x = lv))
  }
  
  names(cpt) <- xi
  
  model = custom.fit(bn, cpt)
  
  return(list(model = model, bn = bn))
}

get_best <- function(pop) {
  mpop <- data.matrix(pop)-1
  fits <- apply(mpop, 1, fitfunction)
  indbest <- which.max(fits)
  
  return(list(best = mpop[indbest,], fit_best = fits[indbest]))
}

tournament <- function(pop, npop, k, n) {
  
  winners <- vector()
  fits <- vector()
  
  nwinners <- ceiling(npop*n)
  
  for (i in 1:nwinners) {
    ind_competitors <- sample(1:npop, k)
    winner <- get_best(pop[ind_competitors,])
    winners <- rbind(winners, as.numeric(winner$best))
    fits <- c(fits, as.numeric(winner$fit_best))
  }
  
  return (list(winners, fits))
}

build_bn <- function(dfpop, size, previous_bn) {
  dfpop <- apply(dfpop, 2, as.character)
  #include levels (at least one 0 e one 1)
  dfpop <- rbind(dfpop, rep("0", size), rep("1", size))
  dfpop <- as.data.frame(dfpop)
  xi <- paste("x", 1:size, sep="")
  colnames(dfpop) <- xi

  bn <- hc(dfpop, score = "k2", start = previous_bn, optimized = TRUE)
  
  model = bn.fit(bn, dfpop, replace.unidentifiable = TRUE)
  return(list(model = model, bn = bn))
}

#algorithm--------------------------------------------------
#parameters of selection
n.tournament <- 1 #percent
k.tournament <- 10

#size of cromossomo must be divisible by qtdebb
size <- 6
qtdebb <<- 3

#size <- 4
#qtdebb <<- 2

#parameters of convergence
maxgen <- 50
npop <- 100
fitfunction <- f.trapk
tol <- 0

#inicializes randomly a BN
best_bn <- generate_initial_model(size)

#generate initial population
pop <- rbn(best_bn$model, n = npop)

#get initial best cromosoma
bestcromo <- get_best(pop)

#storage best fitness by iteration
bestfitness <- rep(NA, maxgen)

for (iter in 1:maxgen) {
  cat(iter)
  cat('\014')
  
  #let them compete
  winners <- tournament(pop, npop, k.tournament, n.tournament)
  
  ## convergence verification criteria
  #tolerance
  if (all( abs(as.vector(winners[[1]]) - mean(as.vector(winners[[1]]))) < tol ) == TRUE)
   break()
  #variance
  if  (iter > 5)
    if (length(unique(c(bestfitness[iter - 1],
                        bestfitness[iter - 2],
                        bestfitness[iter - 3],
                        bestfitness[iter - 4]))) == 1)
      break()

  #building a new BN based on winners
  best_bn <- build_bn(winners[[1]], size, best_bn$bn)

  #update population based on new BN
  pop <- rbn(best_bn$model, n = npop)
  
  #evaluate best against winner
  newbestcromo <- get_best(pop)
  
  if (newbestcromo$fit_best > bestcromo$fit_best) {
    bestcromo <- newbestcromo
  }
  
  bestfitness[iter] <- newbestcromo$fit_best

}

plot(bestfitness, type="l", xlab="Iterations", 
     ylab="Best fitness", 
     main="Convergence", panel.first=grid(), lwd=2)

abline(v=iter-1, col="blue", lty=2)

# print("BEST MODEL BLOCKS:")
# print(best_bn$bn)
# print(best_bn$model)
graph::plot(best_bn$bn)
print("STOP ITERATION:")
print(iter)
print("BEST SOLUTION:")
print(bestcromo$best)
print("FIT BEST SOLUTION:")
print(bestcromo$fit_best)
print("K2 SCORE - BEST MODEL:")
suppressWarnings(print(score(best_bn$bn, pop, type = "k2")))
print("STRING MODEL:")
print(modelstring(best_bn$model))
