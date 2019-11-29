#-----------------------------------------------------------------------------#
#                    EXTENDED COMPACT GENETIC ALGORITHM                       #
#                            BINARY PROBLEMS                                  #
#                                                                             #
# Ciniro Ap. Leite Nametala                                                   #
# ciniro@gmail.com                                                            #
# github.com/ciniro                                                           #
# V0.01 em 09/09/2019                                                         #
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

generate_initial_model <- function(size) {

  pm <- list()
  
  for (i in 1:size) {
    pm[[i]] <- rep(0.5, 2)
  }
  
  return (list(paste(1:size, collapse=","), pm))
}

dec2bin <- function(x){
  tmp <- rev(as.numeric(intToBits(x)))
  id <- seq_len(match(1,tmp,length(tmp))-1)
  tmp[-id]
}

draw_bb <- function(prob_ruler) {
  p <- runif(1)
  ind_partition <- NA
  
  for (i in 1:length(prob_ruler)) {
    if (prob_ruler[i] >= p) {
      ind_partition <- (i - 2)
      break()
    }
  }
  
  return(dec2bin(ind_partition))
}

generate_cromo <- function(best_model, size) {
  
  cromo <- rep(NA, size)

  current_model_desc <- best_model[[1]]
  groups <- strsplit(current_model_desc, ",")[[1]]
  
  current_model <- best_model[[2]]
  
  cont <- 1
  for (percentages in current_model) {
    possibilities <- (0:(length(percentages)-1))
    bb <- dec2bin(sample(possibilities, 1, prob = percentages))
    
    sizepos <- log2(length(possibilities))
    
    if (length(bb) != sizepos) {
      bb <- c(rep(0, sizepos - length(bb)) ,bb)
    }
    
    group <- groups[cont]
    
    if (nchar(group) != 1) {
      positions <- as.numeric(strsplit(groups[1], "-")[[1]])
      
      for (i in 1:length(positions)) {
        cromo[positions[i]] <- bb[i]
      }
      
    } else {
      cromo[as.numeric(as.numeric(group))] <-  bb
    }
    
    cont <- cont + 1
  }
  
  return(cromo)
}

generate_pop <- function(npop, best_model, size) {
  pop <- vector()
  
  for (i in 1:npop) {
    pop <- rbind(pop, generate_cromo(best_model, size))  
  }
  
  return(pop)
}

get_winner <- function(competitors) {
  fits <- apply(competitors, 1, fitfunction)
  ind_winner <- which.max(fits)
  return(competitors[ind_winner,])
}

tournament <- function(pop, npop, k, n) {
  
  winners <- vector()
  nwinners <- ceiling(npop*n)
  
  for (i in 1:nwinners) {
    ind_competitors <- sample(1:npop, k)
    winner <- get_winner(pop[ind_competitors,])
    winners <- rbind(winners, winner)
  }
  
  fits <- apply(winners, 1, fitfunction)
  
  return (list(winners, fits))
}

generate_possible_blocks <- function(nvars, blocks) {
  possibilities <- vector()
  
  #build a table with possibilities binary
  for (i in 0:nvars) {
    bin <- dec2bin(i)
    if (length(bin) < dim(blocks)[2]) {
      zeros <- rep(0, dim(blocks)[2] - length(bin))
      bin <- c(zeros, bin)
    }
    possibilities <- rbind(possibilities, bin)
  }
  
  return(possibilities)
}

count_blocks <- function(possibilities, blocks, nvars) {
  occurrences <- rep(0, (nvars + 1))
  for (i in 1:(nvars + 1)) {
    for (j in 1:dim(blocks)[1]) {
      if (all(possibilities[i,] == blocks[j,]) == TRUE)
        occurrences[i] <- occurrences[i] + 1
    }
  }
  
  return(occurrences)
}

calculate_percentage_blocks <- function(blocks) {
  nvars <- 1
  if (is.null(dim(blocks)[2]) == FALSE)
    nvars <- (2^(dim(blocks)[2])) - 1
  
  percentages <- vector()
  
  if (nvars != 1) {
    #build a table with possibilities binary
    possibilities <- generate_possible_blocks(nvars, blocks)
  
    #count occurrences of each possibility
    occurrences <- count_blocks(possibilities, blocks, nvars)
    
    #convert to percentage
    percentages <- occurrences/sum(occurrences)
  }
  else {
    pones <- (sum(blocks == 1))/length(blocks)
    percentages <- c(1 - pones, pones)
  }
  
  return (percentages)
}

generate_model_desc <- function(indsblock, size) {
  block <- paste(indsblock, collapse = "-")
  nvars <- 1:size

  for (i in nvars) {
    contains <- FALSE
    for (j in indsblock) {
      if (i == j)
        contains <- TRUE
    }
    
    if (contains == FALSE)
      block <- paste(block, ",", i, sep="")
  }
  
  return(block)
}

generate_all_models_desc <- function(size) {
  groups <- unlist(sapply(1:size, function(x) apply(combn(1:size, x), 2, paste, collapse = '')))
  groups <- groups[which(nchar(groups) > 1)]
  
  models_desc <- list()
  cont <- 1
  
  for (tamgroup in 1:size) {
    if (tamgroup == 1) {
      models_desc[[cont]] <- paste(1:size, collapse = ",")
      cont <- cont + 1
    } else {
      partitions <- groups[which(nchar(groups) == tamgroup)]
      for (block in partitions) {
        indsblock <- as.numeric(strsplit(block, "")[[1]])
        models_desc[[cont]] <- generate_model_desc(indsblock, size)
        cont <- cont + 1
      }
    }
  }
  
  return(models_desc)
}

generate_models <- function(size, pop, npop, models_desc) {
  
  models <- list()
  contmodels <- 1
  
  for (model_desc in models_desc) {
    model_percentages <- list()
    contpartitions <- 1
    partitions <- strsplit(model_desc, ",")[[1]]
    
    for (partition in partitions) {
      if (nchar(partition) == 1)
        partition <- as.numeric(partition)
      else
        partition <- as.numeric(strsplit(partition, "-")[[1]])
      
      model_percentages[[contpartitions]] <- calculate_percentage_blocks(pop[,partition])
      contpartitions <- contpartitions + 1
    }
    
    models[[contmodels]] <- model_percentages
    contmodels <- contmodels + 1
  }
    
    return(models)
}

complex.cp <- function(models_percentages, npop) {
  
  models_cp <- vector()
  
  for (model_percentages in models_percentages) {
    summodel <- 0
    
     for (partition_percentages in model_percentages) {
       sumpartition <- 0
       for (percentage in partition_percentages) {
         if (percentage != 0) {
            sumpartition <- sumpartition + ( - percentage * log2(percentage))
         }
       }
       
       summodel <- summodel + sumpartition
     }
    
    models_cp <- c(models_cp, (npop*summodel))
  }

  return (models_cp)
}

complex.cm <- function(models_percentages, npop) {
  models_cm <- vector()
  
  for (model_percentages in models_percentages) {
    summodel <- 0
    for (partition_percentages in model_percentages) {
      summodel <- summodel + ((2^log2(length(partition_percentages))) - 1)
    }
    
    models_cm <- c(models_cm, log2(npop)*summodel)
  }
  
  return(models_cm)
}

complex.cc <- function(models_percentages, npop) {
  cc <- complex.cp(models_percentages, npop) + complex.cm(models_percentages, npop)
  return (cc)
}

select_best_model <- function(models_desc, models_percentages, models_complex) {
  ind_best <- which.min(models_complex)
  best_model_desc <- models_desc[[ind_best]]
  best_model_percentages <- models_percentages[[ind_best]]
  
  return(list(best_model_desc, best_model_percentages))
}

#-----------------------------------------------------------------------------#
#                               ALGORITHM                                     #
#                                                                             #
#                     CINIRO APARECIDO LEITE NAMETALA                         #
#-----------------------------------------------------------------------------#

n.tournament <- 0.5 #percent
k.tournament <- 3
size <- 3           #size of solution (max 9)
maxgen <- 50
npop <- 100
fitfunction <- f.onemax

#inicializes a variable independet partition model
best_model <- generate_initial_model(size)

#generate initial population
pop <- generate_pop(npop, best_model, size)

#generate description of possible models
models_desc <- generate_all_models_desc(size)

#generate initial best cromosoma
best <- generate_cromo(best_model)
fit_best <- fitfunction(best)
best <- list(best, fit_best)

#storage best fitness by iteration
bestfitness <- rep(NA, maxgen)

for (iter in 1:maxgen) {
  cat('\014')
  cat(iter)
  
  #let them compete
  winners <- tournament(pop, npop, k.tournament, n.tournament)

  #generate percentage of models based on winners
  models_percentages <- generate_models(size, winners[[1]], dim(winners[[1]])[1], models_desc)
  
  #calculate complex of models
  models_complex <- complex.cc(models_percentages, npop)
  
  #select best model
  best_model <- select_best_model(models_desc, models_percentages, models_complex)

  #update population based on the less complex model
  pop <- generate_pop(npop, best_model, size)
  
  #evaluate best against winner
  ind_best_winners <- which.max(winners[[2]])

  if (winners[[2]][ind_best_winners] > fit_best) {
    best <- list(winners[[1]][ind_best_winners,], winners[[2]][ind_best_winners])
  }
  
  bestfitness[iter] <- best[[2]]

  #convergence verification
  if  (iter > 5)
    if (all.equal(bestfitness[iter], bestfitness[iter - 1], bestfitness[iter - 2], bestfitness[iter - 3], bestfitness[iter - 4]))
      break()
}

plot(bestfitness, type="l", xlab="Iterations", 
     ylab="Best fitness", 
     main="Convergence", panel.first=grid(), lwd=2)

print("STOP ITERATION:")
print(iter)
print("BEST SOLUTION:")
print(best[[1]])
print("FIT BEST SOLUTION:")
print(best[[2]])
print("BEST MODEL BLOCKS:")
print(best_model[[1]])
print("BEST MODEL PERCENTAGES:")
#print(best_model[[2]])

