#-----------------------------------------------------------------------------#
#                                  BOA                                        #
#                            BINARY PROBLEMS                                  #
#                                                                             #
# Ciniro Ap. Leite Nametala                                                   #
# ciniro@gmail.com                                                            #
# github.com/ciniro                                                           #
# V0.01 em 30/10/2019                                                         #
#-----------------------------------------------------------------------------#

#prepare environment
rm(list=ls())
cat('\014')

#libraries
library(bnlearn)

#external functions
source("boa_functions.R")

#==============================================================================
#learning algorithms
algos <- c("hc","tabu")

#network scores
scores <- c("k2","aic","bic","bds")

#parameters of experiments
nexp <- 1
#==============================================================================

#parameters of selection
n.tournament <- 1 #percent
k.tournament <- 10

#size of cromossomo must be divisible by qtdebb
size <- 12
qtdebb <<- 3

#parameters of convergence
maxgen <- 20
npop <- 100
fitfunction <- f.trapk
tol <- 0

#parameters of plots
plotlast <- TRUE
exportConfigFile <- TRUE

if (exportConfigFile == TRUE) {
  nameExperiment <- "BB_12_3"
  
  title1 <- "BAYESIAN OPTIMIZATION ALGORITHM"
  numberexp <- paste("NAME OF EXPERIMENT: ", nameExperiment, sep="")
  tExps <- paste("TOTAL OF EXPERIMENTS: ", nexp, sep="")
  sep1 <- "================================================="
  title2 <- "PARAMETERS"
  tNTournament <- paste("N TOURNAMENT: ", n.tournament, sep="")
  tKTournament <- paste("K TOURNAMENT: ", k.tournament, sep="")
  tSizeBlock <- paste("SIZE BLOCK: ", size, sep="")
  tQtdeBB <- paste("QT OF BBs: ", qtdebb, sep="")
  tMaxGen <- paste("MAX GENERATIONS: ", maxgen, sep="")
  tNPop <- paste("SIZE OF POPULATION: ", npop, sep="")
  tFitFunction <- "FIT FUNCTION: TRAP K"
  tTolerance <- paste("TOLERANCE: ", tol, sep="")
  sep2 <- "================================================="
  title3 <- "RESULTS"
  
  pathfile <- paste("results/output_", nameExperiment, sep="")
  mkdirs(pathfile)
  
  namefile_config <- "data_experiment.txt"
  allpath_config <- paste(pathfile,"/",namefile_config,sep="")
  
  namefile_geral <- "geral_experiment.csv"
  allpath_geral <- paste(pathfile,"/",namefile_geral,sep="")
  
  writeLines(c(title1,
               numberexp,
               tExps,
               sep1,
               title2,
               tNTournament,
               tKTournament,
               tSizeBlock,
               tQtdeBB,
               tMaxGen,
               tNPop,
               tFitFunction,
               tTolerance,
               sep2,
               title3), 
               allpath_config)
}

dfgeral <- NULL

#==============================================================================
for (n in 1:nexp) {
  namefile_trees <- paste("trees", n, ".txt", sep="")
  namefile_inds <- paste("ind", n, ".csv", sep="")
  namefile_stats <- paste("stats", n, ".csv", sep="")
  
  #inicializes randomly a BN
  best_bn_zero <- generate_initial_model(size)
  
  #generate initial population
  pop_zero <- rbn(best_bn_zero$model, n = npop)
  
  #get initial best cromosoma
  bestcromo_zero <- get_best(pop_zero)
  
  for (alg in algos) {
    for (sco in scores) {
      
      print(paste("Exec: ",alg,"-",sco))
      
      path_trees <- paste(pathfile,"/",alg,"_",sco,sep="")
      path_inds <- paste(pathfile,"/",alg,"_",sco,sep="")
      path_stats <- paste(pathfile,"/",alg,"_",sco,sep="")
      
      mkdirs(path_trees)
      mkdirs(path_inds)
      mkdirs(path_stats)
      
      allpath_trees <- paste(path_trees,"/",namefile_trees,sep="")
      allpath_inds <- paste(path_inds,"/",namefile_inds,sep="")
      allpath_stats <- paste(path_stats,"/",namefile_stats,sep="")

      #inicializes randomly a BN
      best_bn <- best_bn_zero
      
      #generate initial population
      pop <- pop_zero
      
      #get initial best cromosoma
      bestcromo <- bestcromo_zero
      
      #storage stats populations by iteration
      dfstats <- matrix(data=NA,nrow=maxgen,ncol=6)
      colnames(dfstats) <- c("mean","median","sd","best","worst","score")
      
      dfinds <- NULL
      
      meanfitness <- rep(NA, maxgen)
      bestfitness <- rep(NA, maxgen)
      
      for (iter in 1:maxgen) {
        print(iter)
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
        best_bn <- build_bn(winners[[1]], size, best_bn$bn, alg, sco)
      
        #update population based on new BN
        pop <- rbn(best_bn$model, n = npop)
        
        #evaluate best against winner
        newbestcromo <- get_best(pop)
        
        if (newbestcromo$fit_best > bestcromo$fit_best) {
          bestcromo <- newbestcromo
        }
        
        lsstats <- get_stats_pop(pop)
        
        dfstats[iter,"mean"] <- round(lsstats$meanpop,2)
        dfstats[iter,"median"] <- round(lsstats$medianpop,2)
        dfstats[iter,"sd"] <- round(lsstats$sdpop,2)
        dfstats[iter,"best"] <- round(lsstats$bestpop,2)
        dfstats[iter,"worst"] <- round(lsstats$worstpop,2)
        dfstats[iter,"score"] <- round(score(best_bn$bn, pop, type = sco),2)
        
        meanfitness[iter] <- lsstats$meanpop
        bestfitness[iter] <- newbestcromo$fit_best

        dfinds <- rbind(dfinds,newbestcromo$best)
      }
      
      write.table(dfstats[1:(iter-1),], file=allpath_stats, row.names=FALSE, col.names=TRUE, sep = ";")
      write.table(dfinds, file=allpath_inds, row.names=FALSE, col.names=TRUE, sep = ";")
      
      write(paste(alg,"-",sco),
            file=allpath_config,
            append=TRUE)
      write(paste("BEST FIT:",newbestcromo$fit_best),
            file=allpath_config,
            append=TRUE)
      write(paste("STOP:",iter),
            file=allpath_config,
            append=TRUE)
      write("--------------------------------------------",
            file=allpath_config,
            append=TRUE)
      
      results <- c(n, 
                   alg, 
                   sco, 
                   newbestcromo$fit_best,
                   iter,
                   round(lsstats$meanpop,2),
                   round(lsstats$medianpop,2),
                   round(lsstats$sdpop,2))
      
      dfgeral <- rbind(dfgeral, results)
                   
      #==============================================================================
      if (plotlast == TRUE) {
        #MEAN FITNESS CONVERGENCE
        plot(meanfitness, type="l", xlab="Iterations", 
             ylab="Mean fitness", 
             main=paste(alg, "-", sco, ": Convergence - Mean fitness", sep=""), panel.first=grid(), lwd=2)
        
        abline(v=iter-1, col="blue", lty=2)
        
        #BEST FITNESS CONVERGENCE
        # plot(bestfitness, type="l", xlab="Iterations", 
        #      ylab="Best fitness", 
        #      main="Convergence - Best fitness", panel.first=grid(), lwd=2)
        # 
        # abline(v=iter-1, col="blue", lty=2)
        
        # print("BEST MODEL BLOCKS:")
        # print(best_bn$bn)
        # print(best_bn$model)
        graph::plot(best_bn$bn, main=paste(alg, "-", sco, sep=""))
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
  }
    }
  }
}

colnames(dfgeral) <- c("exp","alg","sco","best","iter","mean","median","sd")
write.table(dfgeral, file=allpath_geral, row.names=FALSE, col.names=TRUE, sep = ";")
