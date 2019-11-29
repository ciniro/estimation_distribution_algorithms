#-----------------------------------------------------------------------------#
#                            BOA - FUNCTIONS                                  #
#                            BINARY PROBLEMS                                  #
#                                                                             #
# Ciniro Ap. Leite Nametala                                                   #
# ciniro@gmail.com                                                            #
# github.com/ciniro                                                           #
# V0.01 em 23/10/2019                                                         #
#-----------------------------------------------------------------------------#

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

get_stats_pop <- function(pop) {
  mpop <- data.matrix(pop)-1
  fits <- apply(mpop, 1, fitfunction)
  
  lsstats <- list(meanpop = mean(fits),
                  medianpop = median(fits),
                  sdpop = sd(fits),
                  bestpop = max(fits),
                  worstpop = min(fits))
  
  return(lsstats)
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

build_bn <- function(dfpop, size, previous_bn, alg, sco) {
  dfpop <- apply(dfpop, 2, as.character)
  #include levels (at least one 0 e one 1)
  dfpop <- rbind(dfpop, rep("0", size), rep("1", size))
  dfpop <- as.data.frame(dfpop)
  xi <- paste("x", 1:size, sep="")
  colnames(dfpop) <- xi
  
  bn <- NULL
  if (alg == "hc")
    bn <- hc(dfpop, score = sco, start = previous_bn, optimized = TRUE)
  else if (alg == "tabu")
    bn <- tabu(dfpop, score = sco, start = previous_bn, optimized = TRUE)
  
  model = bn.fit(bn, dfpop, replace.unidentifiable = TRUE)
  
  return(list(model = model, bn = bn))
}

readFile = function(filepath) {
  con = file(filepath, "r")
  value <- NULL
  while ( TRUE ) {
    line = readLines(con, n = 1, warn = FALSE)
    if ( length(line) == 0 ) {
      break
    }
    value <- line
  }
  
  close(con)
  return(value)
}

mkdirs <- function(fp) {
  if(!file.exists(fp)) {
    mkdirs(dirname(fp))
    dir.create(fp)
  }
} 
