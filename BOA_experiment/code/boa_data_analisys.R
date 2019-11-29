#-----------------------------------------------------------------------------#
#                                  BOA                                        #
#                             DATA ANALISYS                                   #
#                                                                             #
# Ciniro Ap. Leite Nametala                                                   #
# ciniro@gmail.com                                                            #
# github.com/ciniro                                                           #
# V0.01 em 31/10/2019                                                         #
#-----------------------------------------------------------------------------#

#prepare environment
rm(list=ls())
cat('\014')

#libraries
library(dplyr)
library(tidyr)
library(tibble)

#LOAD DATA: 50_5 -------------------------------------------------------
dfgeral_50_5 <- read.csv("results/output_BB_50_5/geral_experiment.csv",
                    sep=";",
                    header = TRUE)

#mean of iterations to stop
dfsumm_50_5 <- dfgeral_50_5 %>%
            group_by(alg, sco) %>%
            summarise(iter = mean(iter), 
                      mean = mean(mean), 
                      median = median(median), 
                      sd = sd(sd)) %>%
            arrange(alg, sco)

#counts how many times the optimum has been hit
dfoptim_50_5 <- dfgeral_50_5 %>%
            group_by(alg, sco) %>%
            count(best) %>%
            filter(best == 50) %>%
            arrange(alg, sco)

#aggregate
dfsumm_50_5 <- add_column(dfsumm_50_5, n = dfoptim_50_5$n)
rm(dfoptim_50_5)

#save file
#export df summarized
write.table(as.data.frame(dfsumm_50_5),
            file = "results/output_BB_50_5/summarized.csv", 
            row.names=FALSE, 
            col.names=TRUE, 
            sep = ";")

#------------------------------------------------------------------------
#LOAD DATA: 50_10 -------------------------------------------------------
dfgeral_50_10 <- read.csv("results/output_BB_50_10/geral_experiment.csv",
                         sep=";",
                         header = TRUE)

#mean of iterations to stop
dfsumm_50_10 <- dfgeral_50_10 %>%
  group_by(alg, sco) %>%
  summarise(iter = mean(iter), 
            mean = mean(mean), 
            median = median(median), 
            sd = sd(sd)) %>%
  arrange(alg, sco)

#counts how many times the optimum has been hit
dfoptim_50_10 <- dfgeral_50_10 %>%
  group_by(alg, sco) %>%
  count(best) %>%
  filter(best == 50) %>%
  arrange(alg, sco)

#aggregate
dfsumm_50_10 <- add_column(dfsumm_50_10, n = dfoptim_50_10$n)
rm(dfoptim_50_10)

#save file
#export df summarized
write.table(as.data.frame(dfsumm_50_10),
            file = "results/output_BB_50_10/summarized.csv", 
            row.names=FALSE, 
            col.names=TRUE, 
            sep = ";")

#------------------------------------------------------------------------
#LOAD DATA: 100_10 -------------------------------------------------------
dfgeral_100_10 <- read.csv("results/output_BB_100_10/geral_experiment.csv",
                         sep=";",
                         header = TRUE)

#mean of iterations to stop
dfsumm_100_10 <- dfgeral_100_10 %>%
  group_by(alg, sco) %>%
  summarise(iter = mean(iter), 
            mean = mean(mean), 
            median = median(median), 
            sd = sd(sd)) %>%
  arrange(alg, sco)

#counts how many times the optimum has been hit
dfoptim_100_10 <- dfgeral_100_10 %>%
  group_by(alg, sco) %>%
  count(best) %>%
  filter(best == 100) %>%
  arrange(alg, sco)

#aggregate
dfsumm_100_10 <- add_column(dfsumm_100_10, n = dfoptim_100_10$n)
rm(dfoptim_100_10)

#save file
#export df summarized
write.table(as.data.frame(dfsumm_100_10),
            file = "results/output_BB_100_10/summarized.csv", 
            row.names=FALSE, 
            col.names=TRUE, 
            sep = ";")

#------------------------------------------------------------------------
#LOAD DATA: 100_25 -------------------------------------------------------
dfgeral_100_25 <- read.csv("results/output_BB_100_25/geral_experiment.csv",
                         sep=";",
                         header = TRUE)

#mean of iterations to stop
dfsumm_100_25 <- dfgeral_100_25 %>%
  group_by(alg, sco) %>%
  summarise(iter = mean(iter), 
            mean = mean(mean), 
            median = median(median), 
            sd = sd(sd)) %>%
  arrange(alg, sco)

#counts how many times the optimum has been hit
dfoptim_100_25 <- dfgeral_100_25 %>%
  group_by(alg, sco) %>%
  count(best) %>%
  filter(best == 100) %>%
  arrange(alg, sco)

#aggregate
dfsumm_100_25 <- add_column(dfsumm_100_25, n = dfoptim_100_25$n)
rm(dfoptim_100_25)

#save file
#export df summarized
write.table(as.data.frame(dfsumm_100_25),
            file = "results/output_BB_100_25/summarized.csv", 
            row.names=FALSE, 
            col.names=TRUE, 
            sep = ";")

#------------------------------------------------------------------------
#LOAD DATA: 300_30 -------------------------------------------------------
dfgeral_300_30 <- read.csv("results/output_BB_300_30/geral_experiment.csv",
                         sep=";",
                         header = TRUE)

#mean of iterations to stop
dfsumm_300_30 <- dfgeral_300_30 %>%
  group_by(alg, sco) %>%
  summarise(iter = mean(iter), 
            mean = mean(mean), 
            median = median(median), 
            sd = sd(sd)) %>%
  arrange(alg, sco)

#counts how many times the optimum has been hit
dfoptim_300_30 <- dfgeral_300_30 %>%
  group_by(alg, sco) %>%
  count(best) %>%
  filter(best == 300) %>%
  arrange(alg, sco)

#aggregate
dfsumm_300_30 <- add_column(dfsumm_300_30, n = dfoptim_300_30$n)
rm(dfoptim_300_30)

#save file
#export df summarized
write.table(as.data.frame(dfsumm_300_30),
            file = "results/output_BB_300_30/summarized.csv", 
            row.names=FALSE, 
            col.names=TRUE, 
            sep = ";")

#------------------------------------------------------------------------
rotulos <- c("HC-AIC","HC-BDS","HC-BIC","HC-K2","TABU-AIC","TABU-BDS","TABU-BIC","TABU-K2")
cor <-  c("chartreuse1","yellow","gold1","lavenderblush1")

#_50_5
#OBTENCAO DO OTIMO VS VELOCIDADE DE CONVERGENCIA

dfsumm <- dfsumm_50_5

d <- data.frame(row.names=toupper(rotulos), 
                iter = dfsumm$iter, 
                ns = dfsumm$n)
d <- do.call(rbind, d)

bp <- barplot(d, beside = TRUE, 
              legend.text = c("Iteração média de parada", "Número de vezes que o ótimo foi atingido"), 
              args.legend = list(x = "topleft", 
                                 bty="n", 
                                 cex=0.9), 
              names.arg = toupper(rotulos), 
              xlab="Algoritmo - Métrica de score",
              cex.names = 0.8,
              col=c("grey47","gray86"),
              ylim=c(0,25))

dados <- rep(NA, length(dfsumm$n)*2)
cont <- 1
for (i in 1:length(dados)) {
  dados[cont] <- dfsumm$iter[i]
  dados[cont + 1] <- dfsumm$n[i]
  cont <- cont + 2
}

text(x = bp, 
     y = as.numeric(as.vector(t(dados))), 
     label = round(as.numeric(as.vector(t(dados))),1), 
     pos = 3, 
     cex = 0.9)

mtext("20 execuções - Indivíduo de tamanho 50 - 5 Blocos de tamanho 10", side = 3)

#_50_10
#OBTENCAO DO OTIMO VS VELOCIDADE DE CONVERGENCIA

dfsumm <- dfsumm_50_10

d <- data.frame(row.names=toupper(rotulos), 
                iter = dfsumm$iter, 
                ns = dfsumm$n)
d <- do.call(rbind, d)

bp <- barplot(d, beside = TRUE, 
              legend.text = c("Iteração média de parada", "Número de vezes que o ótimo foi atingido"), 
              args.legend = list(x = "topleft", 
                                 bty="n", 
                                 cex=0.9), 
              names.arg = toupper(rotulos), 
              xlab="Algoritmo - Métrica de score",
              cex.names = 0.8,
              col=c("grey47","gray86"),
              ylim=c(0,25))

dados <- rep(NA, length(dfsumm$n)*2)
cont <- 1
for (i in 1:length(dados)) {
  dados[cont] <- dfsumm$iter[i]
  dados[cont + 1] <- dfsumm$n[i]
  cont <- cont + 2
}

text(x = bp, 
     y = as.numeric(as.vector(t(dados))), 
     label = round(as.numeric(as.vector(t(dados))),1), 
     pos = 3, 
     cex = 0.9)

mtext("20 execuções - Indivíduo de tamanho 50 - 10 Blocos de tamanho 5", side = 3)

#_100_10
#OBTENCAO DO OTIMO VS VELOCIDADE DE CONVERGENCIA

dfsumm <- dfsumm_100_10

d <- data.frame(row.names=toupper(rotulos), 
                iter = dfsumm$iter, 
                ns = dfsumm$n)
d <- do.call(rbind, d)

bp <- barplot(d, beside = TRUE, 
              legend.text = c("Iteração média de parada", "Número de vezes que o ótimo foi atingido"), 
              args.legend = list(x = "topleft", 
                                 bty="n", 
                                 cex=0.9), 
              names.arg = toupper(rotulos), 
              xlab="Algoritmo - Métrica de score",
              cex.names = 0.8,
              col=c("grey47","gray86"),
              ylim=c(0,25))

dados <- rep(NA, length(dfsumm$n)*2)
cont <- 1
for (i in 1:length(dados)) {
  dados[cont] <- dfsumm$iter[i]
  dados[cont + 1] <- dfsumm$n[i]
  cont <- cont + 2
}

text(x = bp, 
     y = as.numeric(as.vector(t(dados))), 
     label = round(as.numeric(as.vector(t(dados))),1), 
     pos = 3, 
     cex = 0.9)

mtext("20 execuções - Indivíduo de tamanho 100 - 10 Blocos de tamanho 10", side = 3)

#_50_5
#OBTENCAO DO OTIMO VS VELOCIDADE DE CONVERGENCIA

dfsumm <- dfsumm_100_25

d <- data.frame(row.names=toupper(rotulos), 
                iter = dfsumm$iter, 
                ns = dfsumm$n)
d <- do.call(rbind, d)

bp <- barplot(d, beside = TRUE, 
              legend.text = c("Iteração média de parada", "Número de vezes que o ótimo foi atingido"), 
              args.legend = list(x = "topleft", 
                                 bty="n", 
                                 cex=0.9), 
              names.arg = toupper(rotulos), 
              xlab="Algoritmo - Métrica de score",
              cex.names = 0.8,
              col=c("grey47","gray86"),
              ylim=c(0,25))

dados <- rep(NA, length(dfsumm$n)*2)
cont <- 1
for (i in 1:length(dados)) {
  dados[cont] <- dfsumm$iter[i]
  dados[cont + 1] <- dfsumm$n[i]
  cont <- cont + 2
}

text(x = bp, 
     y = as.numeric(as.vector(t(dados))), 
     label = round(as.numeric(as.vector(t(dados))),1), 
     pos = 3, 
     cex = 0.9)

mtext("20 execuções - Indivíduo de tamanho 100 - 25 Blocos de tamanho 4", side = 3)

#_300_30
#OBTENCAO DO OTIMO VS VELOCIDADE DE CONVERGENCIA

dfsumm <- dfsumm_300_30

d <- data.frame(row.names=toupper(rotulos), 
                iter = dfsumm$iter, 
                ns = dfsumm$n)
d <- do.call(rbind, d)

bp <- barplot(d, beside = TRUE, 
              legend.text = c("Iteração média de parada", "Número de vezes que o ótimo foi atingido"), 
              args.legend = list(x = "topleft", 
                                 bty="n", 
                                 cex=0.9), 
              names.arg = toupper(rotulos), 
              xlab="Algoritmo - Métrica de score",
              cex.names = 0.8,
              col=c("grey47","gray86"),
              ylim=c(0,25))

dados <- rep(NA, length(dfsumm$n)*2)
cont <- 1
for (i in 1:length(dados)) {
  dados[cont] <- dfsumm$iter[i]
  dados[cont + 1] <- dfsumm$n[i]
  cont <- cont + 2
}

text(x = bp, 
     y = as.numeric(as.vector(t(dados))), 
     label = round(as.numeric(as.vector(t(dados))),1), 
     pos = 3, 
     cex = 0.9)

mtext("20 execuções - Indivíduo de tamanho 300 - 30 Blocos de tamanho 10", side = 3)
#-------------------------------------------------------------------
#DESVIO PADRAO MEDIO


dados <- round(dfsumm$sd, 2)
bp <- barplot(dados,
              ylab="Desvio Padrão Médio",
              xlab="Algoritmo - Métrica de score",
              ylim=c(0,5),
              names.arg = toupper(rotulos),
              cex.names = 0.9,
              col=cor)
text(x=bp, y=dados, label=dados, cex=0.8, pos=3)


#DESVIO PADRAO MEDIO


dados <- round(dfsumm$sd, 2)
bp <- barplot(dados,
              ylab="Desvio Padrão Médio",
              xlab="Algoritmo - Métrica de score",
              ylim=c(0,5),
              names.arg = toupper(rotulos),
              cex.names = 0.9,
              col=cor)
text(x=bp, y=dados, label=dados, cex=0.8, pos=3)