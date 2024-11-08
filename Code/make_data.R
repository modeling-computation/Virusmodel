## Make longitudial & Cross-sectional data using truth ground viral load ##

rm(list=ls())
set.seed(0930)

library(tidyverse)
library(openxlsx)
library(deSolve)
library(dplyr)
library(gridExtra)

environment <- new.env()
source(file = "./functions_make_data.R", local = environment)
attach(environment)

alpha <- 0.05   # statistical significance level
seed <- 1234    # seed for pseudo-random number generation
N <- 10000      # the number of ground truth VL (individuals)
min_time <- 0 ; max_time <- 30 ; time_interval <- 0.01
times<-c(seq(min_time, max_time, time_interval)) # length:3001
# threshold <- (40-49.5932)/(-3.6097) # infectiousness threshold (in log10 form)
threshold <- 1 # censoring 조절을 위해


#### 0. Set the distribution of population parameters.

#### 1, 2. Generate ground truth viral loads and observed data based on the truth VL.
# Scenario 1: T=2, S=3, M=3
# T: Random sampling from the distribution of tau
# S: Interval between measurements
# M: The number of test measurements
# E: Random sampling from the Gaussian distribution based on the error model
# N_long, N_cros: 각각의 분석에서 사용할 true curve 개수

S<-seq(1, 5, by=1)
M<-seq(3, 10, by=1)

# S<-2
# M<-10

N_long <- 100



sim_paras <- function(params_data, S, M, N_long, N_cros){
  
  # population parameters
  pop_params <- parameter_confidence_intervals(params_data)
  pred_pop_params <- get_predicted_parameters(pop_params)
  
  # parameters for 10,000 truth VL
  sim_params <- simulate_population_parameters(pred_pop_params, pop_params, N)
  true_vl_par <- data.frame(ID=1:N, beta=sim_params$beta,
                            gamma=sim_params$gamma,
                            delta=sim_params$delta,
                            tau=sim_params$tau)
  
  ## For longitudinal
  true_vl_long <- data.frame(time=times)
  # Random sampling from true viral load parameters
  true_vl_par_long <- true_vl_par[sample(1:N, N_long, replace = F),]
  true_vl_par_long$ID <- 1:N_long
  for (i in 1:N_long) {
    indiVL_long <- make_indiVL(i,true_vl_par_long)
    true_vl_long[[paste0("VL_", i)]] <- indiVL_long$aV 
  }
  
  # Only integer times
  true_vl_long <- true_vl_long[true_vl_long$time %in% seq(0, 30, by = 1), ]
  
  # Generate observed VLs
  obs_vl_long <- data.frame(ID=rep(1:N_long,each=M), time=as.numeric(NA), VL=as.numeric(NA))
  for (i in 1:N_long){
    repeat {
      T.tau <- rweibull(1, shape=1.28, scale=4.8)
      T.tau <- ceiling(T.tau)
      if ((T.tau + S*M) <= max_time) break
    }
    obs_vl_long[(1+M*(i-1)):(M*(i)),"time"] <- seq(T.tau, by=S, length.out=M)
    obs_vl <- make_obs_vl(true_vl_long, i, T.tau, S, M)
    obs_vl_long[(1+M*(i-1)):(M*(i)),"VL"] <- obs_vl[,1]
    obs_vl_long[(1+M*(i-1)):(M*(i)),"error"] <- obs_vl[,2]
  }
  
  ## For cross-sectional
  true_vl_cros <- data.frame(time=times)
  # Random sampling from true viral load parameters
  true_vl_par_cros <- true_vl_par[sample(1:N, N_cros, replace = F),]
  true_vl_par_cros$ID <- 1:N_cros
  for (i in 1:N_cros) {
    indiVL_cros <- make_indiVL(i,true_vl_par_cros)
    true_vl_cros[[paste0("VL_", i)]] <- indiVL_cros$aV 
  }
  
  # Only integer times
  true_vl_cros <- true_vl_cros[true_vl_cros$time %in% seq(0, 30, by = 1), ]
  
  # Generate observed VLs
  obs_vl_cros <- data.frame(ID=1:N_cros, time=as.numeric(NA), VL=as.numeric(NA))
  
  for (i in 1:N_cros){
    repeat {
      T.tau <- rweibull(1, shape=1.28, scale=4.8)
      T.tau <- ceiling(T.tau)
      if ((T.tau + S*M) <= max_time) break
    }
    obs_vl_cros[i,"time"] <- seq(T.tau, by=S, length.out=M)[1] # = T.tau
    obs_vl <- make_obs_vl(true_vl_cros, i, T.tau, S, M)
    obs_vl_cros[i,"VL"] <- obs_vl[1,1]
    obs_vl_cros[i,"error"] <- obs_vl[1,2]
  }
  
  #### 3. Estimate the population parameters for longitudinal or cross-sectional data.
  
  ## longitudinal data
  obs_vl_long['censored'] <- as.numeric(obs_vl_long$VL < threshold)
  write.csv(obs_vl_long, paste0(getwd(),"/../../Data/longitudinal/longitudinal_M_",M,"_S_",S,".csv"))
  
  ## cross-sectional data
  obs_vl_cros['censored'] <- as.numeric(obs_vl_cros$VL < threshold)
  obs_vl_cros['ID2'] <- 1
  write.csv(obs_vl_cros, paste0("../../Data/crosssectional/crosssectional_M_",M,"_S_",S,".csv"))
  
  return(list(obs_vl_long, obs_vl_cros))
}


# load data
delta_params <- read.table("./Delta_populationParameters.txt", sep = ",", header = TRUE, row.names = 1)
# delta_params <- as.data.frame(t(delta_params))

# simulate data
for (i in 1:length(S)){
  for (j in 1:length(M)){
    N_cros <- N_long*M[j]
    results <- sim_paras(delta_params, S[i], M[j], N_long, N_cros)
    # obs_vl_long[i] <- results[[1]]
    # obs_vl_cros[i] <- results[[2]]
  }
}



#### Figure ####

obs_vl_long <- read.csv("../../Data/longitudinal/longitudinal_M_10_S_2.csv")
obs_vl_cros <- read.csv("../../Data/crosssectional/crosssectional_M_10_S_2.csv")

## longitudinal data
g1 <- ggplot(obs_vl_long, aes(x=time, y=VL)) + geom_point() +
  xlab(xlabel) + ylab(ylabel) +
  scale_x_continuous(breaks=seq(-20,30,by=10),labels = expression(-20,-10,0,10,20,30),limits=c(-20,30)) +
  scale_y_continuous(breaks=seq(1,9,by=2),labels = expression(10^1,10^3,10^5,10^7,10^9),limits=c(0,10)) +
  theme(axis.text = element_text(colour = "black"), axis.ticks = element_line(colour = "black"), axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        legend.position='none', axis.title.y = element_text(size=7,family="Helvetica"), axis.title.x = element_text(size=9,family="Helvetica"))

## cross-sectional data
g2 <- ggplot(obs_vl_cros, aes(x=time, y=VL)) + geom_point() +
  xlab(xlabel) + ylab(ylabel) +
  scale_x_continuous(breaks=seq(-20,30,by=10),labels = expression(-20,-10,0,10,20,30),limits=c(-20,30)) +
  scale_y_continuous(breaks=seq(1,9,by=2),labels = expression(10^1,10^3,10^5,10^7,10^9),limits=c(0,10)) +
  theme(axis.text = element_text(colour = "black"), axis.ticks = element_line(colour = "black"), axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        legend.position='none', axis.title.y = element_text(size=7,family="Helvetica"), axis.title.x = element_text(size=9,family="Helvetica"))

# scale_x_continuous(breaks=seq(min_time, max_time, by=2), limits=c(-2,30)) +
# scale_y_continuous(limits=c(min(obs_vl_cros$VL,obs_vl_long$VL), 10)) +
# ggtitle("Cross-sectional data") +
# theme_minimal() +
# theme(
#   plot.title = element_text(hjust=0.5, face="bold"),
#   axis.text = element_text(colour = "black"),
#   text = element_text(size = 14)
# )

# grid.arrange(g1, g2, ncol=2)
as_ggplot(arrangeGrob(g1, g2, ncol=2))
