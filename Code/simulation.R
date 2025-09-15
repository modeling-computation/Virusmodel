rm(list=ls())

install.packages("C:/ProgramData/Lixoft/MonolixSuite2023R1/connectors/lixoftConnectors.tar.gz", 
                 repos = NULL, type="source", INSTALL_opts ="--no-multiarch")

libraries = c("deSolve", "dplyr","reshape2","openxlsx","tidyverse","lubridate",
              "readxl","ggplot2","viridisLite","colorspace","ggpubr","gridExtra",
              "png","grid","magrittr","scales","RColorBrewer","fBasics",
              "lixoftConnectors", "RJSONIO", "boot")

for(x in libraries) { if (!require(x)) install.packages(x) }
for(x in libraries) { library(x,character.only=TRUE,warn.conflicts=FALSE) }
initializeLixoftConnectors(software = "monolix")

environment <- new.env()
source("./ftns_simulation.R", local = environment)
attach(environment)

seed <- 250121    # seed for pseudo-random number generation
set.seed(seed)

# Basic settings
alpha <- 0.05   # statistical significance level
N <- 10000      # the number of ground truth VL (individuals)
min_time <- 0 ; max_time <- 30 ; time_interval <- 0.01
times<-c(seq(min_time, max_time, time_interval)) # length:3001
detectionlimit <- (40-49.5932)/(-3.6097) # threshold: detection limit (in log10 form) = 2.657617
cros_col <- '#3772ff'
col <- '#ff7f7e'

# Scenario settings
S_list <- c(3)
M_list <- c(1, 2, 3, 4, 5, 8)
N_test_list <- c(200, 500, 1000)
variant="delta" ; variant2="Delta"
K <- 100


# Path
data_path <- paste0("../../Data/",variant,"_scenario/")
data_path_ground_truth <- paste0("../../Data/delta_scenario_ground_truth/")
project_path <- paste0("../Monolix_fitting_Kim/Variant_step2/",variant,"_scenario/")
model_path <- "../Monolix_fitting_Kim/"
monolixPath <- "C:/Users/user/AppData/Roaming/Microsoft/Windows/Start Menu/Programs/Lixoft/MonolixSuite2023R1"


# Ground truth parameters
NBA_params <- read.table(paste0("./",variant2,"_populationParameters.txt"), sep = ",", header = TRUE, row.names = 1)
pop_params <- parameter_confidence_intervals(NBA_params)
pred_pop_params <- get_predicted_parameters(pop_params)

# Ground truth viral load 
pop_fit <- read.csv(paste0("../../Result/",variant,"/NBA_Delta_population_fit.csv"), header = TRUE, row.names = 1)

sim_params <- simulate_population_parameters(pred_pop_params, pop_params, N)
sampled_true_vl_par <- data.frame(ID=1:N, beta=sim_params$beta,
                                  gamma=sim_params$gamma,
                                  delta=sim_params$delta)


h=1
for (N_test in N_test_list){
  r=1 # index for saving to data frame
  colnum=1
  
  # RMSE result data frame
  result_df <- as.data.frame(matrix(0, nrow = length(S_list)*length(M_list)*length(N_test), ncol = 8))
  colnames(result_df) <- c("interval", "indi_test_num", "RMSE_VL", "RMSE_VL_lower", "RMSE_VL_upper",
                           "participant_num", "total_test_num", "variant")
  result_df$interval <- rep(S_list, each = length(M_list)*length(N_test))
  result_df$indi_test_num <- rep(M_list, each = length(N_test)*length(S_list))
  result_df$variant <- rep(variant, length(S_list)*length(M_list)*length(N_test))
  
  # Mean result data frame
  result_df2 <- as.data.frame(matrix(0, nrow = length(S_list)*length(M_list)*length(N_test), ncol = 14))
  colnames(result_df2) <- c("interval", "indi_test_num",
                            "duration", "duration_lower", "duration_upper",
                            "peaksize", "peaksize_lower", "peaksize_upper",
                            "peaktime", "peaktime_lower", "peaktime_upper",
                            "participant_num", "total_test_num", "variant")
  result_df2$interval <- rep(S_list, each = length(M_list)*length(N_test))
  result_df2$indi_test_num <- rep(M_list, each = length(N_test)*length(S_list))
  result_df2$variant <- rep(variant, length(S_list)*length(M_list)*length(N_test))
  
  # Coverage score result data frame
  result_coverage <- data.frame(time = seq(min_time, max_time, by=time_interval), M1=0, M2=0, M3=0, M4=0, M5=0, M8=0)
  
  # Parameter result data frame
  result_beta <- data.frame(index = 1:K, beta_U1=0, beta_U2=0, beta_U3=0, beta_U4=0, beta_U5=0, beta_U8=0)
  result_gamma <- data.frame(index = 1:K, gamma_U1=0, gamma_U2=0, gamma_U3=0, gamma_U4=0, gamma_U5=0, gamma_U8=0)
  result_delta <- data.frame(index = 1:K, delta_U1=0, delta_U2=0, delta_U3=0, delta_U4=0, delta_U5=0, delta_U8=0)
  
  result_parameter <- as.data.frame(matrix(0, nrow = length(S_list)*length(M_list)*length(N_test), ncol = 14))
  colnames(result_parameter) <- c("interval", "indi_test_num",
                                  "beta", "beta_lower", "beta_upper",
                                  "delta", "delta_lower", "delta_upper",
                                  "gamma", "gamma_lower", "gamma_upper",
                                  "participant_num", "total_test_num", "variant")
  result_parameter$interval <- rep(S_list, each = length(M_list)*length(N_test))
  result_parameter$indi_test_num <- rep(M_list, each = length(N_test)*length(S_list))
  result_parameter$variant <- rep(variant, length(S_list)*length(M_list)*length(N_test))

  
  cat("Total data points = ", N_test, "\n")

  
  for (M in M_list) {
    for (S in S_list) {
      
      # Set the number of tests
      N_partici <- floor(N_test/M)
      result_df[result_df$interval == S & result_df$indi_test_num == M, "participant_num"] <- N_partici
      result_df[result_df$interval == S & result_df$indi_test_num == M & result_df$participant_num == N_partici, "total_test_num"] <- N_test
      result_df2[result_df2$interval == S & result_df2$indi_test_num == M, "participant_num"] <- N_partici
      result_df2[result_df2$interval == S & result_df2$indi_test_num == M & result_df2$participant_num == N_partici, "total_test_num"] <- N_test
      result_parameter[result_parameter$interval == S & result_parameter$indi_test_num == M, "participant_num"] <- N_partici
      result_parameter[result_parameter$interval == S & result_parameter$indi_test_num == M & result_parameter$participant_num == N_partici, "total_test_num"] <- N_test
      
      # Save the true features of viral load (true metric)
      true_df <- as.data.frame(matrix(0, nrow = 1, ncol = 10))
      colnames(true_df) <- c("interval", "participant_num", "duration", "peaksize", "peaktime", "beta", "delta", "gamma", "total_test_num", "variant")
      true_df$interval <- S ; true_df$total_test_num <- N_test ; true_df$variant <- variant ; true_df$participant_num <- M
      
      #### 1. True viral load data
      true_duration <- max(pop_fit$time[pop_fit$VL >= detectionlimit], na.rm=TRUE) - min(pop_fit$time[pop_fit$VL >= detectionlimit], na.rm=TRUE)
      true_peaksize <- max(pop_fit$VL, na.rm=TRUE)
      true_peaktime <- pop_fit$time[which.max(pop_fit$VL)]
      true_df$duration <- true_duration ; true_df$peaksize <- true_peaksize ; true_df$peaktime <- true_peaktime
      
      #### 2. Synthetic data & estimate K=100 viral trajectories
      monolix_startup(monolixPath)
      for (i in 1:K){
        # 2-1. Make n ground truth
        true_vl_data <- estimate_groundtruth(sampled_true_vl_par, S, M, N_partici, i)
        true_vl_data <- read.csv(paste0(data_path, "../delta_scenario_ground_truth/", N_test, "truedata_N_",N_partici,"_M_",M,"_S_",S,"_realnbr_", i, ".csv"))
        true_vl_data <- true_vl_data[true_vl_data$time %in% seq(0, 30, by = 1), ] # only integer time

        # 2-2. Generate n synthetic data
        results <- sim_paras(true_vl_data, S, M, N_partici, iter=i)
        
        cat(sprintf("\rMake observed data: %d", i))

        ## Estimate K=100 viral trajectories from synthetic data
        run_monolix_script(N_partici, M, S, iter=i)
      }

      
      #### 3. Comparison
      esti_vl <- list() ; esti_vl_PI <- list()
      rmse_vl <- rep(0, K) ; rmse_duration <- rep(0, K) ; rmse_peaksize <- rep(0, K) ; rmse_peaktime <- rep(0, K)
      esti_duration <- rep(0, K) ; esti_peaksize <- rep(0, K) ; esti_peaktime <- rep(0, K)
      
      for (i in 1:K){
        sim_pop <- read.table(paste0(project_path, N_test, "obsdata_N_",N_partici,"_M_",M,"_S_",S,"_",i,"/populationParameters.txt"), sep = ",", header = TRUE, row.names = 1)
        fit_list <- draw_trajectory(sim_pop)
        
        cat(sprintf("\rEstimated %dth trajectory (m=%d, s=%d, u=%d)", i, N_test, S, M))
        
        fit <- fit_list[[1]]
        fit_PI <- fit_list[[2]] #lower, upper
        
        # Save estimated viral load
        esti_vl[[i]] <- fit$VL
        esti_vl_PI[[i]] <- fit_PI
        
        # RMSE
        rmse_vl[i] <- sqrt(mean((pop_fit$VL - esti_vl[[i]])^2, na.rm=TRUE))
        
        # Metrics
        esti_duration[i] <- max(fit$time[fit$VL >= detectionlimit], na.rm=TRUE) - min(fit$time[fit$VL >= detectionlimit], na.rm=TRUE)
        esti_peaksize[i] <- max(fit$VL, na.rm=TRUE)
        esti_peaktime[i] <- fit$time[which.max(fit$VL)]

        # Parameters
        result_beta[result_beta$index == i, paste0("beta_U", M)] <- sim_pop['beta1_pop','value']
        result_gamma[result_gamma$index == i, paste0("gamma_U", M)] <- sim_pop['gamma_pop','value']
        result_delta[result_delta$index == i, paste0("delta_U", M)] <- sim_pop['delta_pop','value']
      }
      
      ### 4. Comparison
      # Calculate the median and 95% PI
      esti_duration_ci <- vl_95CI(esti_duration)
      esti_peaksize_ci <- vl_95CI(esti_peaksize) # median, lower, upper
      esti_peaktime_ci <- vl_95CI(esti_peaktime)
      
      # Calculate the RMSE and 95% CI
      rmse_ci <- vl_95CI(rmse_vl)

      # Calculate coverage score of estimated viral load
      coverage_vl_ratio <- numeric(length(pop_fit$time))
      for (t in 1:length(pop_fit$time)) {
        count_in_interval <- sum(sapply(esti_vl_PI, function(df) df$lower[t] <= pop_fit$VL[t] & pop_fit$VL[t] <= df$upper[t]))
        coverage_vl_ratio[t] <- count_in_interval / K
      }
      result_coverage[,r+1] <- coverage_vl_ratio
      r <- r+1
      
      ## Calculate the 95% CI for parameters
      beta_ci <- vl_25_75QT(result_beta[,paste0("beta_U", M)] )
      gamma_ci <- vl_25_75QT(result_gamma[,paste0("gamma_U", M)])
      delta_ci <- vl_25_75QT(result_delta[,paste0("delta_U", M)])
      
      # RMSE result
      result_df[result_df$interval == S & result_df$indi_test_num == M & result_df$total_test_num == N_test & result_df$participant_num == N_partici, "RMSE_VL"] <- rmse_ci[1]
      result_df[result_df$interval == S & result_df$indi_test_num == M & result_df$total_test_num == N_test & result_df$participant_num == N_partici, "RMSE_VL_lower"] <- rmse_ci[2]
      result_df[result_df$interval == S & result_df$indi_test_num == M & result_df$total_test_num == N_test & result_df$participant_num == N_partici, "RMSE_VL_upper"] <- rmse_ci[3]

      # Mean result
      result_df2[result_df2$interval == S & result_df2$indi_test_num == M & result_df2$total_test_num == N_test & result_df2$participant_num == N_partici, "duration"] <- esti_duration_ci[1]
      result_df2[result_df2$interval == S & result_df2$indi_test_num == M & result_df2$total_test_num == N_test & result_df2$participant_num == N_partici, "duration_lower"] <- esti_duration_ci[2]
      result_df2[result_df2$interval == S & result_df2$indi_test_num == M & result_df2$total_test_num == N_test & result_df2$participant_num == N_partici, "duration_upper"] <- esti_duration_ci[3]
      
      result_df2[result_df2$interval == S & result_df2$indi_test_num == M & result_df2$total_test_num == N_test & result_df2$participant_num == N_partici, "peaksize"] <- esti_peaksize_ci[1]
      result_df2[result_df2$interval == S & result_df2$indi_test_num == M & result_df2$total_test_num == N_test & result_df2$participant_num == N_partici, "peaksize_lower"] <- esti_peaksize_ci[2]
      result_df2[result_df2$interval == S & result_df2$indi_test_num == M & result_df2$total_test_num == N_test & result_df2$participant_num == N_partici, "peaksize_upper"] <- esti_peaksize_ci[3]
      
      result_df2[result_df2$interval == S & result_df2$indi_test_num == M & result_df2$total_test_num == N_test & result_df2$participant_num == N_partici, "peaktime"] <- esti_peaktime_ci[1]
      result_df2[result_df2$interval == S & result_df2$indi_test_num == M & result_df2$total_test_num == N_test & result_df2$participant_num == N_partici, "peaktime_lower"] <- esti_peaktime_ci[2]
      result_df2[result_df2$interval == S & result_df2$indi_test_num == M & result_df2$total_test_num == N_test & result_df2$participant_num == N_partici, "peaktime_upper"] <- esti_peaktime_ci[3]
      
      # Parameter result
      result_parameter[result_parameter$interval == S & result_parameter$indi_test_num == M & result_parameter$total_test_num == N_test & result_parameter$participant_num == N_partici, "beta"] <- beta_ci[1]
      result_parameter[result_parameter$interval == S & result_parameter$indi_test_num == M & result_parameter$total_test_num == N_test & result_parameter$participant_num == N_partici, "beta_lower"] <- beta_ci[2]
      result_parameter[result_parameter$interval == S & result_parameter$indi_test_num == M & result_parameter$total_test_num == N_test & result_parameter$participant_num == N_partici, "beta_upper"] <- beta_ci[3]
      
      result_parameter[result_parameter$interval == S & result_parameter$indi_test_num == M & result_parameter$total_test_num == N_test & result_parameter$participant_num == N_partici, "delta"] <- delta_ci[1]
      result_parameter[result_parameter$interval == S & result_parameter$indi_test_num == M & result_parameter$total_test_num == N_test & result_parameter$participant_num == N_partici, "delta_lower"] <- delta_ci[2]
      result_parameter[result_parameter$interval == S & result_parameter$indi_test_num == M & result_parameter$total_test_num == N_test & result_parameter$participant_num == N_partici, "delta_upper"] <- delta_ci[3]
      
      result_parameter[result_parameter$interval == S & result_parameter$indi_test_num == M & result_parameter$total_test_num == N_test & result_parameter$participant_num == N_partici, "gamma"] <- gamma_ci[1]
      result_parameter[result_parameter$interval == S & result_parameter$indi_test_num == M & result_parameter$total_test_num == N_test & result_parameter$participant_num == N_partici, "gamma_lower"] <- gamma_ci[2]
      result_parameter[result_parameter$interval == S & result_parameter$indi_test_num == M & result_parameter$total_test_num == N_test & result_parameter$participant_num == N_partici, "gamma_upper"] <- gamma_ci[3]
    }
    

    colnum <- colnum+1
  }
  
  write.csv(result_df, paste0("Rmse_",N_test,"_S",S,".csv"))
  write.csv(result_df2, paste0("Metric_",N_test,"_S",S,".csv"))
  write.csv(true_df, paste0("True_",N_test,"_S",S,".csv"))
  write.csv(result_coverage, paste0("Coveragescore_",N_test,"_S",S,".csv"))
  write.csv(result_parameter, paste0("Parameter_",N_test,"_S",S,".csv"))
  
  h <- h+1
}