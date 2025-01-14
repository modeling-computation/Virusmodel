rm(list=ls())

# install.packages("C:/ProgramData/Lixoft/MonolixSuite2023R1/connectors/lixoftConnectors.tar.gz",
# repos = NULL, type="source", INSTALL_opts ="--no-multiarch")

libraries = c("deSolve", "dplyr","reshape2","openxlsx","tidyverse","lubridate",
              "readxl","ggplot2","viridisLite","colorspace","ggpubr","gridExtra",
              "png","grid","magrittr","scales","RColorBrewer","fBasics",
              "lixoftConnectors", "RJSONIO")
# for(x in libraries) { if (!require(x)) install.packages(x) }
for(x in libraries) { library(x,character.only=TRUE,warn.conflicts=FALSE) }
initializeLixoftConnectors(software = "monolix")

environment <- new.env()
source("./ftns_simulation.R", local = environment)
attach(environment)


# Basic settings
alpha <- 0.05   # statistical significance level
N <- 10000      # the number of ground truth VL (individuals)
min_time <- 0 ; max_time <- 30 ; time_interval <- 0.01
times<-c(seq(min_time, max_time, time_interval)) # length:3001
threshold <- (40-49.5932)/(-3.6097) # threshold: detection limit (in log10 form) = 2.657617
# infectionlimit <- (30-49.5932)/(-3.6097) # infection limit (in log10 form) = 5.4279
infectionlimit <- threshold

# Scenario settings
S_list <- c(3)
M_list <- c(1, 2, 3, 4, 5)
N_test_list <- c(1000) # , 2000, 5000, 10000
variant="delta" ; variant2="Delta"
# variant="omicron" ; variant2="Omicron"
K <- 200

result_df <- as.data.frame(matrix(0, nrow = length(S_list)*length(M_list)*length(N_test_list), ncol = 9))
colnames(result_df) <- c("interval", "indi_test_num", "RMSE", "RMSE_duration", 
                         "RMSE_peaksize", "RMSE_peaktime", "participant_num", "total_test_num", "variant")
result_df$interval <- rep(S_list, each = length(M_list)*length(N_test_list))
result_df$indi_test_num <- rep(M_list, each = length(N_test_list)*length(S_list))
result_df$variant <- rep(variant, length(S_list)*length(M_list)*length(N_test_list))


# Path
data_path <- paste0("../../Data/",variant,"_scenario/")
project_path <- paste0("../Monolix_fitting_Kim/Variant_step2/",variant,"_scenario/")
model_path <- "../Monolix_fitting_Kim/"
monolixPath <- "C:/Users/user/AppData/Roaming/Microsoft/Windows/Start Menu/Programs/Lixoft/MonolixSuite2023R1" #jihyeon's labtop


for (N_test in N_test_list[[1]]){
  for (M in M_list) {
    for (S in S_list){

      N_partici <- ceiling(N_test/M)
      result_df[result_df$interval == S & result_df$indi_test_num == M, "participant_num"] <- N_partici
      result_df[result_df$interval == S & result_df$indi_test_num == M & result_df$participant_num == N_partici, "total_test_num"] <- N_test
      cat("\nS: ", S, "M: ", M, "N_test: ", N_test, "N_participants: ", N_partici, "\n")
      
      obs_vl <- list()

      #### 1. Make data
      #### 2. Run monolix
      groundtrue_params <- read.table(paste0("./",variant2,"_populationParameters.txt"), sep = ",", header = TRUE, row.names = 1) # ground truth
      # monolix_startup(monolixPath)
      # for (i in 1:K){
      #   # Make data
      #   results <- sim_paras(groundtrue_params, S, M, N_partici, iter=i+100)
      #   obs_vl[[i]] <- results
      #   cat(sprintf("\rMake data: %d", i+100))
      # 
      #   # Run monolix
      #   run_monolix_script(N_partici, M, S, iter=i)
      # }



      #### 3. Population curve
      plot <- list() ; rmse_vl <- rep(0, K)
      duration <- rep(0, K) ; peaksize <- rep(0, K)
      peaktime <- rep(0, K) ; coverage <- rep(0, K)
      VL <- list()

      cros_col <- '#3772ff'
      col <- '#ff7f7e'

      min_time <- 0 ; max_time <- 30 ; time_interval <- 0.01
      times<-c(seq(min_time, max_time, time_interval)) # length:3001
      alpha <- 0.05   # statistical significance level
      seed <- 1012    # seed for pseudo-random number generation


      # True population curve
      pop_fit <- read.csv(paste0("../../Result/",variant,"/pop_fit_",variant,".csv"), header = TRUE, row.names = 1)
      pop_duration <- max(pop_fit$time[pop_fit$VL >= infectionlimit]) - min(pop_fit$time[pop_fit$VL >= infectionlimit])
      pop_peaksize <- max(pop_fit$VL)
      pop_peaktime <- pop_fit$time[which.max(pop_fit$VL)]

      for (i in 1:K){
        pop <- read.table(paste0(project_path, N_test, "obsdata_N_",N_partici,"_M_",M,"_S_",S,"_",i,"/populationParameters.txt"), sep = ",", header = TRUE, row.names = 1)
        pop_params <- parameter_confidence_intervals(pop)
        pred_pop_params <- get_predicted_parameters(pop_params)

        fit <- esti_vl(pred_pop_params)

        longplot <- plot_vl(fit, paste0(variant2,"_M_",M,"_S_",S,"",i), col)
        plot[[i]] <- longplot

        # Save viral load values
        VL[[i]] <- fit

        # For comparison
        rmse_vl[i] <- sqrt(mean((fit$aV - pop_fit$VL)^2))
        duration[i] <- max(fit$time[fit$aV >= infectionlimit]) - min(fit$time[fit$aV >= infectionlimit])
        peaksize[i] <- max(fit$aV)
        peaktime[i] <- fit$time[which.max(fit$aV)]

        # coverage[i] <- get_coverage(fit)
      }


      ### 4. Comparison
      # RMSE
      rmse_mean <- mean(rmse_vl)
      cat("RMSE: ", rmse_mean)

      # RMSE for duration of viral shedding
      duration_rmse <- sqrt(mean((duration - pop_duration)^2))
      cat("\nRMSE for duration of viral shedding: ", duration_rmse)

      # RMSE for peak viral load
      peaksize_rmse <- sqrt(mean((peaksize - pop_peaksize)^2))
      cat("\nRMSE for peak viral load: ", peaksize_rmse)

      # RMSE for peak time
      peaktime_rmse <- sqrt(mean((peaktime - pop_peaktime)^2))
      cat("\nRMSE for peak time: ", peaktime_rmse)


      # Save the RMSE
      result_df[result_df$interval == S & result_df$indi_test_num == M & result_df$total_test_num == N_test & result_df$participant_num == N_partici, "RMSE"] <- rmse_mean
      result_df[result_df$interval == S & result_df$indi_test_num == M & result_df$total_test_num == N_test & result_df$participant_num == N_partici, "RMSE_duration"] <- duration_rmse
      result_df[result_df$interval == S & result_df$indi_test_num == M & result_df$total_test_num == N_test & result_df$participant_num == N_partici, "RMSE_peaksize"] <- peaksize_rmse
      result_df[result_df$interval == S & result_df$indi_test_num == M & result_df$total_test_num == N_test & result_df$participant_num == N_partici, "RMSE_peaktime"] <- peaktime_rmse


      # Confidence interval
      confidence <- get_confidence_interval(VL, "longitudinal data")
      confidence_data <- confidence[[1]]
      confidence_plot <- confidence[[2]]

      # Credible interval
      # credible <- get_credible_interval(VL, "longitudinal data")
      # credible_data <- credible[[1]]
      # credible_plot <- credible[[2]]

      # together <- get_plot_together(list(confidence_data, credible_data))
      # together[[1]]
      # together[[2]]
      # ggsave(paste0("../../Figure/",variant,"_scenario/", N_test, "K_", K, "CI_S_",S,"_M_",M,"_N_",N_partici,".png"), plot=together[[1]], width = 12, height = 6)
      # ggsave(paste0("../../Figure/",variant,"_scenario/", N_test, "K_", K, "CrI_S_",S,"_M_",M,"_N_",N_partici,".png"), plot=together[[2]], width = 12, height = 6)


      # Figure
      g1 <- ggplot() +
          geom_line(data = confidence_data, aes(x = time, y = mean_aV, color="Mean Viral load"), linewidth = 1.8, linetype = "solid") +
          geom_ribbon(data = confidence_data, aes(x = time, ymin = lower_CI, ymax = upper_CI, fill="Confidence interval"), alpha = 0.7)
      for (i in 1:K){
        g1 <- g1 + geom_line(data = VL[[i]], aes(x = time, y = aV), color=i, linewidth = 0.3)
      }

      gg1 <- g1 +
        # scale_y_continuous(limits = c(0, 10)) +
        theme_bw() +
        labs(title = paste0("Population curve with 95% CI (test=", M, ", patient=", N_partici, ")"), x = "Days since infection", y = "log10(Viral load) copies/mL") +
        scale_fill_manual(values = c("Confidence interval" = "#3772ff")) +
        scale_color_manual(values = c("Mean Viral load" = "#3772ff")) +
        coord_cartesian(ylim = c(0, 10)) +
        theme(plot.title=element_text(hjust=0.5, size=20, face='bold')) +
        scale_x_continuous(breaks = seq(-5, 25, 5), limits=c(-5, 25)) +
        geom_text(aes(x = 15, y = 7, label = paste("RMSE: ", round(rmse_mean, 2))), size = 5)

      # ggarrange(gg1, ncol = 1, nrow = 1)
      ggsave(paste0("../../Figure/",variant,"_scenario/", N_test, "K_", K, "CI_Population_curve_S_",S,"_M_",M,"_N_",N_partici,".png"), plot=ggarrange(gg1, ncol = 1, nrow = 1), width = 12, height = 6)
    }
  }
}

head(result_df)
write.csv(result_df, "../../Result/delta/RMSE(K_300).csv")
