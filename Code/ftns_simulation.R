#https://monolix.lixoft.com/monolix-api/#ListFunctions


# Scenario 1: S=2, M=10, N_test=100 (# of data: 1,000)
# T: Random sampling from the distribution of tau
# S: Interval between measurements
# M: The number of test measurements
# E: Random sampling from the Gaussian distribution based on the error model
# N_test, N_cros: 각각의 분석에서 사용할 true curve 개수
# K: The number of iterations for one scenario
# threshold: Detection limit of viral load



# Model
Covfun <- function(pars){
  
  derivs<-function(time, y, pars){
    
    beta <- as.numeric(pars["beta"])*(10^(-5))
    r <- as.numeric(pars["gamma"])
    delta <- as.numeric(pars["delta"])
    
    with(as.list(c(pars, y)),{
      df <- -beta*f*V
      dV <- gamma*f*V-delta*V
      
      return(list(c(df,dV)))
    })
  }
  y<-c(f=1, V=0.01)
  
  times<-c(seq(0,30,0.01))
  out<-lsoda(y=y, parms=pars, times=times, func=derivs, rtol=0.00004, atol=0.00000000000001)
  out2<-data.frame(aV=(log10(out[,3])))
  return(out2)
}


## Compute CI of population parameters
parameter_confidence_intervals <- function(population_parameters) {
  # This function computes the CIs for population parameters using standard errors derived from the Fisher information matrix (FIM).
  # 'population_parameters' must be a DataFrame object with the population parameters derived from Monolix.
  
  # A log-normal distribution is assumed for the population parameters.
  # The FIM and variance-covariance matrix are calculated on the transformed normally distributed parameters.
  # *_pop is the median of the estimated population distribution (i.e. exp(mu)), while omega_* is the deviation from the population typical value.
  # Standard errors are back-transformed to Gaussian space for calculation of the confidence interval before re-transformation into log-space.
  population_parameters <- population_parameters %>%
    rownames_to_column(., var = "parameter") %>%
    as.data.frame() %>%
    mutate(se_trans = ifelse((endsWith(parameter, "_pop") | startsWith(parameter, "omega_")),
                             sqrt(log((1 + sqrt(1 + (2 * se_sa / value)^2)) / 2)), NA)) %>%
    mutate(lower = ifelse((endsWith(parameter, "_pop") | startsWith(parameter, "omega_")),
                          exp(log(value) + qnorm(alpha / 2) * se_trans),
                          value + qnorm(alpha / 2) * se_sa),
           upper = ifelse((endsWith(parameter, "_pop") | startsWith(parameter, "omega_")),
                          exp(log(value) + qnorm(1 - (alpha / 2)) * se_trans),
                          value + qnorm(1 - (alpha / 2)) * se_sa)) %>%
    dplyr::select(parameter, value, lower, upper, CV, se_sa, rse_sa) %>%
    column_to_rownames(., var = "parameter") %>%
    suppressWarnings()
  return(population_parameters)
}

# Predicted population parameters 
get_predicted_parameters <- function(population_parameters) {
  # extract the fixed effects
  gamma_pop <- ifelse(is.na(population_parameters["gamma_pop", "value"]), 0, population_parameters["gamma_pop", "value"])
  beta_pop <- ifelse(is.na(population_parameters["beta1_pop", "value"]), 0, population_parameters["beta1_pop", "value"] * (10^(-5)))
  delta_pop <- ifelse(is.na(population_parameters["delta_pop", "value"]), 0, population_parameters["delta_pop", "value"])
  tau_pop <- ifelse(is.na(population_parameters["tau_pop", "value"]), 0, population_parameters["tau_pop", "value"])
  
  gamma_lower <- ifelse(is.na(population_parameters["gamma_pop", "lower"]), 0, population_parameters["gamma_pop", "lower"])
  beta_lower <- ifelse(is.na(population_parameters["beta1_pop", "lower"]), 0, population_parameters["beta1_pop", "lower"] * (10^(-5)))
  delta_lower <- ifelse(is.na(population_parameters["delta_pop", "lower"]), 0, population_parameters["delta_pop", "lower"])
  tau_lower <- ifelse(is.na(population_parameters["tau_pop", "lower"]), 0, population_parameters["tau_pop", "lower"])
  
  gamma_upper <- ifelse(is.na(population_parameters["gamma_pop", "upper"]), 0, population_parameters["gamma_pop", "upper"])
  beta_upper <- ifelse(is.na(population_parameters["beta1_pop", "upper"]), 0, population_parameters["beta1_pop", "upper"] * (10^(-5)))
  delta_upper <- ifelse(is.na(population_parameters["delta_pop", "upper"]), 0, population_parameters["delta_pop", "upper"])
  tau_upper <- ifelse(is.na(population_parameters["tau_pop", "upper"]), 0, population_parameters["tau_pop", "upper"])
  
  # Construct data frames
  predicted_parameters <- data.frame(gamma = gamma_pop,
                                     beta = beta_pop,
                                     delta = delta_pop,
                                     tau = tau_pop)
  predicted_lower <- data.frame(gamma_lower = gamma_lower,
                                beta_lower = beta_lower,
                                delta_lower = delta_lower,
                                tau_lower = tau_lower)
  predicted_upper <- data.frame(gamma_upper = gamma_upper,
                                beta_upper = beta_upper,
                                delta_upper = delta_upper,
                                tau_upper = tau_upper)
  
  # Combine data frames
  predicted_params <- cbind(predicted_parameters, predicted_lower, predicted_upper)
  
  # Rename columns
  colnames(predicted_params) <- c("gamma", "beta", "delta", "tau", 
                                  "gamma_lower", "beta_lower", "delta_lower", "tau_lower",
                                  "gamma_upper", "beta_upper", "delta_upper", "tau_upper")
  
  rownames(predicted_params) <- NULL
  return(predicted_params)
}



# Generate parameters for N ground truth VL (Not randomly sampled tau)
simulate_population_parameters <- function(predicted_parameters, population_parameters, n) {
  # This function randomly samples parameter sets from the estimated parameter distribution.
  
  # Check input.
  if (is.vector(predicted_parameters)) {
    predicted_parameters <- data.frame(t(predicted_parameters))
  } else if (is.matrix(predicted_parameters)) {
    predicted_parameters <- data.frame(predicted_parameters)
  }
  
  if (!all(c("gamma", "beta", "delta", "tau") %in% colnames(predicted_parameters))) {
    stop("Data is of incorrect format. Please check input parameters.", call. = FALSE)
  } else if (nrow(predicted_parameters) > 1) {
    stop("This function only supports a single set of parameters for simulation.", call. = FALSE)
  }
  
  # Obtain the standard deviation of the random effects (omega_*).
  omega_gamma <- population_parameters["omega_gamma", "value"]
  omega_beta <- population_parameters["omega_beta", "value"]
  omega_delta <- population_parameters["omega_delta", "value"]
  omega_tau <- population_parameters["omega_tau", "value"]
  
  # Random sampling of population parameters
  simulated_parameters <- data.frame(gamma = rlnorm(n, meanlog = log(predicted_parameters$gamma), sdlog = omega_gamma),
                                     beta = rlnorm(n, meanlog = log(predicted_parameters$beta), sdlog = omega_beta),
                                     delta = rlnorm(n, meanlog = log(predicted_parameters$delta), sdlog = omega_delta),
                                     tau = rgamma(1, shape=4.09, scale=1.08))
  
  # Original estimated distribution
  # rlnorm(n, meanlog = log(predicted_parameters$tau), sdlog = omega_tau)
  
  # [ref] Delta
  # T.tau <- rgamma(1, shape=4.09, scale=1.08)
  # [ref] Omicron
  # T.tau <- rgamma(1, shape=4.07, scale=1.12)
  colnames(simulated_parameters) <- c("gamma", "beta", "delta", "tau") #, "tau"
  rownames(simulated_parameters) <- NULL
  return(simulated_parameters)
}

# Error model of VL
error_dist <- function(n=1){
  return(rnorm(n, mean=0, sd=1.3**2))
}


# generate ground truth VLs (step1)
make_indiVL <- function(IDnum, true_vl_par_df){
  indipars <- c(beta=true_vl_par_df[IDnum,"beta"], gamma=true_vl_par_df[IDnum,"gamma"],
                delta=true_vl_par_df[IDnum,"delta"])
  fitted <- Covfun(indipars)
  
  return(fitted)
}


# generate observed VLs (step2) - One trial (Simulation)
make_obs_vl <- function(true_vl_df, IDnum, Tnum, Snum, Mnum) {
  # True VL values
  truevl <- true_vl_df[c("time", paste0("VL_",IDnum))]
  colnames(truevl) <- c("time", "VL")
  
  # Generate obs. data
  obsvl <- truevl$VL[truevl$time %in% seq(Tnum, by=Snum, length.out=Mnum)]
  obsvl <- obsvl + error_dist(Mnum) # add error
  obsvl <- cbind(obsvl, error_dist(Mnum))
  return(obsvl)
}

sim_paras <- function(params_data, S, M, N_partici, iter){
  
  # population parameters
  pop_params <- parameter_confidence_intervals(params_data)
  pred_pop_params <- get_predicted_parameters(pop_params)
  
  # parameters for 10,000 truth VL
  sim_params <- simulate_population_parameters(pred_pop_params, pop_params, N)
  true_vl_par <- data.frame(ID=1:N, beta=sim_params$beta,
                            gamma=sim_params$gamma,
                            delta=sim_params$delta) #,tau=sim_params$tau
  
  ## For longitudinal
  true_vl <- data.frame(time=times)
  # Random sampling from true viral load parameters
  true_vl_par <- true_vl_par[sample(1:N, N_partici, replace = F),]
  true_vl_par$ID <- 1:N_partici
  for (i in 1:N_partici) {
    indiVL <- make_indiVL(i,true_vl_par)
    true_vl[[paste0("VL_", i)]] <- indiVL$aV 
  }
  
  # Only integer times
  true_vl <- true_vl[true_vl$time %in% seq(0, 30, by = 1), ]
  
  # Generate observed VLs
  obs_vl <- data.frame(ID=rep(1:N_partici,each=M), time=as.numeric(NA), VL=as.numeric(NA))
  for (i in 1:N_partici){
    repeat {
      # Delta
      T.tau <- rgamma(1, shape=4.09, scale=1.08) # rate = 0.92
      # Omicron
      # T.tau <- rgamma(1, shape=3.93, scale=0.92)      
      T.tau <- ceiling(T.tau)
      if ((T.tau + S*M) <= max_time) break
    }
    if (M==1){
      obs_vl[i,"time"] <- T.tau
      obs_vl_value <- make_obs_vl(true_vl, i, T.tau, 1, 1)
      obs_vl[i,"VL"] <- obs_vl_value[1,1]
      obs_vl[i,"error"] <- obs_vl_value[1,2]
    } else{
      obs_vl[(1+M*(i-1)):(M*(i)),"time"] <- seq(T.tau, by=S, length.out=M)
      obs_vl_value <- make_obs_vl(true_vl, i, T.tau, S, M)
      obs_vl[(1+M*(i-1)):(M*(i)),"VL"] <- obs_vl_value[,1]
      obs_vl[(1+M*(i-1)):(M*(i)),"error"] <- obs_vl_value[,2]
    }
  }

  ## Estimate the population parameters for longitudinal or cross-sectional data.
  if (M==1){
    obs_vl['censored'] <- as.numeric(obs_vl$VL <= threshold)
    obs_vl[obs_vl['censored']==1,'VL'] <- threshold
    obs_vl['ID2'] <- 1
  } else{
    obs_vl['censored'] <- as.numeric(obs_vl$VL <= threshold)
    obs_vl[obs_vl['censored']==1,'VL'] <- threshold
  }
  write.csv(obs_vl, paste0(data_path, N_test, "obsdata_N_",N_partici,"_M_",M,"_S_",S,"_",iter,".csv"))
  
  return(obs_vl)
}


#### 2. Run monolix
## Author: Hoong Kai Chua


################################## FIT DATA IN MONOLIX ##################################

monolix_startup <- function(monolixPath = NULL) {
  # This function sets up the Monolix software in R and ensures that items have been loaded correctly.
  # See documentation: https://monolix.lixoft.com/monolix-api/lixoftconnectors_installation/
  # 'monolixPath' refers to the path to the installation directory of the Monolix suite.
  
  # Check if the lixoftConnectors package is installed.
  if (system.file(package = "lixoftConnectors") == "") {
    stop("The lixoftConnectors package needs to be installed.", call. = FALSE)
  } else {
    require(lixoftConnectors)
  }
  
  # Set the Monolix path if not provided.
  if (is.null(monolixPath)) {
    if (.Platform$OS.type == "windows") {
      monolixPath <- "C:/Users/user/AppData/Roaming/Microsoft/Windows/Start Menu/Programs/Lixoft/MonolixSuite2023R1" #jihyeon's labtop
    } else if (.Platform$OS.type == "unix") {
      monolixPath <- "/Applications/MonolixSuite2023R1.app/"
    }
  }
  
  # Initialize the Monolix API.
  state <- getLixoftConnectorsState(quietly = TRUE)
  if (is.null(state)) {
    initializeLixoftConnectors(software = "monolix", path = monolixPath, force = TRUE)
    state <- getLixoftConnectorsState(quietly = TRUE)
    if (is.null(state)) {
      stop("The lixoftConnectors package failed to load.", call. = FALSE)
    }
  }
  
  # Check if the Monolix software and lixsoftConnectors package versions are identical.
  monolix_version <- state$version
  monolix_version_major <- as.numeric(sub("([0-9]+).*$", "\\1", monolix_version))
  lixoftConnectors_version <- packageVersion("lixoftConnectors")
  lixoftConnectors_version_major <- as.numeric(sub("([0-9]+).*$", "\\1", lixoftConnectors_version))
  
  if (monolix_version_major != lixoftConnectors_version_major) {
    stop(paste0("The major version number for the Monolix software and lixoftConnectors package must be identical.",
                "\nMonolix software version -> ", monolix_version,
                "\nRsmlx package version -> ", lixoftConnectors_version), call. = FALSE)
  }
  return(invisible(TRUE))
}

run_monolix <- function(dataFile, modelFile, headerTypes, covariateModel, initialParameters, projectFile = "monolix_project.mlxtran",
                        dataDirectory = NULL, modelDirectory = NULL, projectDirectory = NULL, monolixPath = NULL, verbose = TRUE,
                        runEBE = TRUE, runCDS = TRUE, runFIM = TRUE, runLL = TRUE, runPlots = TRUE) {
  # This function automatically fits the data in Monolix with pre-set settings.
  # See documentation: https://monolix.lixoft.com/monolix-api/#ListFunctions
  # 'dataFile' refers to the absolute or relative path to the data file.
  # 'modelFile' refers to the absolute or relative path to the model file.
  # 'headerTypes' refers to a collection of header types found in the same order as the dataFile.
  # 'covariateModel' is a list of comma separated pairs: parameterName = c(covariateName = TRUE/FALSE)
  # 'initialParameters' is a list-like object containing 'initialValue' and 'method' for each parameter.
  # 'projectFile' is the name of the Monolix project file which will be saved.
  # 'dataDirectory' is the main directory of the data file (if relative path was specified).
  # 'modelDirectory' is the main directory of the model file (if relative path was specified).
  # 'projectDirectory' is the main directory of the project file (if relative path was specified).
  # 'monolixPath' refers to the path to the installation directory of the Monolix suite.
  # 'verbose' is a Boolean indicating whether to print messages.
  # 'runEBE' is a Boolean indicating whether to run empirical Bayes' estimation.
  # 'runCDS' is a Boolean indicating whether conditional distribution sampling should be performed.
  # 'runFIM' is a Boolean indicating whether the Fisher information matrix should be computed.
  # 'runLL' is a Boolean indicating whether log likelihood estimates should be computed.
  # 'runPlots' is a Boolean indicating whether diagnostic plots should be computed and exported.
  
  # Run Monolix startup.
  state <- monolix_startup(monolixPath = monolixPath)
  if (state != TRUE) {
    stop("Monolix run has been terminated due to failure in startup.", call. = FALSE)
  }
  
  # Error handling on data directory.
  dataDirectory <- ifelse(is.null(dataDirectory), getwd(), dataDirectory)
  modelDirectory <- ifelse(is.null(modelDirectory), getwd(), modelDirectory)
  projectDirectory <- ifelse(is.null(projectDirectory), getwd(), projectDirectory)
  dataFile <- ifelse(endsWith(dataDirectory, "/"), paste0(dataDirectory, dataFile),  paste0(dataDirectory, "/", dataFile))
  modelFile <- ifelse(endsWith(modelDirectory, "/"), paste0(modelDirectory, modelFile), paste0(modelDirectory, "/", modelFile))
  projectFile <- ifelse(endsWith(projectDirectory, "/"), paste0(projectDirectory, projectFile), paste0(projectDirectory, "/", projectFile))
  
  # Create a new Monolix project.
  newProject(data = list(dataFile = dataFile,
                         headerTypes = headerTypes,
                         observationTypes = list("VL" = "continuous")),
             modelFile = modelFile)
  
  # Save the project.
  saveProject(projectFile = projectFile)
  
  # Set the structural model.
  setStructuralModel(modelFile = modelFile)
  
  # Set the observational model and the error model type to be used.
  # Observations are normally distributed and the error is constant.
  setObservationDistribution(VL = "normal")
  setErrorModel(VL = "constant")
  
  # Set parameters to a log-normal distribution to ensure positiveness.
  # Allow random effects for all parameters.
  setIndividualParameterDistribution(gamma = "logNormal", beta1 = "logNormal", delta = "logNormal", tau = "logNormal")
  setIndividualParameterVariability(gamma = TRUE, beta1 = TRUE, delta = TRUE, tau = TRUE)
  
  # Transform the continuous age to set the reference to the weighted mean age of the dataset.
  # addContinuousTransformedCovariate(tAge = "Age - 42")
  
  # Set the covariate model.
  setCovariateModel(covariateModel)
  
  # Set the initial population parameters and estimation methods.
  setPopulationParameterInformation(initialParameters)
  saveProject(projectFile = projectFile)
  
  # Set higher maximum iterations as convergence may not be reached with larger sample sizes.
  setPopulationParameterEstimationSettings(method = "saem",
                                           exploratoryAutoStop = TRUE, nbExploratoryIterations = 1000,
                                           smoothingAutoStop = TRUE, nbSmoothingIterations = 200)
  
  # The probability of a given individual parameter cannot be directly calculated (no closed form).
  # Sampling from the conditional distribution has to be done by the Metropolis-Hastings algorithm.
  # Set the minimum number of simulated parameters for computing the prediction intervals.
  setConditionalDistributionSamplingSettings(nbMinIterations = 100, nbSimulatedParameters = 100)
  
  # Indicator for verbosity.
  if (verbose == TRUE) {
    setConsoleMode(mode = "complete")
  } else {
    setConsoleMode(mode = "none")
  }
  
  # Set the tasks and run scenario. By default, SAEM must be run first.
  scenario <- getScenario()
  scenario$tasks <- c(populationParameterEstimation = TRUE, conditionalModeEstimation = runEBE,
                      conditionalDistributionSampling = runCDS, standardErrorEstimation = runFIM,
                      logLikelihoodEstimation = runLL, plots = runPlots)
  setScenario(scenario)
  runScenario()
  
  # Save the project again.
  saveProject(projectFile = projectFile)
  
  # Compute and export charts data.
  if (runPlots == TRUE) {
    computeChartsData()
    exportChartDataSet(type = "vpc")
    exportChartDataSet(type = "indfits")
  }
}


run_monolix_script <- function(N_partici, M, S, iter) {
  # Set the path to the data and model files
  dataFile <- paste0(data_path, N_test, "obsdata_N_",N_partici,"_M_",M,"_S_",S,"_",iter,".csv")
  if(M==1){
    headerTypes <- c("ignore", "ignore", "time", "observation", "ignore", "cens", "id")
  } else{
    headerTypes <- c("ignore", "id", "time", "observation", "ignore", "cens")
  }
  projectFile = paste0( N_test,"obsdata_N_",N_partici,"_M_",M,"_S_",S,"_",iter,".mlxtran")
  modelFile <- paste0(model_path,"Model.txt")
  
  # covariate model
  covariateModel <- list(
    "CL" = c("WT" = FALSE, "AGE" = FALSE),
    "V" = c("WT" = FALSE, "AGE" = FALSE)
  )
  
  # initial value of parameters
  initialParameters <- list(
    "gamma" = list(initialValue = 20.0, method = "estimate"),
    "beta1" = list(initialValue = 0.01, method = "estimate"),
    "delta" = list(initialValue = 1.0, method = "estimate"),
    "tau" = list(initialValue = 5.0, method = "estimate")
  )
  
  # Run Monolix
  run_monolix(
    dataFile = dataFile,
    modelFile = modelFile,
    headerTypes = headerTypes,
    covariateModel = covariateModel,
    initialParameters = initialParameters,
    projectFile = projectFile,
    dataDirectory = data_path,
    modelDirectory = model_path,
    projectDirectory = project_path,
    monolixPath = monolixPath,
    verbose = TRUE,
    runEBE = TRUE,
    runCDS = TRUE,
    runFIM = TRUE,
    runLL = TRUE,
    runPlots = TRUE
  )
  
  cat("iteration=",i,", Monolix run completed!")
}

#### 3. Population curve fitting
# Model
Covfun2 <- function(pars){
  
  beta <- as.numeric(pars["beta"])
  r <- as.numeric(pars["gamma"])
  delta <- as.numeric(pars["delta"])
  
  derivs<-function(time, y, pars){
    with(as.list(c(pars, y)),{
      df <- -beta*f*V
      dV <- gamma*f*V-delta*V
      
      return(list(c(df,dV)))
    })
  }
  y<-c(f=1, V=0.01)
  
  times<-c(seq(0,30,0.01))
  out<-lsoda(y=y, parms=pars, times=times, func=derivs, rtol=0.00004, atol=0.00000000000001)
  # out2<-cbind(time=out[,1],aV=((log10(out[,3]))))
  out2<-data.frame(aV=(log10(out[,3])))
  return(out2)
}


# ## Compute CI of population parameters
# parameter_confidence_intervals <- function(population_parameters) {
#   population_parameters <- population_parameters %>%
#     rownames_to_column(., var = "parameter") %>%
#     as.data.frame() %>%
#     mutate(se_trans = ifelse((endsWith(parameter, "_pop") | startsWith(parameter, "omega_")),
#                              sqrt(log((1 + sqrt(1 + (2 * se_sa / value)^2)) / 2)), NA)) %>%
#     mutate(lower = ifelse((endsWith(parameter, "_pop") | startsWith(parameter, "omega_")),
#                           exp(log(value) + qnorm(alpha / 2) * se_trans),
#                           value + qnorm(alpha / 2) * se_sa),
#            upper = ifelse((endsWith(parameter, "_pop") | startsWith(parameter, "omega_")),
#                           exp(log(value) + qnorm(1 - (alpha / 2)) * se_trans),
#                           value + qnorm(1 - (alpha / 2)) * se_sa)) %>%
#     dplyr::select(parameter, value, lower, upper, CV, se_sa, rse_sa) %>%
#     column_to_rownames(., var = "parameter") %>%
#     suppressWarnings()
#   return(population_parameters)
# }




# generate ground truth VLs (step1)
esti_vl <- function(params){
  
  infectiontime <- params[['tau']]
  # newtimes <- c(seq(-infectiontime, 30-infectiontime, 0.01))
  newtimes <- c(seq(0, max_time, 0.01))
  
  sR_pars <- c(beta = params[['beta']], gamma = params[['gamma']], delta = params[['delta']])
  sR_pars_lower <- c(beta = params[['beta_lower']], gamma = params[['gamma_lower']], delta = params[['delta_lower']])
  sR_pars_upper <- c(beta = params[['beta_upper']], gamma = params[['gamma_upper']], delta = params[['delta_upper']])
  
  # sR <- Covfun2(sR_pars)
  # sR_lower <- Covfun2(sR_pars_lower)
  # sR_upper <- Covfun2(sR_pars_upper)
  
  fitted <- Covfun2(sR_pars)
  infectiontime <- params[['tau_lower']]
  fitted_lower <- Covfun2(sR_pars_lower)
  infectiontime <- params[['tau_upper']]
  fitted_upper <- Covfun2(sR_pars_upper)
  
  fitted1 <- data.frame(time = newtimes, aV = fitted$aV, aV_lower = fitted_lower$aV, aV_upper = fitted_upper$aV)
  
  return(fitted1)
}


plot_vl <- function(fit, title, curvecolor){
  p <- ggplot(fit, aes(x = time, y = aV)) +
    geom_line(color = curvecolor, linewidth = 1.5) +
    labs(title = title, x = "Time (day)", y = "Viral load (log10 copies/mL)") +
    theme_bw() +
    scale_y_continuous(limits=c(0, 9)) +
    scale_x_continuous(limits=c(-2, 25)) +
    theme(plot.title = element_text(hjust = 0.5, size = 15),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "none")
  return(p)
}


#### 4. Comparison
# get_duration <- function(fit, infectionlimit){
#   duration <- numeric(length(fit))
#   for (i in 1:length(fit)){
#     duration[i] <- max(fit[i]$time[fit[i]$aV > infectionlimit]) - min(fit[i]$time[fit[i]$aV > infectionlimit])
#   }
#   return(duration)
# }
# 
# get_peaksize <- function(fit){
#   peaksize <- numeric(length(fit))
#   for (i in 1:length(fit)){
#     peaksize[i] <- max(fit[i]$aV)
#   }
#   return(peaksize)
# }
# 
# get_peaktime <- function(fit){
#   peaktime <- numeric(length(fit))
#   for (i in 1:length(fit)){
#     peaktime[i] <- fit[i]$time[fit[i]$aV == max(fit[i]$aV)]
#   }
#   return(peaktime)
# }


get_confidence_interval <- function(fit_df, dataname){
  
  if (dataname=="cross-sectional data"){
    plot_col = "#3772ff"
  } else { plot_col = "#ff7f7e"}
  
  combined_data <- bind_rows(fit_df, .id = "simulation_id")

  # confidence_data <- combined_data %>%
  #   group_by(time) %>%
  #   summarise(
  #     mean_aV = mean(aV),
  #     lower_CI = mean_aV - qt(1-alpha/2, df = K - 1) * sd(aV) / sqrt(K),
  #     upper_CI = mean_aV + qt(1-alpha/2, df = K - 1) * sd(aV) / sqrt(K)
  #   )
  
  # 2.5th, 97.5th percentile
  confidence_data <- combined_data %>%
    group_by(time) %>%
    summarise(
      mean_aV = mean(aV),
      lower_CI = quantile(aV, 0.025),
      upper_CI = quantile(aV, 0.975)
    )
  
  confidence_plot <- ggplot(confidence_data, aes(x = time)) +
    geom_line(aes(y = mean_aV, color="Mean Viral load"), linewidth=1.8, linetype = "solid") +
    geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI, fill="Confidence interval"), alpha = 0.2) +
    labs(title = paste0("Confidence Interval of ", dataname), y = "Viral RNA loads", x = "Time") +
    scale_fill_manual(values = c("Confidence interval" = plot_col)) +
    scale_color_manual(values = c("Mean Viral load" = plot_col)) +
    coord_cartesian(ylim = c(0, 10)) +
    scale_x_continuous(limits = c(0, 20)) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=20, face='bold'),
          legend.position = "top", legend.title = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 14),
          legend.text = element_text(size = 14))
  
  return(list(confidence_data, confidence_plot))
}

get_credible_interval <- function(fit_df, dataname){
  
  if (dataname=="cross-sectional data"){
    plot_col = "#3772ff"
  } else { plot_col = "#ff7f7e"}
  
  combined_data <- bind_rows(fit_df, .id = "simulation_id")
  
  credible_data <- combined_data %>%
    group_by(time) %>%
    summarise(
      median_aV = median(aV),
      mean_aV = mean(aV),
      lower_CrI = quantile(aV, 0.025), # 2.5% quantile
      upper_CrI = quantile(aV, 0.975)  # 97.5% quantile
    )
  
  credible_plot <- ggplot(credible_data, aes(x = time)) +
    geom_line(aes(y = mean_aV, color="Mean Viral load"), linewidth = 1.8, linetype = "solid") +
    geom_ribbon(aes(ymin = lower_CrI, ymax = upper_CrI, fill="Credible interval"), alpha = 0.2) +
    labs(title = paste0("Credible Interval of ", dataname), y = "Viral RNA loads", x = "Time") +
    scale_fill_manual(values = c("Credible interval" = plot_col)) +
    scale_color_manual(values = c("Mean Viral load" = plot_col)) +
    coord_cartesian(ylim = c(0, 10)) +
    scale_x_continuous(limits = c(0, 20)) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=20, face='bold'),
          legend.position = "top", legend.title = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 14),
          legend.text = element_text(size = 14))
  
  return(list(credible_data, credible_plot))
}


get_plot_together <- function(df_list){
  
  df_confidence <- df_list[[1]] ;  df_credible <- df_list[[2]]
  
  
  confidence_plot <- ggplot() +
    geom_line(data = df_confidence, aes(x = time, y = mean_aV, color = "Mean Viral load (longitudinal)"), linewidth=1.8) +
    geom_ribbon(data = df_confidence, aes(x = time, ymin = lower_CI, ymax = upper_CI, fill = "Confidence interval (longitudinal)"), alpha = 0.2) +
    labs(title = "Confidence Interval", y = "Viral RNA loads", x = "Time") +
    scale_fill_manual(values = c("Confidence interval" = "#ff7f7e")) +
    scale_color_manual(values = c("Mean Viral load" = "#ff7f7e")) +
    coord_cartesian(ylim = c(0, 10)) +
    scale_x_continuous(limits = c(0, 25)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
          legend.position = c(0.8,0.8), legend.title = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.box.background = element_rect(color = "black"))
  
  credible_plot <- ggplot() +
    geom_line(data = df_credible, aes(x = time, y = mean_aV, color = "Mean Viral load (longitudinal)"), linewidth=1.8) +
    geom_ribbon(data = df_credible, aes(x = time, ymin = lower_CrI, ymax = upper_CrI, fill = "Credible interval (longitudinal)"), alpha = 0.2) +
    labs(title = "Credible Interval", y = "Viral RNA loads", x = "Time") +
    scale_fill_manual(values = c( "Credible interval" = "#ff7f7e")) +
    scale_color_manual(values = c("Mean Viral load" = "#ff7f7e")) +
    coord_cartesian(ylim = c(0, 10)) +
    scale_x_continuous(limits = c(0, 25)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
          legend.position = c(0.8,0.8), legend.title = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.box.background = element_rect(color = "black"))
  
  
  return(list(confidence_plot, credible_plot))
}
