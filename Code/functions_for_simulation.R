#https://monolix.lixoft.com/monolix-api/#ListFunctions


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
                                     tau = rgamma(1, shape=4.09, scale=1.08)) # from reference, mean=4.43, sd=2.19
  colnames(simulated_parameters) <- c("gamma", "beta", "delta", "tau")
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


###
estimate_groundtruth <- function(true_vl_par_10000, S, M, N_partici, iter){
  
  true_vl <- data.frame(time=times)
  
  # Random sampling from true viral load parameters
  true_vl_par <- true_vl_par_10000[sample(1:N, N_partici, replace = F),]
  true_vl_par$ID <- 1:N_partici
  
  for (i in 1:N_partici) {
    indiVL <- make_indiVL(i,true_vl_par)
    true_vl[[paste0("VL_", i)]] <- indiVL$aV 
  }
  
  write.csv(true_vl, paste0(data_path_ground_truth, N_test, "truedata_N_",N_partici,"_M_",M,"_S_",S,"_realnbr_", iter, ".csv"))
  # Only integer times
  true_vl <- true_vl[true_vl$time %in% seq(0, 30, by = 1), ]
  
  # Fit population VL of true_vl
  true_vl_long <- true_vl %>%
    pivot_longer(cols = starts_with("VL_"), names_to = "ID", values_to = "VL") %>%
    mutate(ID = as.integer(sub("VL_", "", ID)),
           censored = 0)
  true_vl_long['censored'] <- as.numeric(true_vl_long$VL <= detectionlimit)
  true_vl_long[true_vl_long['censored']==1,'VL'] <- detectionlimit
  write.csv(true_vl_long, paste0(data_path_ground_truth, N_test, "truedata_N_",N_partici,"_M_",M,"_S_",S,"_", iter, ".csv"))
  run_monolix_script_ground_truth(N_partici, M, S) # estimated trajectory of ground truth VL
  
  return(true_vl)
}



##

sim_paras <- function(true_vl, S, M, N_partici, iter){
  
  # Generate observed VLs
  obs_vl <- data.frame(ID=rep(1:N_partici,each=M), time=as.numeric(NA), VL=as.numeric(NA))
  for (i in 1:N_partici){
    repeat {
      # Delta
      T.tau <- rgamma(1, shape=4.09, scale=1.08) # rate = 0.92
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
    obs_vl['censored'] <- as.numeric(obs_vl$VL <= detectionlimit)
    obs_vl[obs_vl['censored']==1,'VL'] <- detectionlimit
    obs_vl['ID2'] <- 1
  } else{
    obs_vl['censored'] <- as.numeric(obs_vl$VL <= detectionlimit)
    obs_vl[obs_vl['censored']==1,'VL'] <- detectionlimit
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
  
  # Set the structural model.
  setStructuralModel(modelFile = modelFile)
  
  # Set the observational model and the error model type to be used.
  # Observations are normally distributed and the error is constant.
  setObservationDistribution(VL_ = "normal")
  setErrorModel(VL_ = "constant")
  
  # Set parameters to a log-normal distribution to ensure positiveness.
  # Allow random effects for all parameters.
  setIndividualParameterDistribution(gamma = "logNormal", beta1 = "logNormal", delta = "logNormal", tau = "logNormal")
  setIndividualParameterVariability(gamma = TRUE, beta1 = TRUE, delta = TRUE, tau = TRUE)
  
  # Set the initial population parameters and estimation methods.
  setPopulationParameterInformation(gamma_pop = list(initialValue = 20.0),
                                    beta1_pop = list(initialValue = 0.01),
                                    delta_pop = list(initialValue = 1.0),
                                    tau_pop = list(initialValue = 5.0))
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
  
  # Save the project.
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
  projectFile = paste0(N_test,"obsdata_N_",N_partici,"_M_",M,"_S_",S,"_",iter,".mlxtran")
  modelFile <- paste0(model_path,"Model_simulation.txt")
  
  # initial value of parameters
  initialParameters <- c(
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
  out2<-data.frame(aV=(log10(out[,3])))
  return(out2)
}


# generate ground truth VLs (step1)
esti_vl <- function(params){
  
  newtimes <- c(seq(0, max_time, 0.01))
  
  sR_pars <- c(beta = params[['beta']], gamma = params[['gamma']], delta = params[['delta']])
  sR_pars_lower <- c(beta = params[['beta_lower']], gamma = params[['gamma_lower']], delta = params[['delta_lower']])
  sR_pars_upper <- c(beta = params[['beta_upper']], gamma = params[['gamma_upper']], delta = params[['delta_upper']])
  
  
  fitted <- Covfun2(sR_pars)
  fitted_lower <- Covfun2(sR_pars_lower)
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
vl_onestd <- function(data){
  mean_val <- mean(data)
  se <- sd(data)
  lower_bound <- mean_val - se
  upper_bound <- mean_val + se
  return(c(mean_val, lower_bound, upper_bound))
}

# 25%, 75% quantile
vl_25_75QT <- function(data) {
  lower_bound <- quantile(data, alpha/2)
  upper_bound <- quantile(data, 1-alpha/2)
  return(c(median(data), lower_bound, upper_bound))
}


draw_trajectory <- function(pop_params) {
  parms <- get_parameters(population_parameters = pop_params)
  
  # The population fit.
  pop_fit <- fit_model(parameters = parms)
  
  # Sample from the parameter distribution.
  sim_parms <- simulate_population_parameters(predicted_parameters = parms, population_parameters = pop_params, n = 100)
  
  # Get prediction intervals.
  sim_fit <- fit_simulations(simulated_parameters = sim_parms) # fitted VL
  intervals <- simulated_prediction_intervals(predictions = sim_fit)
  return(list(pop_fit, intervals))
}


get_parameters <- function(population_parameters) {
  
  # Extract the fixed effects (*_pop).
  gamma_pop <- population_parameters["gamma_pop", "value"]
  beta_pop <- population_parameters["beta1_pop", "value"]
  delta_pop <- population_parameters["delta_pop", "value"]
  tau_pop <- population_parameters["tau_pop", "value"]
  
  # Compute the predicted population parameters based on the covariate model.
  predicted_parameters <- data.frame(gamma = gamma_pop ,
                                     beta = beta_pop *  10**(-5),
                                     delta = delta_pop,
                                     tau = tau_pop)
  colnames(predicted_parameters) <- c("gamma", "beta", "delta", "tau")
  rownames(predicted_parameters) <- NULL
  return(predicted_parameters)
}


fit_model <- function(parameters) {
  # Check input.
  if (is.vector(parameters)) {
    parameters <- data.frame(t(parameters))
  } else if (is.matrix(parameters)) {
    parameters <- data.frame(parameters)
  }
  
  if (!all(c("gamma", "beta", "delta", "tau") %in% colnames(parameters))) {
    stop("Data is of incorrect format. Please check input parameters.", call. = FALSE)
  } else if (nrow(parameters) > 1) {
    stop("This function only supports model fitting for a single parameter set.", call. = FALSE)
  }
  
  # Define the parameters.
  gamma <- as.numeric(parameters$gamma) # maximum rate constant for viral replication
  beta <- as.numeric(parameters$beta)   # rate constant for virus infection
  delta <- as.numeric(parameters$delta) # death rate of infected cells
  tau <- as.numeric(parameters$tau)     # time interval from infection to diagnosis
  y <- c(f = 1, V = 0.01)               # initial state values for the ODE system
  
  # Times at which estimates for y are desired.
  # The time scale here is days after infection (i.e. negative times do not make sense). # time: 감염 이후의 날
  times <- seq(from = 0, to = 30, by = 0.01)
  
  # Derivatives of y with respect to time.
  derivatives <- function(times, y, parms) {
    with(as.list(c(y, parms)), {
      df <- -beta * f * V             # fraction of uninfected cells
      dV <- gamma * f * V - delta * V # amount of virus (copies/mL)
      return(list(c(df, dV)))
    })
  }
  
  output <- lsoda(y = y, times = times, func = derivatives, parms = parameters, rtol = 1e-6, atol = 1e-12)
  # output[, 3] <- replace(output[, 3], which(output[, 3] < 0), NA)
  output <- cbind(time = output[, 1], VL = log10(output[, 3]))
  return(data.frame(output))
}





fit_simulations <- function(simulated_parameters) {
  # This function fits the viral dynamics model for each of the simulated parameter sets.
  # 'simulated_parameters' must be a DataFrame object with simulated parameters.
  
  # Check input.
  if (is.vector(simulated_parameters)) {
    simulated_parameters <- data.frame(t(simulated_parameters))
  } else if (is.matrix(simulated_parameters)) {
    simulated_parameters <- data.frame(simulated_parameters)
  }
  
  if (!all(c("gamma", "beta", "delta", "tau") %in% colnames(simulated_parameters))) {
    stop("Data is of incorrect format. Please check input parameters.", call. = FALSE)
  } else if (nrow(simulated_parameters) == 1) {
    message("[WARNING] Only a single set of parameters found. Consider using the fit_model() function.")
  }
  
  # Initialize the predictions matrix.
  predictions <- matrix(NA,
                        nrow = length(seq(from = 0, to = max_time, by = time_interval)),
                        ncol = nrow(simulated_parameters))
  
  # Fit each simulated patient with the randomly sampled set of parameters.
  for (i in 1:nrow(simulated_parameters)) {
    parameters_i <- data.frame(gamma = simulated_parameters$gamma[i],
                               beta = simulated_parameters$beta[i],
                               delta = simulated_parameters$delta[i],
                               tau = simulated_parameters$tau[i])
    trajectory_i <- fit_model(parameters = parameters_i)
    predictions[, i] <- trajectory_i$VL
  }
  
  # Bind the time points and output the predictions matrix.
  predictions <- data.frame(cbind(seq(from = 0, to = max_time, by = time_interval), predictions))
  colnames(predictions) <- c("time", paste0("VL", 1:nrow(simulated_parameters)))
  return(predictions)
}

simulated_prediction_intervals <- function(predictions) {
  # This function determines the PIs for the simulated fits.
  # 'predictions' refers to the prediction matrix obtained from the fit_simulations() function.
  
  # Omit the time variable if present so that it is not calculated into the PIs.
  if ("time" %in% colnames(predictions)) {
    time <- predictions %>% dplyr::select("time")
    predictions <- predictions %>% dplyr::select(-c("time"))
  } else {
    time <- seq(from = 0, to = max_time, by = time_interval)
    if (length(time) != nrow(predictions)) {
      stop("Fitting of the time dimension is incorrect. Please check the prediction matrix.", call. = FALSE)
    }
  }
  
  # Define the PIs, which represent the (1 - alpha) * 100% estimates across simulations at each particular time point (row).
  lower <- apply(predictions, MARGIN = 1, FUN = function(x) (quantile(x, probs = alpha / 2, na.rm = TRUE)))
  upper <- apply(predictions, MARGIN = 1, FUN = function(x) (quantile(x, probs = 1 - (alpha / 2), na.rm = TRUE)))
  
  # Add the time variable back.
  intervals <- data.frame(cbind(time, lower, upper))
  colnames(intervals) <- c("time", "lower", "upper")
  return(intervals)
}

