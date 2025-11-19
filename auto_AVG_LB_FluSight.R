#' ES_AVERAGE 
#' This function utilizes ensembles and single automatic ARIMAX models which have mean cases by AVERAGE states as exogenous variables.
#' The function fits rolling windows of N weeks for the state under analysis and rolling windows of the same size with 1 week-lag for the exogenous variables to generate forecasts.
#' It return some metrics that evaluate the performance of the models:
#' target_end_date, abs_error, cases, forecast, 'N_of_models", weighted interval score (WIS), predictive quantiles (%)  
#' The user can choose single best automatic ARIMAXs (auto=TRUE), or ensembles of 27 permutations of 0,1,2 pdq's (ES27=TRUE) or 64 permutations of 0,1,2,3 pdq's (ES64=TRUE).
#' The user also chooses the number of weeks ahead for each forecast, and the size of the rolling window.
#' 
#' @param current_state_tb A list containing Influenza Like Illness tibbles for one or many U.S. states.
#' @param auto A logical value indicating whether to use AUTO ARIMA. Default is \code{FALSE}.
#' @param ES27 A logical value indicating whether to use ensembles of 27 models. Default is \code{TRUE}.
#' @param ES64 A logical value indicating whether to use ensembles of 64 models. Default is \code{FALSE}.
#' @param n_weeks_ahead An integer specifying the number of weeks ahead for each forecast. Default is \code{1}.
#' @param week_lag An integer specifying the week lag between the exogenous variables and the ILI cases in the state under analysis. Default is \code{1}.
#' @param list_of_states A list containing Influenza Like Illness tibbles for ALL U.S. states.
#' @param window An integer specifying the size of the rolling window in number of weeks. Default is \code{104}.
#'
#' @return A list containing the forecast results and performance metrics.

ES_AVERAGE<-function(current_state_tb, auto=FALSE, n_weeks_ahead=1, week_lag=1, ES27=TRUE, ES64=FALSE, window=104, list_of_states=list_of_states){
  
  set.seed(1)
  # The model run the ES27 or the ES64.
  ES27=!ES64 
  # Empty list that will contain forecasts, predictive quantiles and number of models.
  results<-listenv()   
  if(ES27){
    pdq=c(0,1,2) # Possible ARIMA pdq's.
    my_order_params<-permutations(3,3,pdq, repeats.allowed = TRUE) # Create 27 permutations of [0,1,2]
  }
  if(ES64){
    pdq=c(0,1,2,3) # Possible ARIMA pdq's.
    my_order_params<-permutations(4,3,pdq, repeats.allowed = TRUE) # Create 64 permutations of [0,1,2,3]
  }
  # Apply the ROLLING_ARIMA function to get results. Window set to 104 weeks (2 years).
  results[[1]]<- ROLLING_ARIMAX_AVERAGE2(current_state_tb, n_weeks_ahead=n_weeks_ahead, week_lag=week_lag, window = window, order_params=my_order_params, auto=auto, list_of_states=list_of_states) # %packages% "forecast" 

  # Put forecasts, prediction intervals and number of models into separate lists. 
  list_all_pred<-list() # Empty list for forecasts
  list_all_pred_quantiles<- list() # Empty list for predictive quantiles
  list_number_of_models<- list() # Empty list for number of models
  # Get forecasts and dates
  list_all_pred[[1]]<- results[[1]][[1]][[1]]
  # Get prediction intervals
  list_all_pred_quantiles[[1]]<-results[[1]][[2]][[1]]
  
  quantiles_by_date <- data.frame()
  for (i in 1:length(list_all_pred_quantiles[[1]])) {
    # Extract the date (name of the data frame) and the quantiles
    date <- as.Date(names(list_all_pred_quantiles[[1]][i]))
    quantiles <- t(list_all_pred_quantiles[[1]][[i]][2])
    quantiles<-expm1(quantiles)
    # Create a temporary data frame with the quantiles and date
    temp_df <- data.frame(date = date, quantiles)
    # Bind the temporary data frame to the final data frame
    quantiles_by_date <- rbind(quantiles_by_date, temp_df)
  }
  colnames(quantiles_by_date)<-c("target_end_date","0.010", "0.025", "0.050", "0.100", "0.150", "0.200", "0.250", "0.300", "0.350", "0.400", "0.450", "0.500", "0.550",
                                 "0.600", "0.650", "0.700", "0.750", "0.800", "0.850", "0.900", "0.950", "0.975", "0.990")

  results_and_quantiles<-quantiles_by_date
  return(results_and_quantiles)
}

#' ROLLING_ARIMAX_AVERAGE2
#'
#' This function works together with the ES_AVERAGE function.
#' It fits the models, generate the forecasts and counts the number of models utilized in each ensemble.
#' 
#' @param current_state_tb A tibble containing the ILI data for the current U.S. state.
#' @param n_weeks_ahead An integer specifying the number of weeks ahead for each forecast. Default is \code{1}.
#' @param window An integer specifying the number of weeks to look back for the rolling window. Default is \code{104}.
#' @param order_params A list specifying the permutation of pdq's for the ensembles. Default is \code{NULL}.
#' @param week_lag An integer specifying the week lag between the exogenous variables and the ILI cases in the state under analysis. Default is \code{1}.
#' @param auto A logical value indicating whether to use an automatic best fitted ARIMA instead of an ensembles. Default is \code{FALSE}.
#' @param list_of_states A list containing Influenza Like Illness tibbles for ALL U.S. states.
#' 
#' @return A list containing the forecast, predictive quantiles and number of models in each ensemble.
#'

ROLLING_ARIMAX_AVERAGE2 <- function(current_state_tb, week_lag=week_lag, n_weeks_ahead=1, window = 104, order_params=NULL, auto=FALSE, list_of_states=list_of_states) {
  
  set.seed(1)
  ###############################
  # MEAN CASES BY AVERAGE STATE without current state
  current_state_name <- unique(current_state_tb$State)
  print(current_state_name)
  
  # Initialize sum_cases
  sum_cases <- 0
  count <- 0
  
  # Loop through without current state
  # sum the values by the same weeks
  # calculate the mean dividing by number of states without current state
  
  # loop for summing cases
  for (state_df in list_of_states) {
    given_state_name <- unique(state_df$State)
    
    # Skip the current state
    if (given_state_name != current_state_name) {
      sum_cases <- sum_cases + state_df$cases
      count <- count + 1
    }
  }
  
  # Calculate mean
  mean_cases_AVERAGE_states <- sum_cases / count
  
  #################################
  # SOME LISTS AND VARIABLES
  # All models in the ensemble
  N_of_models<-c() 
  # Final predictions list
  prediction<-list() 
  # Predictions and dates data frame 
  prediction_df<- data.frame("predicted_date"= NULL, "Prediction" = NULL) 
  # Predictive quantiles lists
  prediction_quantile<-list() 
  prediction_quantile_ls<- list() 
  
  #################################################################
  # Iterations over the dataset adapted to the number of week_lags
  for(iter in  (1+week_lag):(NROW(current_state_tb$cases)-(window))){ 
    
    # rolling window for current state ILI data  
    current_state_rolling_window<- (iter):(window+iter-1)  
    # rolling window for the exogenous data with N week_lags
    exog_rolling_window<- (iter-week_lag):(window+iter-1-week_lag)    
    
    # list that will get our ARIMA models
    fitted_models<-list()
    # list that will get the AIC scores 
    model_aic_scores<-c() 
    # Model id, utilized in the loop
    model_id<-1
    
    ##################################################
    # AVERAGE states time series for each iteration #
    # Selection of the calculated mean cases based on the exogenous window, which has 1 week lag
    AVERAGE_states_dataset= data.frame(mean_cases_AVERAGE_states[exog_rolling_window])
    #########################################
    # Exogenous variable for each iteration #
    # mean cases in the date the forecasting is being made
    exog_var<-c((mean_cases_AVERAGE_states[104+(iter-1)]), 
                (mean_cases_AVERAGE_states[104+(iter-1)]), 
                (mean_cases_AVERAGE_states[104+(iter-1)]), 
                (mean_cases_AVERAGE_states[104+(iter-1)]))
    
    ##########
    # Fitting 
    ##########
    
    # Start if we have 104 elements in the rolling window.
    if(length(current_state_tb$cases[current_state_rolling_window])==window){ 
      
      # run 1 time for the auto.arima or run 27 or 64 times for the ensembles
      for(j in 1:nrow(order_params)){
        fit<- NULL # start with fit as NULL
        # try to fit an ARIMA model
        tryCatch(
          expr = {
            if(!auto){
              # if auto = FALSE, run the ensembles of 27 or 64
              set.seed(1)
              # fit ensembles of ARIMAs on log1p of the data
              fit <- Arima(
                ts(log1p(current_state_tb$cases[current_state_rolling_window]), frequency = 52),
                xreg = log1p(AVERAGE_states_dataset[, 1]),
                order = order_params[j,],
                seasonal = list(order = c(1, 0, 1), period = 52),
                method = "CSS-ML")
            }
            
            # if auto = TRUE, run auto.arima
            else{
              set.seed(1)
              fit<-invisible(auto.arima(log1p(current_state_tb$cases[current_state_rolling_window]), xreg=log1p(AVERAGE_states_dataset[,1]) ,stepwise=TRUE))
              }
            
            # save each fitted ARIMA in fitted_models[[j]]
            fitted_models[[j]]<-fit
            # save the AIC of each fitted model
            model_aic_scores[model_id]<- fit$aic
          }
          # fit will be NULL if there is an error in the fitting process
          ,error = function(e){
          }
        )#end tryCatch
        
        # If fit == NULL, save the fitted model and the AIC as NAN
        if(is.null(fit) || is.null(fitted_models[[j]])){ 
          fitted_models[[j]]<-NA
          model_aic_scores[model_id]<- NA
        }
        # if auto==TRUE break the model on the first run
        # since we just need one result and not 27 or 64
        if(auto)
          break
        # model_ids are important for the ensembles
        # for the auto.arima it will be == 1
        model_id<-model_id+1 
      }
      
      ##############
      # Forecasting 
      ##############
      
      # general initial variables
      predicted_value<- 0 # predicted values 
      m<- numeric(n_weeks_ahead) # mean forecast value
      s<- numeric(n_weeks_ahead) # standard deviation
      sims<-c() # simulations for the mixture of gaussians  
      
      # Ensemble weights initial variables
      model_weights<- c()
      min_aic<- min(model_aic_scores, na.rm = TRUE) # min models' aic
      total_aic<-sum(exp(-.5*(model_aic_scores-min_aic)), na.rm =TRUE ) # sum of aics without nan values
      
      # Counts the number of models utilized in each forecast  
      my_n_models<-0 
      for(my_model in fitted_models){ # for each valid model in fitted_models sum 1 to my_n_models
        if(length(my_model)>0 && !is.na(my_model[1]) && !(is.na(my_model$aic))){
          my_n_models<-my_n_models+1 
        }}
      
      ######################################################
      # Save the number of models utilized in each iteration  
      N_of_models<-append(N_of_models,my_n_models)
      
      ###############################
      # Generate the target_end_dates
      # Generates a sequence of dates based on the last date by n weeks ahead
      weekly_dates<- current_state_tb$target_end_date[current_state_rolling_window] # current data inside the 104 weeks window
      last_date <- max(weekly_dates) # last date of this window
      my_predicted_dates <- seq.Date(from = last_date + 7 , by = "week", length.out = n_weeks_ahead) 
      predicted_date<-my_predicted_dates[n_weeks_ahead]
      
      ##########################################
      # run the models stored on fitted models #
      for(my_model in fitted_models){
        ######################
        # if a model is valid
        if(length(my_model)>0 && !is.na(my_model[1]) && !(is.na(my_model$aic))){ 
          # calculate the weights based on the total AIC previously calculated
          model_weights_<- exp(-.5*(my_model$aic - min_aic))/total_aic
          
          set.seed(1)
          #######################################################################################
          # simulate values for each model based on its weights to build the mixture of Gaussians
          new.sims<-c() # list of simulations for each mixture of Gaussians                  
          fc <- forecast(my_model, h=n_weeks_ahead, xreg=log1p(exog_var[1:n_weeks_ahead]), level=99)
          m <- fc$mean[n_weeks_ahead]  ## mean forecast value 
          s <- ((fc$upper[n_weeks_ahead]-fc$lower[n_weeks_ahead])/2.576/2)  # 99% confidence interval
          n<-ceiling(model_weights_*1e6) # number of simulations based on the weights
          new.sims <- rnorm(n, m=m, sd=s) # simulate values for each weighted model in as a gaussian
          sims <- c(sims, new.sims) # combine simulated values for each model
          
          #####################################
          # calculate the ensemble prediction #
          predicted_value <- model_weights_*m + predicted_value ### m = mean forecast for 99% confidence
          ###########################
          # if predicted value is NA 
          if(is.na(predicted_value)){
            print("predicted_value is na")
          }
        }
      }
      
      # get the prediction dataset
      single_prediction<- data.frame("predicted_date"= predicted_date, "Prediction" = predicted_value)# get the forecast for that date
      #rbind all forecasts and dates
      prediction_df<-rbind(prediction_df, single_prediction) 
      # Define the 23 quantiles
      probabilities <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99) 
      # get 23 predictive quantiles based on the mixture of gaussians 
      my_quantiles <- quantile(sims, probs=probabilities)
      
      ########################
      # Predictive quantiles #
      ########################
      # Creating predictive quantiles levels index
      pi_level<-c(0.01, 0.025, 0.05, 0.10, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
                  0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)
      
      # reset the initial prediction_df_quantile for each iteration
      prediction_df_quantile<- data.frame("pi_level"= NULL, "quantile"= NULL, "point_forecast" = NULL)# empity dataframe of prediction intervals
      # save the 23 predictive quantiles, the upper and lower bounds for the given week ahead, and the ensemble forecasts into a data frame
      for(j in 1:23){ 
        single_df_quantile<- data.frame("pi_level"= pi_level[j],"quantile"= my_quantiles[j], "point_forecast" = predicted_value)# fill prediction intervals dataframe with correct values for each week ahead
        prediction_df_quantile<-rbind(prediction_df_quantile, single_df_quantile)
      }
      # put it into a list named with its forecast target_end_date
      prediction_quantile_ls[[toString(predicted_date)]]<-prediction_df_quantile 
    }
    # If don't have 104 weeks of data
    else
      print(paste0("Not enough values"))
  }
  # after everything is done
  print("Complete.")
  
  prediction[[1]]<-prediction_df # save the list of forecasts
  prediction_quantile[[1]]<- prediction_quantile_ls # save the list of predictive quantiles
  df_N_of_models<-data.frame(N_of_models) # save the number of models utilized in each forecast
  
  return(list("Point_ForeCast "=prediction, "Quantiles"=prediction_quantile, "Number_of_models"=df_N_of_models))
}


