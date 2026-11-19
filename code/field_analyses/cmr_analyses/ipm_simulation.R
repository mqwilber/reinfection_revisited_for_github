trajectory_simulation = function(steps, params){
  # Simulate an individual infection trajectory
  #
  # Parameters
  # ----------
  # steps : int
  #   Number of time steps to simulate the model
  # params : list
  #   List of parameters used in the IPM
  #
  # Returns
  # -------
  # : list
  #   'infection': Bd infection trajectories on the natural scale
  #   'alive': Vector of whether or not a host is alive: 1=alive, 0=dead
  
  infection_load = array(NA, dim=steps + 1) # array to store load results
  alive = array(NA, dim=steps + 1) # Array to store alive results
  infection_load[1] = 0 # Start as uninfected (however, keep in mind that we are on the log scale so we will need to convert!)
  alive[1] = 1 # Host starts as alive (1). 0 Indicates dead
  
  # Loop through the time steps
  for(i in 1:steps){
    
    load_now = infection_load[i]
    alive_now = alive[i]
    
    # Is host alive at the start of the time step?
    if(alive_now == 1){
      
      # Is host infected?
      if(load_now == 0){
        # Host is alive but uninfected
        
        # s0: An uninfected hostcould live or die
        still_alive = rbinom(1, 1, params$s0)
        
        if(still_alive == 1){
          
          alive[i + 1] = 1
          
          # phi: Does the host get infected?
          inf = rbinom(1, 1, params$phi)
          
          # G_0(x'): What is the initial load if a host gets infected?
          # Remember, we are working on the log scale so let's translate back 
          # to the natural scale with the exp.
          init = exp(rnorm(1, mean=params$init_int, sd=params$init_sigma))
          
          # Think about this: if inf = 0, you are not infected and the infection load
          # at time i + 1 = 0.  If inf = 1 then infection load at time i + 1 = init.
          infection_load[i + 1] = inf * init
          
        } else{
          
          # Uninfected host died
          infection_load[i + 1] = NA
          alive[i + 1] = 0
          
        }
        
      } else {
        # Host is alive and infected
        
        
        # s(x): Survival probability given infection
        # surv_prob = (1 / (1 + exp(-(params$surv_int + params$surv_slope*log(load_now)))))
        surv_prob = exp(-exp(params$surv_int + params$surv_slope*log(load_now + 1)))
        still_alive = rbinom(1, 1, surv_prob)
        
        if(still_alive == 1){
          
          alive[i + 1] = 1
          
          # Does the host lose infection?
          loss = rbinom(1, 1, 1 / (1 + exp(-(params$loss_int + params$loss_slope*log(load_now)))))
          
          # G(x', x): How does Bd grow?
          # Remember, this function is defined on the log scale
          new_load = exp(rnorm(1, mean=params$growth_int + params$growth_slope*log(load_now), sd=params$growth_sigma))
          
          # If loss is 1, then infection load in the next time step is 0
          # If loss is 0, the infection load in the next time step is new_load
          infection_load[i + 1] = (1 - loss) * new_load
          
        }else{
          
          # Infected host died
          infection_load[i + 1] = NA
          alive[i + 1] = 0
        }
      }
    } else {
      
      # Host dead at the start of interval. Host stays dead
      infection_load[i + 1] = NA
      alive[i + 1] = 0
    }
  }
  
  return(list(infection=infection_load, alive=alive))
  
}

set_discretized_values = function(min_size, max_size, bins){

  # Calculates the necessary parameters to use the midpoint rule to evaluate
  # the IPM model

  # Parameters
  # ----------
  # min_size : The lower bound of the integral
  # max_size : The upper bound of the integral
  # bins : The number of bins in the discretized matrix

  # Returns
  # -------
  # list
  # min_size, max_size, bins, bnd (edges of discretized kernel), y (midpoints),
  # h (width of cells)


  # Set the edges of the discretized kernel
  bnd = min_size+c(0:bins)*(max_size-min_size) / bins

  # Set the midpoints of the discretizing kernel. Using midpoint rule for evaluation
  y = 0.5 * (bnd[1:bins] + bnd[2:(bins + 1)])

  # Width of cells
  h = y[2] - y[1]

  return(list(min_size=min_size, max_size=max_size, bins=bins, bnd=bnd, y=y,
                h=h))

}

survival_function = function(x, p){
  # Survival function
  
  prob = 1 / (1 + exp(-(p$surv_int + p$surv_slope*x)))
  return(prob)
}

loss_function = function(x, p){
  # Loss of infection function
  
  prob = 1 / (1 + exp(-(p$loss_int + p$loss_slope*x)))
  return(prob)
}

transmission_function = function(p){
  # Transmission function
  
  prob = p$phi
  return(prob)
}

init_infection_function = function(x, p, h){
  # Initial infection function
  
  mu = p$init_int
  prob = pnorm(x + 0.5*h, mean=mu, sd=p$init_sigma) - pnorm(x - 0.5*h, mean=mu, sd=p$init_sigma)
  return(prob)
}

growth_function = function(xnext, xnow, p, h){
  # Growth function
  
  mu = p$growth_int + p$growth_slope * xnow
  sigma = p$growth_sigma
  prob = pnorm(xnext + 0.5*h, mu, sigma) - pnorm(xnext - 0.5*h, mu, sigma)
  return(prob)
}


build_P_matrix = function(lower, upper, bins, p){
	# Transition matrix from I to I
	#
	# Parameters
	# ----------
	# lower : lower bound on kernel
	# upper : upper bound on kernel
	# bins : number of bins in kernel
	#
	# Return
	# ------
	# : I to I transition matrix

	ipm_bnds = set_discretized_values(lower, upper, bins)
	x = ipm_bnds$y # midpoints
	h = ipm_bnds$h
  
	# Apply the mid-point rule
	G = sapply(x, function(i) growth_function(x, i, p, h))
	S = survival_function(x, p)
	L = loss_function(x, p)

	P = G %*% diag(S * (1 - L))
	return(list(P=P, G=G))
}


build_full_P_matrix = function(P, lower, upper, bins, p){
	# Build the full transition matrix with S and I classes

	ipm_bnds = set_discretized_values(lower, upper, bins)
	x = ipm_bnds$y # midpoints
	h = ipm_bnds$h

	fullP = array(0, dim=c(bins + 1, bins + 1))
  
	inf_prob = transmission_function(p)
	
	StoS = p$s0*(1 - inf_prob) # Assume trivial natural mortality over the experiment

	G0_init = init_infection_function(x, p, h)
	StoI = p$s0*inf_prob*G0_init

	S = survival_function(x, p)
	L = loss_function(x, p)
	ItoS = S*L

	# Fill in the matrix
	fullP[2:(bins + 1), 2:(bins + 1)] = P
	fullP[1, 1] = StoS
	fullP[2:(bins + 1), 1] = StoI
	fullP[1, 2:(bins + 1)] = ItoS

	return(fullP)
}

format_for_heat = function(mat, row_nm, col_nm, var_names){
  # Function converts a 2D matrix for plotting in ggplot

  tmat = mat
  colnames(tmat) = col_nm
  rownames(tmat) = row_nm
  long_mat = melt(tmat)
  colnames(long_mat) = var_names
  
  return(long_mat)

}

plot_kernel = function(P, x){
  # Visualize the matrix
  
  myPalette = colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
  long_P = format_for_heat(t(P), x, x, c('size', 'sizeNext', 'values'))
  sg_plot = ggplot(long_P, aes(x=size, y=sizeNext, fill=values)) + geom_tile() + labs(title="Survival/growth kernel")
  sg_plot = sg_plot + scale_fill_gradientn(colours = myPalette(100), name="Probability")
  sg_plot = sg_plot + theme_bw() + xlab("Log zoospore load at t") + ylab("Log zoospore load at t + 1")
  sg_plot = sg_plot + geom_abline(intercept=0, slope=1)
  sg_plot = sg_plot + theme(legend.text=element_text(size=5),
                            legend.title=element_text(size=5),
                            plot.title=element_text(size=12))
  
  return(sg_plot)
}


add_measurement_error = function(x_true, sigma_measure){
  # Add measurement error to x
  
  x_obs = rnorm(length(x_true), x_true, sigma_measure)
  return(x_obs)
}

bin_load_data = function(data, min_size, max_size, bins){
  # Take in a data.table with a column load (load on the natural scale)
  # and place them in the log load bins

  vals = set_discretized_values(min_size, max_size, bins)
  bounds = vals$bnd

  # Untransform log loads + 1
  data$logload = log(data$load)

  data[, c("state"):=-1]
  data[logload == -Inf, state:=1] # Uninfected has state 1
  data[is.na(logload), state:=(bins + 2)]

  # Put loads into their appropriate bins
  data$state[data$state == -1] = findInterval(data[state == -1]$logload, bounds) + 1

  return(data)

}

set_up_data_for_model = function(data, p_detect, bins, num_surveys=5){
  # Format the true data and generate the observed data for
  # fitting to the model

  # Calculate the true states
  unq_ids = unique(data$id)
  unq_times = unique(data$time)
  true_states = array(NA, dim=c(length(unq_ids), length(unq_times) - 1))
  load_states = array(NA, dim=c(length(unq_ids), length(unq_times) - 1))

  for(i in 1:length(unq_ids)){
    tid = unq_ids[i]
    tdat = data[id == tid]
    tobs = tdat$state
    tobs = tobs[2:length(tobs)]
    true_states[i, ] = tobs
    load_states[i, ] = tdat$logload[2:nrow(tdat)]
  }

  # Get the observed states given some detection probability
  secondary_periods = array(num_surveys, dim=c(nrow(true_states), ncol(true_states)))

  # Detection error
  secondary_obs = rbinom(nrow(true_states)*ncol(true_states), size=as.vector(secondary_periods), prob=p_detect)
  secondary_obs = array(secondary_obs, dim=dim(true_states))

  # If frog is dead, it can't be observed
  secondary_obs[true_states == bins + 2] = 0

  obs_states = true_states

  # If no frog was observed, set to unobserved state
  obs_states[secondary_obs == 0] = bins + 2

  return(list(true_states=true_states, obs_states=obs_states, 
              secondary_obs=secondary_obs, secondary_periods=secondary_periods,
              load_states=load_states))

}


load_field_data = function(site){
  # Load and format the CMR field data from a specified site.

  capture = fread(paste0("../../../data/field_data/", site, "/captures_", site, ".csv"))
  survey = fread(paste0("../../../data/field_data/", site, "/surveys_", site, "_mrmr.csv"))
  initial = fread(paste0("../../../data/field_data/", site, "/reintroductions_", site, "_mrmr.csv"))

  if(site != c("50785")){
    colnames(capture)[colnames(capture) == 'visit_date'] = 'survey_date'
  }

  if(site == c("50785")){
    # Drop three dead individuals at Rocky Basin
    capture = capture[capture_animal_state == "healthy"]
  }

  if(site %in% c("50785")){
    initial = initial[treatment != "none_occurred"]
  }

  capture = merge(capture, survey, by='survey_date', all.x=TRUE)
  capture = merge(capture, initial, by.x='pit_tag_id', by.y='pit_tag_id', all.x=TRUE)

  # Drop all individuals that were not in the original exposure
  capture = capture[pit_tag_id %in% initial$pit_tag_id]

  # How many recaptures in different primary periods?
  recaptures = capture[, .(primary_period=unique(primary_period)), by=.(pit_tag_id)]
  counts = recaptures[, .(count=.N), by=.(pit_tag_id)]
  return(list(capture=capture, survey=survey, initial=initial))

}


format_field_data_for_fitting = function(scale, bins, min_size, max_size, time_step, capture, survey, initial, site){
  # Format the CMR field data for fitting with the Stan model


  #### Step 1: Convert capture times to time steps in the model ####

  release_dates = sort(unique(capture$release_date))
  survey = survey[secondary_period != 0] # Remove all of the "release surveys"
  survey_dates = sort(unique(survey$survey_date))

  # Get the range model time steps
  step_by_date = seq(release_dates[1] + 1, survey_dates[length(survey_dates)] + time_step, by=time_step)
  step_by_date_df = data.table(date=step_by_date)
  step_by_date_df$step = 0:(length(step_by_date) - 1)

  # Put the survey dates on the model time steps
  capture$time_interval = findInterval(capture$survey_date, step_by_date_df$date)

  # Put release times on the model time steps 
  capture$release_time = findInterval(capture$release_date, step_by_date_df$date)
  initial$release_time = findInterval(initial$release_date, step_by_date_df$date)

  #### Step 2: Convert Bd load data to discrete intervals ####

  vals = set_discretized_values(min_size, max_size, bins)
  bounds = vals$bnd
  capture$logload = log(capture$bd_load / scale)
  capture[, c("state"):=-1]
  capture[logload == -Inf, state:=1]

  # Add one because these states are between 2 and bins + 1
  capture$state[capture$state == -1] = findInterval(capture[state == -1]$logload, bounds) + 1

  #### Step 3: Make a state trajectory matrix for all individuals over all time steps in the model ####

  # Individual were only swabbed once even if they were captured repeatedly
  # Drop all NA swab values
  capture_no_na = capture[!is.na(capture$state)]

  individual_ids = initial$pit_tag_id

  N = length(individual_ids)
  T = length(step_by_date) - 1

  obs_states = array(NA, dim=c(N, T))
  load_states = array(NA, dim=c(N, T))
  secondary_periods = array(-1, dim=c(N, T))
  secondary_obs = array(-1, dim=c(N, T))

  for(n in 1:N){
    
    unq_id = individual_ids[n]
    init_dat = initial[pit_tag_id == unq_id]
    release_t = init_dat$release_time[1] # When that individual was released
    
    for(t in 1:T){
      
      # If animal has been released, classify state
      if(t >= release_t){
        tdat_time = capture_no_na[time_interval == t & pit_tag_id == unq_id]
      
        if(nrow(tdat_time) > 0){
          obs_states[n, t] = tdat_time$state[1]
          load_states[n, t] = exp(tdat_time$logload[1])
        } else{
          obs_states[n, t] = bins + 2 # Undetected
          load_states[n, t] = NA
        }
      } else{
        # Animal has not yet been released so don't consider it yet
        obs_states[n, t] = -1
        load_states[n, t] = NA
      }
    }
  }

  #### Step 4:  Set-up the primary periods and secondary periods ####

  # When could the individual be seen?
  could_be_seen = survey_dates[survey_dates >= release_dates[1] + 1]
  could_be_seen_interval = findInterval(could_be_seen, step_by_date)

  # Which periods were surveyed?
  primary_period = array(0, dim=ncol(obs_states))
  primary_period[unique(could_be_seen_interval)] = 1

  # How many surveys per primary period?
  survey_trun = survey[survey_date > release_dates[1]]
  survey_trun$survey_time = findInterval(survey_trun$survey_date, step_by_date)

  # Loop through primary periods
  unq_primary_periods = unique(survey_trun$primary_period)
  secondary_periods = array(-1, dim=length(primary_period))
  secondary_obs = array(-1, dim=c(N, T))
  secondary_periods = array(-1, dim=c(N, T))
  pp_track = array(-1, dim=length(primary_period))

  for(p in 1:length(unq_primary_periods)){
    tpp = unq_primary_periods[p]
    
    pdat = survey_trun[primary_period == tpp]
    
    unq_time_steps = unique(pdat$survey_time)
    pp_track[unq_time_steps] = tpp
    
    # Maps primary and secondary periods onto model times and does not 
    # double count captures within a primary period. The trick here
    # is that some primary periods span multiple model steps and we 
    # want to account for that
    for(i in 1:length(individual_ids)){
      tid = individual_ids[i]
      
      id_dat = capture[pit_tag_id == tid & primary_period == tpp]
      tvals = obs_states[i, unq_time_steps]
      
      # Where was the individual detected if it was?
      ind = which(tvals != (bins + 2)) # 52 means undetected
      
      if(length(ind) > 0){
        
        if(all(tvals == -1)){
          # Just for consistency. Catches if individuals have not yet been released
          secondary_obs[i, unq_time_steps[ind[1]]] = nrow(id_dat)
          secondary_periods[i, unq_time_steps[ind[1]]] = nrow(pdat)
        } else{
          secondary_obs[i, unq_time_steps[ind]] = nrow(id_dat)
          secondary_periods[i, unq_time_steps[ind]] = nrow(pdat)
        }
      } else{
        secondary_obs[i, unq_time_steps[1]] = nrow(id_dat)
        secondary_periods[i, unq_time_steps[1]] = nrow(pdat)
      }

    }
  }

  #### Step 5: Clean up a couple of things in the observed states ####

  # 1. After the initial release, no animals were swabbed if they were seen again.  
  # If they were detected, assume they were uninfected so we can account for detection probability 
  # in the first interval
  # 2. The primary period seven was just a release event, so don't count it for tyndall. 
  # Same for 5 for rangersta.
  # Note that the observed states are
  # 1: Uninfected
  # 2: Infected
  # 3: Undetected

  # If we captured the frog again immediately after release, we assume it is uninfected
  obs_states[secondary_obs[, 1] > 0, 1] = 1

  # Set-up the obs_state with no load information
  obs_states_noload = obs_states
  load_states[is.na(load_states)] = -1
  obs_states_noload[(obs_states > 1) & (obs_states < bins + 2)] = 2 # Infected
  obs_states_noload[(obs_states == bins + 2)] = 3 # Undetected state
  exposure = as.numeric(initial$treatment == "exposed")


  return(list(obs_states=obs_states, 
              load_states=load_states,
              obs_states_noload=obs_states_noload,
              exposure=exposure,
              secondary_periods=secondary_periods, 
              secondary_obs=secondary_obs,
              capture=capture,
              pp_track=pp_track))

}


average_probit = function(a, beta, mu, sigma){
  return(pnorm((a + beta*mu) / (sqrt(1 + (beta*sigma)^2))))
}